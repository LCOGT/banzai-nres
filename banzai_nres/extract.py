import numpy as np

from banzai.stages import Stage
from banzai_nres.frames import Spectrum1D, NRESObservationFrame
from banzai_nres.utils.trace_utils import get_trace_region
import logging
from banzai.data import ArrayData


logger = logging.getLogger('banzai')


class WeightedExtract(Stage):
    def do_stage(self, image: NRESObservationFrame):
        if image.weights is None:
            logger.error('Extraction weights are missing. Rejecting image.', image=image)
            return None
        if image.wavelengths is None:
            logger.error('Wavelengths are missing. Rejecting image.', image=image)
            return None
        # consider adding a method to image, i.e. image.extracted_spectrum_shape.
        flux = np.zeros((image.num_traces, image.data.shape[1]), dtype=float)
        wavelength = np.zeros_like(flux, dtype=float)
        variance = np.zeros_like(flux, dtype=float)
        mask = np.zeros_like(flux, dtype=np.uint8)

        trace_ids = np.arange(1, image.num_traces + 1)
        for i, trace_id in enumerate(trace_ids):
            this_trace = get_trace_region(np.isclose(image.traces, trace_id))
            # get the horizontal (x) extent of the trace. Consider making this a get_extent function.
            x_extent = slice(np.min(this_trace[1]), np.max(this_trace[1]) + 1)
            flux[i, x_extent] = self.extract_order(image.data[this_trace], image.weights[this_trace])
            variance[i, x_extent] = self.extract_order(image.uncertainty[this_trace] ** 2,
                                                       image.weights[this_trace] ** 2)
            # get the average wavelength: Sum wavelengths weighted by 1 over the vertical width of the trace (e.g. 1/10)
            wavelength[i, x_extent] = self.extract_order(image.wavelengths[this_trace],
                                                         weights=1/image.wavelengths[this_trace].shape[0])
            mask[i, x_extent] = image.weights[this_trace].sum(axis=0) == 0.0

        image.spectrum = Spectrum1D({'id': trace_ids, 'order': image.fibers['order'],
                                     'fiber': image.fibers['fiber'], 'wavelength': wavelength,
                                     'flux': flux, 'uncertainty': np.sqrt(variance), 'blaze': image.blaze['blaze'],
                                     'blaze_error': image.blaze['blaze_error'], 'mask': mask})
        return image

    @staticmethod
    def extract_order(values, weights=1):
        """
        Performs a weighted sum of values vertically (i.e. along the rows (along axis 0)), weighted by weights.
        :param values: ndarray, shape N, M. The values to extract by summing along the rows.
        :param weights: float, or ndarray with the same shape as values: N, M.
        weights are the extraction weights to be applied to each value in values.
        If weights is float, every value in values will have that weight.
        :return: out, ndarray.
                 out has shape M. Meaning values has been summed vertically (across the rows).
        """
        return np.sum(values * weights, axis=0)


class GetOptimalExtractionWeights(WeightedExtract):
    def do_stage(self, image: NRESObservationFrame):
        if image.profile is None:
            logger.error('Profile missing. Rejecting image.', image=image)
            return None
        image['WEIGHTS'] = ArrayData(np.zeros_like(image.data, dtype=float), name='WEIGHTS')
        trace_ids = np.arange(1, image.num_traces + 1)
        for trace_id in trace_ids:
            yx = get_trace_region(np.isclose(image.traces, trace_id))
            image.weights[yx] = self.weights(image.profile[yx], image.uncertainty[yx]**2, image.mask[yx])
        return image

    def weights(self, profile_im, var_im, mask):
        """
        Calculates the optimal extraction weights for a single order
        per Horne (1986, 1986PASP...98..609H) Equation 8 .
        The optimally extracted spectrum is then the sum of data * weights
        The associated variance is by definition var * weights**2 , which agrees with Horne (1986)
        equation 9.

        Note: Our conventions differ slightly with that of Horne 1986.
        The weights by which the spectrum is extracted in Equation 8 of Horne are:
        .. math::
            W_{x, \\lambda} = \frac{M_{x, \\lambda}  P_{x, \\lambda} / V_{x, \\lambda},
            \\sum_x M_{x, \\lambda} P^2_{x, \\lambda} / V_{x, \\lambda} }
        These quantities are labelled here as:
        .. math:: profile_im = P_{x, \\lambda}
        .. math:: var_im = V_{x, \\lambda}
        .. math:: mask = M_{x, \\lambda}

        Note that in Horne, x is the cross-dispersed direction, which for us is the
        vertical (y) pixels.

        :param profile_im: ndarray. 2d image of the profile.
        :param var_im: ndarray. 2d image of the variance per pixel
        :param mask: ndarray. 2d image where 0 is a good pixel and 1 is a masked, bad pixel.
        :return: ndarray. Has same shape as profile_im, var_im and mask.
        """
        invmask = np.invert(mask.astype(bool))  # array such that pixels to ignore have value 0
        normalization = self.extract_order(invmask * profile_im ** 2 / var_im)
        weights = np.divide(invmask * profile_im / var_im, normalization)
        weights[~np.isfinite(weights)] = 0  # remove infinities and nans from division by zero.
        return weights


class BoxExtractionWeights(GetOptimalExtractionWeights):
    def weights(self, profile_im, var_im, mask):
        return 1.0
