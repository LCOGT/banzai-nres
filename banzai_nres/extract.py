import numpy as np
from astropy.table import Table
import abc

from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.utils.extract_utils import get_trace_yx_positions


class SpectrumExtractor(Stage):
    def do_stage(self, image: EchelleSpectralCCDData):
        flux = np.zeros((image.num_traces, image.traces.shape[1]), dtype=float)
        variance = np.zeros_like(flux)
        trace_ids, trace_yxpos = get_trace_yx_positions(image.traces)

        for i, yx in enumerate(trace_yxpos):
            x_extent = slice(np.min(yx[1]), np.max(yx[1]) + 1)  # get the horizontal (x) extent of the trace.
            weights = self.weights(yx, image)
            flux[i, x_extent] = self.extract_order(image.data[yx], weights)
            variance[i, x_extent] = self.extract_order(image.uncertainty[yx] ** 2, weights ** 2)

        image.spectrum = Table({'flux': flux, 'uncertainty': np.sqrt(variance), 'id': trace_ids})
        return image

    @staticmethod
    def extract_order(values, weights=1):
        """
        :param values: ndarray. values to extract by summing along the rows.
        :param weights: float, or ndarray with same shape as values.
        weights are the extraction weights to be applied to each value in values.
        If weights is float, every value in values will have that weight.
        :return:
        """
        return np.sum(values * weights, axis=0)

    def weights(self, yx, image):
        """
        :param yx: list of tuples or ndarrays suitable for numpy's fancy indexing.
        yx needs to be such that image.data[yx] returns
        an ndarray of flux of shape N,M where M is e.g. 4096 (num of x pixels in the image data)
        and N is the vertical height of that trace, e.g. 20.
        :param image: EchelleSpectralCCDData.
        :return: float or ndarray. See SpectrumExtractor.extract_order
        """
        return self._weights(yx, image.profile, image.uncertainty**2, image.mask)

    def _weights(self, yx, profile_im, var_im, mask):
        """
        Calculates the optimal extraction weights for a single order
        per Horne (1986, 1986PASP...98..609H) Equation 8 .
        The optimally extracted spectrum is then the sum of data * weights
        The associated variance is by definition var * weights**2 , which agrees with Horne (1986)
        equation 9.
        :param yx: list of tuples or ndarrays suitable for numpy's fancy indexing.
        yx needs to be such that image.data[yx] returns
        an ndarray of flux of shape N,M where M is e.g. 4096 (num of x pixels in the image data)
        and N is the vertical height of that trace, e.g. 20.
        :param profile_im: ndarray. 2d image of the profile.
        :param var_im: ndarray. 2d image of the variance per pixel
        :param mask: ndarray. 2d image where 0 is a good pixel and 1 is a masked, bad pixel.
        :return: ndarray. Has same shape as yx.
        """
        invmask = np.invert(mask[yx].astype(bool))  # array such that pixels to ignore have value 0
        normalization = self.extract_order(invmask * profile_im[yx] ** 2 / var_im[yx])
        weights = np.divide(invmask * profile_im[yx] / var_im[yx], normalization)
        weights[~np.isfinite(weights)] = 0  # remove infinities and nans from division by zero.
        return weights


class BoxExtractor(SpectrumExtractor):
    def weights(self, yx, image):
        return 1
