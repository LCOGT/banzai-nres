import numpy as np
from astropy.table import Table
import abc

from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.utils.extract_utils import get_trace_xy_positions


class SpectrumExtractor(Stage):
    def do_stage(self, image: EchelleSpectralCCDData):
        flux = np.zeros((image.num_traces, image.traces.shape[1]), dtype=float)
        variance = np.zeros_like(flux)
        trace_ids, trace_xypos = get_trace_xy_positions(image.traces)

        for i, xy in enumerate(trace_xypos):
            weights = self.weights(xy, image)
            flux[i] = self.extract_order(image.data[xy], weights)
            variance[i] = self.extract_order(image.uncertainty[xy] ** 2, weights ** 2)

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

    @abc.abstractmethod
    def weights(self, xy, image):
        """
        :param xy: ndarray. ndarray of indices such that image.data[xy] returns
        an ndarray of flux of shape N,M where M is e.g. 4096 (num of x pixels in the image data)
        and N is the vertical height of that trace, e.g. 20.
        :param image: EchelleSpectralCCDData.
        :return: float or ndarray. See SpectrumExtractor.extract_order
        """
        return self._weights(xy, image.profile, image.uncertainty**2, image.mask)

    def _weights(self, xy, profile_im, var_im, mask):
        """
        Calculates the optimal extraction weights for a single order
        per Horne (1986, 1986PASP...98..609H) Equation 8 .
        The optimally extracted spectrum is then the sum of data * weights
        The associated variance is by definition var * weights**2 , which agrees with Horne (1986)
        equation 9.
        :param xy: ndarray. ndarray of indices such that image.data[xy] returns
        an ndarray of flux of shape N,M where M is e.g. 4096 (num of x pixels in the image data)
        and N is the vertical height of that trace, e.g. 20.
        :param profile_im: ndarray. 2d image of the profile.
        :param var_im: ndarray. 2d image of the variance per pixel
        :param mask: ndarray. 2d image where 0 is a good pixel and 1 is a masked, bad pixel.
        :return: ndarray. Has same shape as xy.
        """
        invmask = np.invert(mask[xy].astype(bool))  # array such that pixels to ignore have value 0
        normalization = self.extract_order(invmask * profile_im[xy] ** 2 / var_im[xy])
        return np.divide(invmask * profile_im[xy] / var_im[xy], normalization,
                         out=np.zeros_like(xy), where=~np.isclose(normalization, 0))


class BoxExtractor(SpectrumExtractor):
    def weights(self, xy, image):
        return 1
