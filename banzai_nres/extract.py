import numpy as np
from astropy.table import Table
import abc

from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.utils.extract_utils import index_traces


class SpectrumExtractor(Stage):
    def do_stage(self, image: EchelleSpectralCCDData):
        flux = np.zeros((image.num_traces, image.traces.shape[1]), dtype=float)
        variance = np.zeros_like(flux)
        trace_ids, indexed_traces = index_traces(image.traces)

        for i, itrace in enumerate(indexed_traces):
            weights = self.weights(itrace, image)
            flux[i] = self.extract_order(image.data[itrace], weights)
            variance[i] = self.extract_order(image.uncertainty[itrace] ** 2, weights ** 2)

        image.spectrum = Table({'flux': flux, 'uncertainty': np.sqrt(variance), 'id': trace_ids})
        return image

    @staticmethod
    def extract_order(values, weights):
        return np.sum(values * weights, axis=0)

    @abc.abstractmethod
    def weights(self, itrace, image):
        """
        :param itrace: ndarray. ndarray of indices such that image.data[itrace] returns
        an ndarray of flux of shape N,M where M is something like 4096 pixels (num of x pixels in the image data)
        and N is the width of that trace, e.g. 20.
        :param image:
        :return: ndarray. Has same shape as itrace.
        """
        return self._weights(itrace, image.profile, image.uncertainty**2, image.mask)

    def _weights(self, itrace, profile_im, var_im, mask):
        """
        Calculates the optimal extraction weights for a single order
        per Horne (1986, 1986PASP...98..609H) Equation 8 .
        The optimally extracted spectrum is then the sum of data * weights
        The associated variance is by definition var * weights**2 , which agrees with Horne (1986)
        equation 9.
        :param itrace: ndarray. ndarray of indices such that image.data[itrace] returns
        an ndarray of flux of shape N,M where M is something like 4096 pixels (num of x pixels in the image data)
        and N is the width of that trace, e.g. 20.
        :param profile_im: ndarray. 2d image of the profile.
        :param var_im: ndarray. 2d image of the variance per pixel
        :param mask: ndarray. 2d image where 0 is a good pixel and 1 is a masked, bad pixel.
        :return: ndarray. Has same shape as itrace.
        """
        if np.allclose(mask, 1):
            return np.zeros_like(itrace)
        invmask = np.invert(mask[itrace].astype(bool))  # array such that pixels to ignore have value 0
        return invmask * profile_im[itrace] / var_im[itrace] / self.extract_order(1, invmask *
                                                                                  profile_im[itrace] ** 2 / var_im[itrace])



class BoxExtractor(SpectrumExtractor):
    def weights(self, trace, image):
        return np.ones_like(trace)
