import numpy as np
from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from astropy.table import Table
import abc


class SpectrumExtractor(Stage):
    def do_stage(self, image: EchelleSpectralCCDData):
        flux = np.zeros((67, 4000), dtype=float)
        variance = np.zeros((67, 4000), dtype=float)

        for i, trace in enumerate(image.traces):
            weights = self.weights(trace, image)
            flux[i] = self.extract_order(image.data[trace], weights)
            variance[i] = self.extract_order(image.uncertainty[trace] ** 2, weights ** 2)

        image.spectrum = Table({'flux': flux, 'uncertainty': np.sqrt(variance)})
        return image

    @staticmethod
    def extract_order(values, weights):
        return np.sum(values * weights, axis=0)

    @abc.abstractmethod
    def weights(self, trace, image):
        """
        :param trace: ndarray
        :param image:
        :return: ndarray. Has same shape as trace.
        """
        return self._weights(trace, image.profile, image.uncertainty**2, image.mask)

    def _weights(self, trace, profile_im, var_im, mask):
        """
        Calculates the optimal extraction weights for a single order
        per Horne (1986, 1986PASP...98..609H) Equation 8 .
        The optimally extracted spectrum is then the sum of data * weights
        The associated variance is by definition var * weights**2 , which agrees with Horne (1986)
        equation 9.
        :param trace: ndarray
        :param profile_im: ndarray. 2d image of the profile.
        :param var_im: ndarray. 2d image of the variance per pixel
        :param mask: ndarray. 2d image where 0 is a good pixel and 1 is a masked, bad pixel.
        :return: ndarray. Has same shape as trace.
        """
        # TODO calculate the weights.


class BoxExtractor(SpectrumExtractor):
    def weights(self, trace, image):
        return np.ones_like(trace)
