from banzai.stages import Stage
import numpy as np
from scipy.interpolate import SmoothBivariateSpline


class BackgroundSubtractor(Stage):
    def do_stage(self, image):

        indices_to_interpolate = np.logical_or(image.traces > 0, image.mask)

        x1d, y1d = np.arange(image.shape[1]), np.arange(image.shape[0])
        X, Y = np.meshgrid(x1d, y1d)

        # Note that from the docs, the weights should be about 1 over the standard deviation for the smoothing to
        # work well
        # TODO: add reweighting
        # because we can control the weights we can add iterative reweighting easily to do a better job of
        # of rejecting outliers
        weights = 1.0 / image.primary_hdu.uncertainty
        estimator = SmoothBivariateSpline(X[~indices_to_interpolate], Y[~indices_to_interpolate],
                                          image.data[~indices_to_interpolate], w=weights[~indices_to_interpolate])
        image.background = estimator(x1d, y1d).T
        return image
