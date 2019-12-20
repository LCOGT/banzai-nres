from banzai.stages import Stage
from scipy.interpolate import LSQBivariateSpline
import numpy as np


class BackgroundSubtractor(Stage):
    def do_stage(self, image):
        # TODO: Add iteratively reweighted LSQ
        indices_to_interpolate = np.logical_or(image.trace > 0, image.mask)
        x1d, y1d = np.arange(image.shape[1]), np.arange(image.shape[0])
        X, Y = np.meshgrid(x1d, y1d)
        weights = image.primary_hdu.uncertainty ** -2.0
        estimator = LSQBivariateSpline(X[~indices_to_interpolate], Y[~indices_to_interpolate],
                                       image.data[~indices_to_interpolate], w=weights[~indices_to_interpolate],
                                       tx=np.arange(self.runtime_context.background_smoothing_scale, image.shape[1], 
                                                    self.runtime_context.background_smoothing_scale),
                                       ty=np.arange(self.runtime_context.background_smoothing_scale, image.shape[0], 
                                                    self.runtime_context.background_smoothing_scale),
                                                    kx=self.runtime_context.background_spline_order, 
                                                    ky=self.runtime_context.background_spline_order)
                                       
        image.background = estimator(x1d, y1d).T
        return image
