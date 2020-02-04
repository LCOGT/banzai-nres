from banzai.stages import Stage
import numpy as np
from banzai_nres.fitting import fit_smooth_spline


class BackgroundSubtractor(Stage):
    def do_stage(self, image):
        indices_to_interpolate = np.logical_or(image.traces > 0, image.mask)
        estimator = fit_smooth_spline(image.data, error=image.error, mask=indices_to_interpolate)
        image.background = estimator(np.arange(image.shape[0], dtype=np.float),
                                     np.arange(image.shape[1], dtype=np.float))
        return image
