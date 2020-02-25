from banzai.stages import Stage
import numpy as np
import sep


class BackgroundSubtractor(Stage):
    def do_stage(self, image):
        indices_to_interpolate = np.logical_or(image.traces > 0, image.mask)
        image.background = sep.Background(image.data, mask=indices_to_interpolate).back()
        return image
