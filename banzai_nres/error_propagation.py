from banzai.stages import Stage
import numpy as np


class InitializeInverseVariances(Stage):
    """
    This initializes the inverse variance matrix (per pixel 1/Var) for the image. If a pixel x,y, has N counts,
    then the associated variance is simply N + read_noise^2
    """
    def __init__(self, pipeline_context):
        super(InitializeInverseVariances, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'filling_variances'

    def do_stage(self, images):
        for image in images:
            if not hasattr(image, 'ivar'):
                setattr(image, 'ivar', None)
            read_var = 100 #image.header['RDNOISE']**2
            # the image is converted to electrons at this point, is the read noise already converted as well?
            shot_var = image.data.astype(np.float32)
            vars = shot_var + read_var
            del shot_var
            vars[vars < read_var] = read_var
            image.ivar = np.reciprocal(vars)
            del vars
        return images