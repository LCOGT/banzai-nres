from banzai import Stage
import numpy as np


class InitializeInverseVariances(Stage):
    """
    This initializes the inverse variance matrix (per pixel 1/Var) for the image. If a pixel x,y, has N counts,
    then the associated variance is simply N + read_noise^2
    """
    def __init__(self, pipeline_context):
        super(InitializeInverseVariances, self).__init__(pipeline_context)

    def do_stage(self, images):
        for image in images:
            if not hasattr(image, 'ivar'):
                setattr(image, 'ivar', None)
            read_var = 100 #image.header['RDNOISE']**2
            shot_var = image.data.astype(np.float32)
            vars = shot_var + read_var
            vars[vars < read_var] = read_var
            image.ivar = np.reciprocal(vars)