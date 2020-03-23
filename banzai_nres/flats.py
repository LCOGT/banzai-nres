from banzai.calibrations import CalibrationStacker, CalibrationUser
from banzai.stages import Stage
import numpy as np


class FlatStacker(CalibrationStacker):
    def __init__(self, runtime_context):
        super(FlatStacker, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'LAMPFLAT'


class FlatLoader(CalibrationUser):
    @property
    def calibration_type(self):
        return 'LAMPFLAT'

    def on_missing_master_calibration(self, image):
        if image.obstype.upper() == 'LAMPFLAT':
            return image
        else:
            super().on_missing_master_calibration(image)

    def apply_master_calibration(self, image, master_calibration_image):
        image.traces = master_calibration_image.traces
        image.profile = master_calibration_image.profile
        image.blaze = master_calibration_image.blaze
        return image
