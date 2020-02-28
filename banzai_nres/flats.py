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
        if not image.is_master():
            master_calibration_image.primary_hdu.name = 'LAMPFLAT'
            image.append(master_calibration_image.primary_hdu)
        image.traces = master_calibration_image.traces
        # TODO the below is a place holder until we load the profile from the flat.
        image.profile = np.ones_like(image.data)


class FlatDivider(Stage):
    pass
