from banzai.calibrations import CalibrationStacker, CalibrationUser
from banzai.images import MasterCalibrationFrame


class FlatStacker(CalibrationStacker):
    def __init__(self, runtime_context):
        super(FlatStacker, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'LAMPFLAT'


class FlatLoader(CalibrationUser):
    def calibration_type(self):
        return 'LAMPFLAT'

    def apply_master_calibration(self, image, master_calibration_image):
        if not isinstance(image, MasterCalibrationFrame):
            master_calibration_image.primary_hdu.name = 'LAMPFLAT'
            image.append(master_calibration_image.primary_hdu)
        image.trace = master_calibration_image.trace
