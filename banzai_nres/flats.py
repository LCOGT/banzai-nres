from banzai.calibrations import CalibrationStacker


class FlatStacker(CalibrationStacker):
    def __init__(self, runtime_context):
        super(FlatStacker, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'LAMPFLAT'
