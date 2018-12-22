from banzai.calibrations import CalibrationStacker


class FlatStacker(CalibrationStacker):
    def __init__(self, pipeline_context):
        super(FlatStacker, self).__init__(pipeline_context)

    @property
    def calibration_type(self):
        return 'LAMPFLAT'
