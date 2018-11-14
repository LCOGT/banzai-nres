from banzai.calibrations import CalibrationStacker


class FlatStacker(CalibrationStacker):
    def __init__(self, pipeline_context):
        super(FlatStacker, self).__init__(pipeline_context)

    @property
    def group_by_attributes(self):
        return ['ccdsum', 'fibers_state']

    @property
    def calibration_type(self):
        return 'LAMPFLAT'

    @property
    def min_images(self):
        return 5
