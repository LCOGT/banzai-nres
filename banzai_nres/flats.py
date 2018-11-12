from banzai.calibrations import CalibrationMaker


class FlatStacker(CalibrationMaker):
    def __init__(self, pipeline_context):
        super(CalibrationMaker, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum', 'fibers_state']

    @property
    def calibration_type(self):
        return 'LAMPFLAT'

    @property
    def min_images(self):
        return 5

    def make_master_calibration_frame(self, images, image_config):
        pass
