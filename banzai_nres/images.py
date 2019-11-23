from banzai_nres.fibers import fiber_states_from_header
from banzai.images import LCOFrameFactory, LCOObservationFrame, LCOMasterCalibrationFrame, LCOCalibrationFrame
import logging

logger = logging.getLogger('banzai')


class NRESFrame:
    def __init__(self, header):
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(header)

    def num_lit_fibers(self):
        return 1 * self.fiber0_lit + 1 * self.fiber1_lit + 1 * self.fiber2_lit


class NRESObservationFrame(LCOObservationFrame, NRESFrame):
    def __init__(self, hdu_list: list, file_path: str):
        LCOObservationFrame.__init__(self, hdu_list, file_path)
        NRESFrame.__init__(self, self.meta)


class NRESCalibrationFrame(LCOCalibrationFrame, NRESFrame):
    def __init__(self, hdu_list: list, file_path: str, grouping_criteria: list = None):
        LCOCalibrationFrame.__init__(self, hdu_list, file_path, grouping_criteria)
        NRESFrame.__init__(self, self.meta)


class NRESMasterCalibrationFrame(LCOMasterCalibrationFrame, NRESFrame):
    def __init__(self, images: list, file_path: str, grouping_criteria: list = None):
        LCOMasterCalibrationFrame.__init__(self, images, file_path, grouping_criteria)
        NRESFrame.__init__(self, self.meta)


class NRESFrameFactory(LCOFrameFactory):
    observation_frame_class = NRESObservationFrame
    calibration_frame_class = NRESCalibrationFrame
