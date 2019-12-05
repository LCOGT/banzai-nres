from banzai_nres.fibers import fiber_states_from_header
from banzai.images import LCOFrameFactory
from banzai.images import ObservationFrame, LCOObservationFrame, LCOMasterCalibrationFrame, LCOCalibrationFrame
import logging
from typing import Optional

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

    @classmethod
    def open(cls, path, runtime_context) -> Optional[ObservationFrame]:
        image = super(cls).open(path, runtime_context)
        # Currently we can't distinguish between the NRES composite data products and the raw frames off the telescope
        # As such we have to check an an extra header keyword. We put nres01 etc in TELESCOP in the composite data
        # products to distinguish.
        if image is None or 'nres' not in image.meta.get('TELESCOP', '').lower():
            return None
        else:
            return image
