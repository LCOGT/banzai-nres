from banzai_nres.fibers import fiber_states_from_header
from banzai.utils.fits_utils import to_fits_image_extension
from banzai.images import LCOFrameFactory
from banzai.images import ObservationFrame, LCOObservationFrame, LCOMasterCalibrationFrame, LCOCalibrationFrame, CCDData
import logging
from typing import Optional
import numpy as np
from astropy.table import Table
from typing import Union
from astropy.io import fits

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
        self.trace = None

    @property
    def traces(self):
        return self.primary_hdu.traces

    @traces.setter
    def traces(self, value):
        self.primary_hdu.traces = value

    @property
    def background(self):
        return self.primary_hdu.background

    @background.setter
    def background(self, value):
        self.primary_hdu.background = value


class NRESCalibrationFrame(LCOCalibrationFrame, NRESFrame):
    def __init__(self, hdu_list: list, file_path: str, grouping_criteria: list = None):
        LCOCalibrationFrame.__init__(self, hdu_list, file_path, grouping_criteria)
        NRESFrame.__init__(self, self.meta)


class NRESMasterCalibrationFrame(LCOMasterCalibrationFrame, NRESFrame):
    def __init__(self, images: list, file_path: str, grouping_criteria: list = None):
        LCOMasterCalibrationFrame.__init__(self, images, file_path, grouping_criteria)
        NRESFrame.__init__(self, self.meta)


class EchelleSpectralCCDData(CCDData):
    def __init__(self, data: Union[np.array, Table], meta: fits.Header,
                 mask: np.array = None, name: str = '', uncertainty: np.array = None,
                 background: np.array = None,  traces: np.array = None, memmap=True):
        super().__init__(data=data, meta=meta, mask=mask, name=name, memmap=memmap, uncertainty=uncertainty)
        if traces is None:
            self._traces = None
        else:
            self.traces = traces
        self.background = background

    @property
    def traces(self):
        return self._traces

    @traces.setter
    def traces(self, value):
        self._traces = self._init_array(value)

    @property
    def background(self):
        return self._background

    @background.setter
    def background(self, value):
        self._background = self._init_array(value)

    def to_fits(self, context):
        hdu_list = super().to_fits(context)
        hdu_list.append(to_fits_image_extension(self.traces, self.extension_name, 'TRACES', context,
                                                extension_version=self.meta.get('EXTVER')))
        hdu_list.append(to_fits_image_extension(self.background, self.extension_name, 'BACKGROUND', context,
                                                extension_version=self.meta.get('EXTVER')))
        return hdu_list


class NRESFrameFactory(LCOFrameFactory):
    observation_frame_class = NRESObservationFrame
    calibration_frame_class = NRESCalibrationFrame
    data_class = EchelleSpectralCCDData
    associated_extensions = LCOFrameFactory().associated_extensions + [{'FITS_NAME': 'TRACES', 'NAME': 'traces'},
                                                                       {'FITS_NAME': 'BACKGROUND', 'NAME': 'background'}]

    @classmethod
    def open(cls, path, runtime_context) -> Optional[ObservationFrame]:
        image = super().open(path, runtime_context)
        # Currently we can't distinguish between the NRES composite data products and the raw frames off the telescope
        # As such we have to check an an extra header keyword. We put nres01 etc in TELESCOP in the composite data
        # products to distinguish.
        if image is None or 'nres' not in image.meta.get('TELESCOP', '').lower():
            return None
        else:
            return image
