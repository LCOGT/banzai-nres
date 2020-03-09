from banzai_nres.fibers import fiber_states_from_header
from banzai.utils.fits_utils import to_fits_image_extension
from banzai.lco import LCOFrameFactory, LCOObservationFrame, LCOMasterCalibrationFrame, LCOCalibrationFrame
from banzai.frames import ObservationFrame
from banzai.data import CCDData
import logging
from typing import Optional
import numpy as np
from astropy.table import Table
from typing import Union
from astropy.io import fits

logger = logging.getLogger('banzai')


class NRESObservationFrame(LCOObservationFrame):
    def __init__(self, hdu_list: list, file_path: str):
        LCOObservationFrame.__init__(self, hdu_list, file_path)
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.meta)

    def num_lit_fibers(self):
        return 1 * self.fiber0_lit + 1 * self.fiber1_lit + 1 * self.fiber2_lit

    @property
    def traces(self):
        return self.primary_hdu.traces

    @traces.setter
    def traces(self, value):
        self.primary_hdu.traces = value

    @property
    def num_traces(self):
        return self.primary_hdu.num_traces

    @property
    def background(self):
        return self.primary_hdu.background

    @background.setter
    def background(self, value):
        self.primary_hdu.background = value

    @property
    def profile(self):
        return self.primary_hdu.profile

    @profile.setter
    def profile(self, value):
        self.primary_hdu.profile = value

    @property
    def weights(self):
        return self.primary_hdu.weights

    @weights.setter
    def weights(self, value):
        self.primary_hdu.weights = value

    @property
    def spectrum(self):
        return self.primary_hdu.spectrum

    @spectrum.setter
    def spectrum(self, value):
        self.primary_hdu.spectrum = value


class NRESCalibrationFrame(LCOCalibrationFrame, NRESObservationFrame):
    def __init__(self, hdu_list: list, file_path: str, grouping_criteria: list = None):
        LCOCalibrationFrame.__init__(self, hdu_list, file_path, grouping_criteria)
        NRESObservationFrame.__init__(self, hdu_list, file_path)


class NRESMasterCalibrationFrame(LCOMasterCalibrationFrame, NRESCalibrationFrame):
    def __init__(self, images: list, file_path: str, grouping_criteria: list = None):
        NRESCalibrationFrame.__init__(self, images, file_path, grouping_criteria)
        LCOMasterCalibrationFrame.__init__(self, images, file_path, grouping_criteria)


class EchelleSpectralCCDData(CCDData):
    def __init__(self, data: Union[np.array, Table], meta: fits.Header,
                 mask: np.array = None, name: str = '', uncertainty: np.array = None,
                 background: np.array = None,  traces: np.array = None,
                 profile: np.array = None, weights: np.array = None,
                 spectrum: Table = None, memmap=True):
        super().__init__(data=data, meta=meta, mask=mask, name=name, memmap=memmap, uncertainty=uncertainty)
        if traces is None:
            self._traces = None
        else:
            self.traces = traces
        if background is None:
            self._background = None
        else:
            self.background = background
        if profile is None:
            self._profile = None
        else:
            self.profile = profile
        if weights is None:
            self._weights = None
        else:
            self.weights = weights
        if spectrum is None:
            self._spectrum = None
        else:
            self.spectrum = spectrum

    @property
    def traces(self):
        return self._traces

    @traces.setter
    def traces(self, value):
        self._traces = self._init_array(value)

    @property
    def num_traces(self):
        """
        Counts the number of illuminated orders on the detector by taking
        the largest label present in the trace image.
        :return: int
        """
        return int(np.max(self.traces))

    @property
    def profile(self):
        return self._profile

    @profile.setter
    def profile(self, value):
        self._profile = self._init_array(value)

    @property
    def spectrum(self):
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):
        self._spectrum = value

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = self._init_array(value)

    @property
    def background(self):
        return self._background

    @background.setter
    def background(self, value):
        self.data -= value
        self._background = self._init_array(value)

    def to_fits(self, context):
        hdu_list = super().to_fits(context)
        if self.traces is not None:
            hdu_list.append(to_fits_image_extension(self.traces, self.extension_name, 'TRACES', context,
                                                    extension_version=self.meta.get('EXTVER')))
        if self.background is not None:
            hdu_list.append(to_fits_image_extension(self.background, self.extension_name, 'BACKGROUND', context,
                                                    extension_version=self.meta.get('EXTVER')))
        if self.weights is not None:
            hdu_list.append(to_fits_image_extension(self.weights, self.extension_name, 'BACKGROUND', context,
                                                    extension_version=self.meta.get('EXTVER')))
        if self.spectrum is not None:
            hdu_list.append(fits.BinTableHDU(self.spectrum, name=self.extension_name,
                                             header=fits.Header({'EXTNAME': self.extension_name + '1DSPEC'})))
        return hdu_list


class NRESFrameFactory(LCOFrameFactory):
    observation_frame_class = NRESObservationFrame
    calibration_frame_class = NRESCalibrationFrame
    data_class = EchelleSpectralCCDData
    associated_extensions = LCOFrameFactory().associated_extensions + [{'FITS_NAME': 'TRACES', 'NAME': 'traces'},
                                                                       {'FITS_NAME': 'BACKGROUND', 'NAME': 'background'}]

    def open(self, path, runtime_context) -> Optional[ObservationFrame]:
        image = super().open(path, runtime_context)
        # Currently we can't distinguish between the NRES composite data products and the raw frames off the telescope
        # As such we have to check an an extra header keyword. We put nres01 etc in TELESCOP in the composite data
        # products to distinguish.
        if image is None or 'nres' not in image.meta.get('TELESCOP', '').lower():
            return None
        else:
            return image
