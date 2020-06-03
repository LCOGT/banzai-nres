from banzai_nres.fibers import fiber_states_from_header
from banzai.utils.fits_utils import to_fits_image_extension
from banzai.lco import LCOFrameFactory, LCOObservationFrame, LCOCalibrationFrame
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
    def __init__(self, hdu_list: list, file_path: str, frame_id: int = None, hdu_order: list = None):
        LCOObservationFrame.__init__(self, hdu_list, file_path, frame_id=frame_id, hdu_order=hdu_order)
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

    @property
    def blaze(self):
        return self.primary_hdu.blaze

    @blaze.setter
    def blaze(self, value):
        self.primary_hdu.blaze = value

    @property
    def features(self):
        return self.primary_hdu.features

    @features.setter
    def features(self, value):
        self.primary_hdu.features = value

    @property
    def line_list(self):
        return self.primary_hdu.line_list

    @line_list.setter
    def line_list(self, value):
        self.primary_hdu.line_list = value

    @property
    def wavelengths(self):
        return self.primary_hdu.wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        self.primary_hdu.wavelengths = value


class NRESCalibrationFrame(LCOCalibrationFrame, NRESObservationFrame):
    def __init__(self, hdu_list: list, file_path: str, frame_id: int = None, grouping_criteria: list = None,
                 hdu_order: list = None):
        LCOCalibrationFrame.__init__(self, hdu_list, file_path,  grouping_criteria=grouping_criteria)
        NRESObservationFrame.__init__(self, hdu_list, file_path, frame_id=frame_id, hdu_order=hdu_order)


class EchelleSpectralCCDData(CCDData):
    def __init__(self, data: Union[np.array, Table], meta: fits.Header,
                 mask: np.array = None, name: str = '', uncertainty: np.array = None,
                 background: np.array = None,  traces: np.array = None, wavelengths: np.array = None,
                 profile: np.array = None, weights: np.array = None, line_list=None,
                 spectrum: Table = None, blaze: Table = None, memmap=True, features: Table = None):
        super().__init__(data=data, meta=meta, mask=mask, name=name, memmap=memmap, uncertainty=uncertainty)
        if traces is None:
            self._traces = None
        else:
            self.traces = traces
        if wavelengths is None:
            self._wavelengths = None
        else:
            self.wavelengths = wavelengths
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

        self.spectrum = spectrum
        self.blaze = blaze
        self.features = features
        self.line_list = line_list


    @property
    def traces(self):
        return self._traces

    @traces.setter
    def traces(self, value):
        self._traces = self._init_array(value)

    @property
    def traces(self):
        return self._traces

    @traces.setter
    def traces(self, value):
        self._traces = self._init_array(value)

    @property
    def wavelengths(self):
        return self._wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        # Initialize wavelengths to zeros. See note in banzai_nres.wavelength.WavelengthCalibrate() for explanation.
        self._wavelengths = self._init_array(value)

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
        if self.profile is not None:
            hdu_list.append(to_fits_image_extension(self.profile, self.extension_name, 'PROFILE', context,
                                                    extension_version=self.meta.get('EXTVER')))
        if self.weights is not None:
            hdu_list.append(to_fits_image_extension(self.weights, self.extension_name, 'WEIGHTS', context,
                                                    extension_version=self.meta.get('EXTVER')))
        if self.wavelengths is not None:
            hdu_list.append(to_fits_image_extension(self.wavelengths, self.extension_name, 'WAVELENGTH', context,
                                                    extension_version=self.meta.get('EXTVER')))
        if self.spectrum is not None:
            extname = '1DSPEC'
            hdu_list.append(fits.BinTableHDU(self.spectrum, name=extname, header=fits.Header({'EXTNAME': extname})))
        if self.blaze is not None:
            extname = 'BLAZE'
            hdu_list.append(fits.BinTableHDU(self.blaze, name=extname, header=fits.Header({'EXTNAME': extname})))
        if self.features is not None:
            extname = 'FEATURES'
            hdu_list.append(fits.BinTableHDU(self.features, name=extname, header=fits.Header({'EXTNAME': extname})))

        return hdu_list


class NRESFrameFactory(LCOFrameFactory):

    @property
    def observation_frame_class(self):
        return NRESObservationFrame

    @property
    def calibration_frame_class(self):
        return NRESCalibrationFrame

    @property
    def data_class(self):
        return EchelleSpectralCCDData

    @property
    def associated_extensions(self):
        return LCOFrameFactory().associated_extensions + [{'FITS_NAME': 'TRACES', 'NAME': 'traces'},
                                                          {'FITS_NAME': 'BACKGROUND', 'NAME': 'background'},
                                                          {'FITS_NAME': 'PROFILE', 'NAME': 'profile'},
                                                          {'FITS_NAME': 'BLAZE', 'NAME': 'blaze'},
                                                          {'FITS_NAME': 'WAVELENGTH', 'NAME': 'wavelengths'}]

    def open(self, path, runtime_context) -> Optional[ObservationFrame]:
        image = super().open(path, runtime_context)
        # Currently we can't distinguish between the NRES composite data products and the raw frames off the telescope
        # As such we have to check an an extra header keyword. We put nres01 etc in TELESCOP in the composite data
        # products to distinguish.
        if image is None or 'nres' not in image.meta.get('TELESCOP', '').lower():
            return None
        else:
            return image
