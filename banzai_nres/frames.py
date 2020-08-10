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
import os


logger = logging.getLogger('banzai')


class Spectrum1D:
    def __init__(self, data):
        self._table = Table(data)
        if 'mask' not in self._table.colnames:
            self._table['mask'] = np.zeros_like(self._table['wavelength'], dtype=np.uint8)
        self._table['mask'] |= np.array([row['wavelength'] == 0.0 for row in self._table], dtype=np.uint8)

    def __getitem__(self, item):
        """
        Get the spectrum given a fiber and order
        :param item: tuple of fiber, order
        :return: dict
        """
        fiber, order = item
        row = self._table[np.logical_and(self._table['fiber'] == fiber, self._table['order'] == order)][0]
        good_pixels = row['mask'] == 0
        return {column: row[column][good_pixels] if isinstance(row[column], np.ndarray) else row[column]
                for column in row.colnames if column != 'mask'}

    def __setitem__(self, key, value):
        fiber, order, column_name = key
        if column_name not in self._table.colnames:
            self._table.add_column([np.zeros_like(row['flux']) for row in self._table], name=column_name)
        correct_row = np.where(np.logical_and(self._table['fiber'] == fiber, self._table['order'] == order))[0][0]
        good_pixels = self._table[correct_row]['mask'] == 0
        # call column name first then the logical bool array for the astropy table to actually set the values
        self._table[column_name][correct_row, good_pixels] = value


    @property
    def fibers_and_orders(self):
        return self._table['fiber'], self._table['order']

    def to_fits(self, extname):
        return fits.BinTableHDU(self._table, name=extname, header=fits.Header({'EXTNAME': extname}))


class NRESObservationFrame(LCOObservationFrame):
    def __init__(self, hdu_list: list, file_path: str, frame_id: int = None, hdu_order: list = None):
        LCOObservationFrame.__init__(self, hdu_list, file_path, frame_id=frame_id, hdu_order=hdu_order)
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.meta)

    def num_lit_fibers(self):
        return 1 * self.fiber0_lit + 1 * self.fiber1_lit + 1 * self.fiber2_lit

    def get_output_directory(self, runtime_context) -> str:
        return os.path.join(runtime_context.processed_path, self.instrument.site,
                            self.instrument.name, self.epoch, 'processed')

    @property
    def science_fiber(self):
        if self.fiber0_lit:
            return 0
        elif self.fiber2_lit:
            return 2

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

    @property
    def fibers(self):
        return self.primary_hdu.fibers

    @fibers.setter
    def fibers(self, value):
        self.primary_hdu.fibers = value

    @property
    def ccf(self):
        return self.primary_hdu.ccf

    @ccf.setter
    def ccf(self, value):
        self.primary_hdu.ccf = value


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
                 spectrum: Table = None, blaze: Table = None, memmap=True, features: Table = None,
                 fibers: np.array = None, ccf: Table = None):
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
        self.fibers = fibers
        self.ccf = ccf

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

    @property
    def fibers(self):
        return self._fibers

    @fibers.setter
    def fibers(self, value):
        self._fibers = value

    @property
    def ccf(self):
        return self._ccf

    @ccf.setter
    def ccf(self, value):
        self._ccf = value

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
            hdu_list.append(self.spectrum.to_fits(extname))
        if self.blaze is not None:
            extname = 'BLAZE'
            hdu_list.append(fits.BinTableHDU(self.blaze, name=extname, header=fits.Header({'EXTNAME': extname})))
        if self.features is not None:
            extname = 'FEATURES'
            hdu_list.append(fits.BinTableHDU(self.features, name=extname, header=fits.Header({'EXTNAME': extname})))

        if self.fibers is not None:
            extname = 'FIBERS'
            hdu_list.append(fits.BinTableHDU(self.fibers, name=extname, header=fits.Header({'EXTNAME': extname})))

        if self.ccf is not None:
            extname = 'CCF'
            hdu_list.append(fits.BinTableHDU(self.ccf, name=extname, header=fits.Header({'EXTNAME': extname})))
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
                                                          {'FITS_NAME': 'WAVELENGTH', 'NAME': 'wavelengths'},
                                                          {'FITS_NAME': 'FIBERS', 'NAME': 'fibers'},
                                                          {'FITS_NAME': 'CCF', 'NAME': 'ccf'}]

    def open(self, path, runtime_context) -> Optional[ObservationFrame]:
        image = super().open(path, runtime_context)
        # Currently we can't distinguish between the NRES composite data products and the raw frames off the telescope
        # As such we have to check an an extra header keyword. We put nres01 etc in TELESCOP in the composite data
        # products to distinguish.
        if image is None or 'nres' not in image.meta.get('TELESCOP', '').lower():
            return None
        else:
            return image
