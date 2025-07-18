from banzai_nres.fibers import fiber_states_from_header
from banzai.lco import LCOFrameFactory, LCOObservationFrame, LCOCalibrationFrame
from banzai.frames import ObservationFrame
from banzai.data import DataProduct, ArrayData, HeaderOnly, DataTable
from banzai.logs import get_logger
from banzai.utils import fits_utils
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot
from io import BytesIO
from typing import Optional
import numpy as np
from astropy.table import Table
import os
from astropy.coordinates import Angle
from astropy import units


logger = get_logger()


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
    def mask(self):
        return self._table['mask']

    @mask.setter
    def mask(self, value):
        self._table['mask'] = value

    @property
    def table(self):
        # return the full spectrum, unmasked.
        return self._table

    @property
    def fibers_and_orders(self):
        return self._table['fiber'], self._table['order']


class NRESObservationFrame(LCOObservationFrame):
    def __init__(self, hdu_list: list, file_path: str, frame_id: int = None, hdu_order: list = None):
        LCOObservationFrame.__init__(self, hdu_list, file_path, frame_id=frame_id, hdu_order=hdu_order)
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.meta)
        self.classification = None

        self._traces = None
        self._profile = None
        self._blaze = None
        self._weights = None
        self._fibers = None
        self._features = None
        self._wavelengths = None
        self._spectrum = None
        self._ccf = None
        self.summary_figures = []

    def get_output_data_products(self, runtime_context):
        if self.obstype != 'TARGET':
            return super().get_output_data_products(runtime_context)
        else:
            filename_1d = self.get_output_filename(runtime_context)
            filename_1d = filename_1d.replace('.fits', '-1d.fits')
            filename_2d = filename_1d.replace('-1d.fits', '-2d.fits')
            filename_summary = filename_1d.replace('-1d.fits', '-summary.pdf').replace('.fz', '')

            frame_1d_meta = self.meta.copy()
            # Remove any name from the header of the primary hdu
            if 'EXTNAME' in frame_1d_meta:
                frame_1d_meta.pop('EXTNAME')

            if 'CCF' in self:
                frame_1d = LCOObservationFrame([HeaderOnly(frame_1d_meta, None), self['SPECTRUM1D'], self['CCF']],
                                               os.path.join(self.get_output_directory(runtime_context), filename_1d))
            else:
                frame_1d = LCOObservationFrame([HeaderOnly(frame_1d_meta, None), self['SPECTRUM1D']],
                                               os.path.join(self.get_output_directory(runtime_context), filename_1d))
            fits_1d = frame_1d.to_fits(runtime_context)
            fits_1d['SPECTRUM1D'].name = 'SPECTRUM'
            fits_1d[0].header['L1ID2D'] = filename_2d
            output_product_1d = DataProduct.from_fits(fits_1d, filename_1d, self.get_output_directory(runtime_context))

            frame_2d = LCOObservationFrame([hdu for hdu in self._hdus if hdu.name not in ['SPECTRUM1D', 'CCF']],
                                           os.path.join(self.get_output_directory(runtime_context), filename_2d))
            fits_2d = frame_2d.to_fits(runtime_context)
            fits_2d[0].header['L1ID1D'] = filename_1d
            output_product_2d = DataProduct.from_fits(fits_2d, filename_2d, self.get_output_directory(runtime_context))

            summary_pdf_meta = fits_utils.sanitize_header(self.meta.copy())
            for keyword in ['NAXIS1', 'NAXIS2', 'EXTNAME', 'EXTVER']:
                if keyword in summary_pdf_meta:
                    summary_pdf_meta.pop(keyword)

            summary_buffer = BytesIO()
            pp = PdfPages(summary_buffer, keep_empty=False)
            for fig in self.summary_figures:
                fig.tight_layout()
                pp.savefig(fig)
                pyplot.close(fig)
            pp.close()
            summary_buffer.seek(0)

            output_summary = DataProduct(summary_buffer, filename_summary, self.get_output_directory(runtime_context),
                                         meta=dict(summary_pdf_meta))

            return [output_product_1d, output_product_2d, output_summary]

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
        if 'TRACES' in self:
            return self['TRACES'].data
        else:
            return self._traces

    @traces.setter
    def traces(self, value):
        self._traces = value

    @property
    def background(self):
        return self['BACKGROUND'].data

    @background.setter
    def background(self, value):
        self.add_or_update(ArrayData(value, name='BACKGROUND'))

    @property
    def profile(self):
        if 'PROFILE' in self:
            return self['PROFILE'].data
        else:
            return self._profile

    @profile.setter
    def profile(self, value):
        self._profile = value

    @property
    def weights(self):
        if 'WEIGHTS' in self:
            return self['WEIGHTS'].data
        else:
            return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = value

    @property
    def spectrum(self):
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):
        self._spectrum = value
        self.add_or_update(DataTable(value.table, name='SPECTRUM1D'))

    @property
    def blaze(self):
        if 'BLAZE' in self:
            return self['BLAZE'].data
        else:
            return self._blaze

    @blaze.setter
    def blaze(self, value):
        self._blaze = value

    @property
    def wavelengths(self):
        if 'WAVELENGTH' in self:
            return self['WAVELENGTH'].data
        else:
            return self._wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        self._wavelengths = value

    @property
    def fibers(self):
        if 'FIBERS' in self:
            return self['FIBERS'].data
        else:
            return self._fibers

    @fibers.setter
    def fibers(self, value):
        self._fibers = value

    @property
    def features(self):
        if 'FEATURES' in self:
            return self['FEATURES'].data
        else:
            return self._features

    @features.setter
    def features(self, value):
        self._features = value

    @property
    def ccf(self):
        if 'CCF' in self:
            return self['CCF'].data
        else:
            return self._ccf

    @ccf.setter
    def ccf(self, value):
        self._ccf = value

    @property
    def ra(self):
        if self.primary_hdu.meta['RA'] == 'N/A':
            return np.nan
        return Angle(self.primary_hdu.meta['RA'], units.hourangle).deg

    @ra.setter
    def ra(self, value):
        if isinstance(value, float) and np.isnan(value):
            self.primary_hdu.meta['RA'] = 'N/A'
        else:
            a = Angle(value, units.deg)
            self.primary_hdu.meta['RA'] = a.to_string(unit=units.hourangle, sep=':', precision=3, pad=True)

    @property
    def dec(self):
        if self.primary_hdu.meta['DEC'] == 'N/A':
            return np.nan
        return Angle(self.primary_hdu.meta['DEC'], units.deg).deg

    @dec.setter
    def dec(self, value):
        if isinstance(value, float) and np.isnan(value):
            self.primary_hdu.meta['DEC'] = 'N/A'
        else:
            a = Angle(value, units.deg)
            self.primary_hdu.meta['DEC'] = a.to_string(unit=units.deg, sep=':', precision=2, pad=True)

    @property
    def pm_ra(self):
        if self.primary_hdu.meta['PM-RA'] == 'N/A':
            return np.nan
        # Proper motion is stored in arcseconds/year but we always use it in mas/year
        # Note that the RA proper motion has the cos dec term included both in the header and when we use it
        return self.primary_hdu.meta['PM-RA'] * 1000.0

    @pm_ra.setter
    def pm_ra(self, value):
        if isinstance(value, str):
            self.primary_hdu.meta['PM-RA'] = value
        elif isinstance(value, float) and np.isnan(value):
            self.primary_hdu.meta['PM-RA'] = 'N/A'
        else:
            # Proper motion is stored in arcseconds/year but we always use it in mas/year
            # Note that the RA proper motion has the cos dec term included both in the header and when we use it
            self.primary_hdu.meta['PM-RA'] = value / 1000.0
        self.primary_hdu.meta.comments['PM-RA'] = 'RA proper motion from Gaia [mas/yr * cos(Dec)]'

    @property
    def pm_dec(self):
        if self.primary_hdu.meta['PM-DEC'] == 'N/A':
            return np.nan
        # Proper motion is stored in arcseconds/year but we always use it in mas/year
        return self.primary_hdu.meta['PM-DEC'] * 1000.0

    @pm_dec.setter
    def pm_dec(self, value):
        if isinstance(value, str):
            self.primary_hdu.meta['PM-DEC'] = value
        elif isinstance(value, float) and np.isnan(value):
            self.primary_hdu.meta['PM-DEC'] = 'N/A'
        else:
            # Proper motion is stored in arcseconds/year but we always use it in mas/year
            self.primary_hdu.meta['PM-DEC'] = value / 1000.0
        self.primary_hdu.meta.comments['PM-DEC'] = 'Dec proper motion from Gaia [mas/yr]'

    @property
    def classification(self):
        return self._classification

    @classification.setter
    def classification(self, value):
        self._classification = value
        if value is not None:
            self.meta['TEFF'] = value.T_effective, 'Estimated stellar effective temperature [K]'
            self.meta['LOGG'] = value.log_g, 'Estimated stellar surface gravity [cgs]'
            self.meta['FEH'] = value.metallicity, 'Estimated stellar metallicity [dex]'
            self.meta['ALPHA'] = value.alpha, 'Estimated stellar alpha abundance [dex]'
        else:
            self.meta['TEFF'] = ''
            self.meta['LOGG'] = ''
            self.meta['FEH'] = ''
            self.meta['ALPHA'] = ''

    @property
    def num_traces(self):
        """
        Counts the number of illuminated orders on the detector by taking
        the largest label present in the trace image.
        :return: int
        """
        return int(np.max(self.traces))


class NRESCalibrationFrame(LCOCalibrationFrame, NRESObservationFrame):
    def __init__(self, hdu_list: list, file_path: str, frame_id: int = None, grouping_criteria: list = None,
                 hdu_order: list = None):
        LCOCalibrationFrame.__init__(self, hdu_list, file_path,  grouping_criteria=grouping_criteria)
        NRESObservationFrame.__init__(self, hdu_list, file_path, frame_id=frame_id, hdu_order=hdu_order)

    def write(self, runtime_context):
        LCOCalibrationFrame.write(self, runtime_context)


class NRESFrameFactory(LCOFrameFactory):

    @property
    def observation_frame_class(self):
        return NRESObservationFrame

    @property
    def calibration_frame_class(self):
        return NRESCalibrationFrame

    @staticmethod
    def is_empty_coordinate(coordinate):
        return any(coordinate.lower() in value for value in ['nan', 'n/a', 'unknown'])

    def open(self, path, runtime_context) -> Optional[ObservationFrame]:
        image = super().open(path, runtime_context)
        # Currently we can't distinguish between the NRES composite data products and the raw frames off the telescope
        # As such we have to check an an extra header keyword. We put nres01 etc in TELESCOP in the composite data
        # products to distinguish.
        if image is None or 'nres' not in image.meta.get('TELESCOP', '').lower():
            return None

        # Not all NRES CDPs have all the extensions like bias frames and BPMs so we have to check
        if 'TELESCOPE_1' in image:

            # Fix the RA and DEC keywords to be the requested values for all of our measurements
            objects = image.meta['OBJECTS'].lower().split('&')
            if objects[0] != 'none':
                telescope_num = 1
            else:
                telescope_num = 2

            if self.is_empty_coordinate(image[f'TELESCOPE_{telescope_num}'].meta['CAT-RA']):
                ra_dec_keyword = ''
            else:
                ra_dec_keyword = 'CAT-'

            if self.is_empty_coordinate(image[f'TELESCOPE_{telescope_num}'].meta[f'{ra_dec_keyword}RA']):
                image.ra = np.nan
                image.dec = np.nan
            else:
                image.ra = Angle(image[f'TELESCOPE_{telescope_num}'].meta[f'{ra_dec_keyword}RA'], units.hourangle).deg
                image.dec = image[f'TELESCOPE_{telescope_num}'].meta[f'{ra_dec_keyword}DEC']
            if image[f'TELESCOPE_{telescope_num}'].meta['PM-RA'] == 'N/A':
                image.pm_ra = np.nan
                image.pm_dec = np.nan
            else:
                # Convert to mas / yr from arcsec / year
                image.pm_ra = image[f'TELESCOPE_{telescope_num}'].meta['PM-RA'] * 1000.0
                image.pm_dec = image[f'TELESCOPE_{telescope_num}'].meta['PM-DEC'] * 1000.0
            image.meta['RADESYS'] = image[f'TELESCOPE_{telescope_num}'].meta['RADESYS']
        return image
