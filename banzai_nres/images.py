import abc
from astropy.table import Table
from astropy import fits

from banzai.images import Image
from banzai_nres.fibers import fiber_states_from_header
import banzai_nres.settings as nres_settings  # import to override banzai settings
from banzai import settings
from banzai_nres.utils import fits_utils, db_utils
from banzai import dbs


class NRESImage(Image):
    def __init__(self, runtime_context, filename=None, data=None, data_tables=None,
                 header=None, extension_headers=None, bpm=None):
        super(NRESImage, self).__init__(runtime_context, filename=filename, data=data, data_tables=data_tables,
                                        header=header, extension_headers=extension_headers, bpm=bpm)
        self.trace = None
        self.rectified_2d_spectrum = None
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.header)

    def num_lit_fibers(self):
        return 1 * self.fiber0_lit + 1 * self.fiber1_lit + 1 * self.fiber2_lit


class TableImage(object):
    """
    Image like object with .write and .load methods. For use in Trace and Blaze calibrations
    """
    def __init__(self, data=None, table_name=None, num_centers_per_trace=0, filepath=None,
                 header=None, image=None, obstype='TRACE'):
        if data is None and num_centers_per_trace <= 0:
            raise ValueError('Trace object instantiated but no trace data given and num_centers_per_trace is not > 0')
        if data is None:
            data = Table([Column(name='id'), Column(name='centers', shape=(num_centers_per_trace,))])
            data['id'].description = 'Identification tag for trace'
            data['centers'].description = 'Vertical position of the center of the' \
                                          ' trace as a function of horizontal pixel'
            data['centers'].unit = 'pixel'
        if header is None:
            header = {}
        self.header = header
        self.filepath = filepath
        self.data = Table(data)
        self.table_name = table_name

    def write(self, runtime_context=None, update_db=True):
        hdu = fits.BinTableHDU(self.data, name=self.table_name, header=fits.Header(self.header))
        hdu_list = fits.HDUList([fits.PrimaryHDU(), hdu])
        self._update_filepath(runtime_context)
        fits_utils.writeto(hdu_list=hdu_list, filepath=self.filepath,
                           fpack=getattr(runtime_context, 'fpack', False),
                           overwrite=True, output_verify='fix+warn')
        if update_db:
            dbs.save_calibration_info(self.filepath, image=self,
                                      db_address=runtime_context.db_address)
            if runtime_context.post_to_archive:
                db_utils.post_to_archive(self.filepath)

    def _update_filepath(self, runtime_context):
        if getattr(runtime_context, 'fpack', False) and not self.filepath.endswith('.fz'):
            self.filepath += '.fz'

    @staticmethod
    @abc.abstractmethod
    def load(path, table_name):
        pass
