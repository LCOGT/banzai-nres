"""
blaze.py: Driver script for making blaze calibration files for echelle spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from astropy.table import Table
from scipy import ndimage
import numpy as np
import os

from banzai.calibrations import Stage, create_master_calibration_header, ApplyCalibration
from banzai.utils import file_utils
from banzai_nres.images import ImageBase
from banzai.images import DataTable
import banzai_nres.settings as nres_settings
import banzai.settings as banzai_settings

import logging

logger = logging.getLogger(__name__)


class BlazeMaker(Stage):
    def __init__(self, runtime_context):
        super(BlazeMaker, self).__init__(runtime_context)
        self.runtime_context = runtime_context
        self.blaze_table_name = nres_settings.BLAZE_TABLE_NAME

    @property
    def calibration_type(self):
        return 'BLAZE'

    def do_stage(self, image):
        master_header = create_master_calibration_header(old_header=image.header, images=[image])
        master_header['OBSTYPE'] = self.calibration_type
        master_filename = banzai_settings.CALIBRATION_FILENAME_FUNCTIONS[self.calibration_type](image)
        master_filepath = self._get_filepath(self.runtime_context, image, master_filename)
        logger.info('Making blaze file', image=image)
        lampflat_spectrum = as_astropy_table(image.data_tables[nres_settings.BOX_SPECTRUM_EXTNAME])
        # TODO fix the following issue: the .write method only works if ImageBase.data is an AstropyTable,
        # not a banzai.images.DataTable. Make DataTable in banzai inherit from astropy Table?
        lampflat_spectrum['flux'] = ndimage.median_filter(lampflat_spectrum['flux'].data, (1, 101))
        lampflat_spectrum = self._clip_spectrum(lampflat_spectrum, minval=10 * 100)  # 10 times the read noise.
        lampflat_spectrum = self._normalize_max_value(lampflat_spectrum)
        blaze = ImageBase(data=lampflat_spectrum, filepath=master_filepath, header=master_header,
                          image=image, table_name=self.blaze_table_name, obstype=self.calibration_type)
        logger.info('Created master blaze file', image=image, extra_tags={'calibration_type': self.calibration_type,
                                                                          'output_path': master_filepath,
                                                                          'calibration_obstype': master_header['OBSTYPE']})
        return blaze

    @staticmethod
    def _get_filepath(runtime_context, lampflat_image, master_filename):
        output_directory = file_utils.make_output_directory(runtime_context, lampflat_image)
        return os.path.join(output_directory, os.path.basename(master_filename))

    @staticmethod
    def _clip_spectrum(spec, minval):
        spec['flux'] = spec['flux'].data * (spec['flux'].data > minval) + minval * (spec['flux'].data <= minval)
        return spec

    @staticmethod
    def _normalize_max_value(spec):
        spec['flux'] /= np.max(spec['flux'].data, axis=1).reshape(-1, 1)
        return spec


class ApplyBlaze(ApplyCalibration):
    """
    Loads the blaze spectrum from file and divides any arc fibers by the blaze.
    """
    def __init__(self, runtime_context):
        super(ApplyBlaze, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'BLAZE'

    def apply_master_calibration(self, image, master_calibration_path):
        blaze = ImageBase.load(master_calibration_path, extension_name=nres_settings.TRACE_TABLE_NAME)
        master_filename = os.path.basename(master_calibration_path)
        image.header['L1IDBLAZ'] = (master_filename, 'ID of blaze file')
        logger.info('Loaded blaze file', image=image,  extra_tags={'L1IDBLAZ': image.header['L1IDBLAZ']})
        spectrum = image.data_tables[nres_settings.BOX_SPECTRUM_EXTNAME]
        if not np.allclose(spectrum['id'].data, blaze['id'].data):
            logger.error('Trace IDs from blaze spectrum and spectrum do not agree. Aborting '
                         'blaze correction.', image=image)
        else:
            # TODO only divide the lit_arc_lamps.
            spectrum['flux'] = spectrum['flux'].data / blaze['flux'].data
            image.data_tables[nres_settings.BOX_SPECTRUM_EXTNAME] = spectrum
        return image

    def do_stage(self, image):
        master_calibration_path = self.get_calibration_filename(image)
        if master_calibration_path is None:
            self.on_missing_master_calibration(image)
            return image
        return self.apply_master_calibration(image, master_calibration_path)


def as_astropy_table(table):
    #TODO make DataTable behave with fits.BinTableHDU, so that we do not need this munge function.
    if isinstance(table, DataTable):
        return Table(table._data_table)
    else:
        return Table(table)
