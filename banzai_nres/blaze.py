"""
blaze.py: Driver script for making blaze calibration files for echelle spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from astropy.table import Table
from scipy import ndimage
import numpy as np
import os

from banzai.calibrations import Stage, create_master_calibration_header
from banzai.utils import file_utils
from banzai_nres.images import ImageBase
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
        lampflat_spectrum = Table(image.data_tables[nres_settings.BOX_SPECTRUM_EXTNAME]._data_table)
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
