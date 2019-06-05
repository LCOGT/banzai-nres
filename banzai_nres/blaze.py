"""
blaze.py: Driver script for making blaze calibration files for echelle spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from astropy.table import Table
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
        lampflat_spectrum = image.data_tables[nres_settings.BOX_SPECTRUM_EXTNAME]
        blaze = Blaze(data=None, filepath=master_filepath, header=master_header, image=image,
                      table_name=self.blaze_table_name,
                      obstype=self.calibration_type)
        logger.info('Created master blaze file', image=image, extra_tags={'calibration_type': self.calibration_type,
                                                                          'output_path': master_filepath,
                                                                          'calibration_obstype': master_header['OBSTYPE']})
        return blaze

    @staticmethod
    def _get_filepath(runtime_context, lampflat_image, master_filename):
        output_directory = file_utils.make_output_directory(runtime_context, lampflat_image)
        return os.path.join(output_directory, os.path.basename(master_filename))


class Blaze(ImageBase):
    """
    :param data = {'id': ndarray, 'centers': ndarray}. 'centers' gives a 2d array, where
    the jth row are the y centers across the detector for the trace with identification trace_centers['id'][j]
    """
    def __init__(self, data=None, table_name=None, filepath=None,
                 header=None, image=None, obstype='BLAZE'):
        super(Blaze, self).__init__(data=data, table_name=table_name, filepath=filepath,
                                    header=header, image=image, obstype=obstype)
        self.data = Table(data)
