"""
blaze.py: Driver script for making blaze calibration files for echelle spectrographs.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from astropy.table import Table, Column
import os
import logging

from banzai.calibrations import CalibrationMaker, create_master_calibration_header
from banzai.utils import file_utils
from banzai_nres.images import ImageBase
import banzai_nres.settings as nres_settings
import banzai.settings as banzai_settings


logger = logging.getLogger(__name__)


class BlazeMaker(CalibrationMaker):
    def __init__(self, runtime_context):
        super(BlazeMaker, self).__init__(runtime_context)
        self.runtime_context = runtime_context
        self.trace_table_name = nres_settings.BLAZE_TABLE_NAME

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        blazes = []
        for image in images:
            master_header = create_master_calibration_header(old_header=image.header, images=[image])
            master_header['OBSTYPE'] = self.calibration_type
            master_filename = banzai_settings.CALIBRATION_FILENAME_FUNCTIONS[self.calibration_type](image)
            master_filepath = self._get_filepath(self.runtime_context, image, master_filename)
            logger.info('Making blaze file', image=image)
            lampflat_spectrum = image.data_tables[nres_settings.BOX_SPECTRUM_EXTNAME]
            blazes.append(blazes)
            logger.info('Created master blaze file', image=image, extra_tags={'calibration_type': self.calibration_type,
                                                                              'output_path': master_filepath,
                                                                              'calibration_obstype': master_header['OBSTYPE']})
        return blazes

    @staticmethod
    def _get_filepath(runtime_context, lampflat_image, master_filename):
        output_directory = file_utils.make_output_directory(runtime_context, lampflat_image)
        return os.path.join(output_directory, os.path.basename(master_filename))

    def do_stage(self, images):
        """
        :param images: list of images to run TraceMaker on.
        :return: [first_trace, second_trace,...] etc. This is a munge function which fixes the following problem.
        BlazeMaker.make_master_calibration_frame returns [first_blaze, second_blaze,...] and then BANZAI's do_stage
        wraps that list in another list, e.g. do_stage natively would return [[first_blaze, second_blaze,...]]. This
        rewrite of do_stage simply returns [[first_blaze, second_blaze,...]][0] = [first_blaze, second_blaze,...]
        """
        master_calibrations = super(BlazeMaker, self).do_stage(images)
        if len(master_calibrations) > 0:
            master_calibrations = master_calibrations[0]
        return master_calibrations


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
