from __future__ import absolute_import, division, print_function, unicode_literals

import os.path

import numpy as np

from banzai.images import Image
from banzai import logs
from banzai.stages import CalibrationMaker
from banzai.utils import stats, fits_utils

from scipy.ndimage import filters



class BiasMaker(CalibrationMaker):

    def __init__(self, pipeline_context):
        super(BiasMaker, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'BIAS'

    @property
    def min_images(self):
        return 5

    def make_master_calibration_frame(self, images, image_config, logging_tags):

        bias_data = np.zeros((image_config.ny, image_config.nx, len(images)), dtype=np.float32)
        bias_mask = np.zeros((image_config.ny, image_config.nx, len(images)), dtype=np.uint8)
        bias_level_array = np.zeros(len(images), dtype=np.float32)

        master_bias_filename = self.get_calibration_filename(image_config)
        logs.add_tag(logging_tags, 'master_bias', os.path.basename(master_bias_filename))
        for i, image in enumerate(images):
            bias_level_array[i] = stats.sigma_clipped_mean(image.data, 3.5, mask=image.bpm)

            logs.add_tag(logging_tags, 'filename', os.path.basename(image.filename))
            logs.add_tag(logging_tags, 'BIASLVL', float(bias_level_array[i]))
            self.logger.debug('Calculating bias level', extra=logging_tags)
            # Subtract the bias level for each image
            bias_data[:, :, i] = image.data[:, :] - bias_level_array[i]
            bias_mask[:, :, i] = image.bpm[:, :]

        mean_bias_level = stats.sigma_clipped_mean(bias_level_array, 3.0)

        master_bias = stats.sigma_clipped_mean(bias_data, 3.0, axis=2, mask=bias_mask, inplace=True)

        del bias_data
        del bias_mask

        master_bpm = np.array(master_bias == 0.0, dtype=np.uint8)

        header = fits_utils.create_master_calibration_header(images)

        header['BIASLVL'] = (mean_bias_level, 'Mean bias level of master bias')
        master_bias_image = Image(self.pipeline_context, data=master_bias, header=header)
        master_bias_image.filename = master_bias_filename
        master_bias_image.bpm = master_bpm

        logs.pop_tag(logging_tags, 'master_bias')
        logs.add_tag(logging_tags, 'filename', os.path.basename(master_bias_image.filename))
        logs.add_tag(logging_tags, 'BIASLVL', mean_bias_level)
        self.logger.debug('Average bias level in ADU', extra=logging_tags)

        return [master_bias_image]