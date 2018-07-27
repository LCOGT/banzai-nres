import os.path

import numpy as np

from banzai_nres.images import Image
from banzai import logs
from banzai.dark import DarkMaker as BanzaiDarkMaker
from banzai.utils import stats, fits_utils


logger = logs.get_logger(__name__)


class DarkMaker(BanzaiDarkMaker):

    def __init__(self, pipeline_context):
        super(DarkMaker, self).__init__(pipeline_context)

    @property
    def min_images(self):
        return 3

    def make_master_calibration_frame(self, images, image_config, logging_tags):
        dark_data = np.zeros((images[0].ny, images[0].nx, len(images)), dtype=np.float32)
        dark_mask = np.zeros((images[0].ny, images[0].nx, len(images)), dtype=np.uint8)

        master_dark_filename = self.get_calibration_filename(images[0])

        logs.add_tag(logging_tags, 'master_dark', os.path.basename(master_dark_filename))
        for i, image in enumerate(images):
            logs.add_tag(logging_tags, 'filename', os.path.basename(image.filename))
            self.logger.debug('Combining dark', extra=logging_tags)

            dark_data[:, :, i] = image.data[:, :]
            dark_mask[:, :, i] = image.bpm[:, :]

        master_dark = stats.sigma_clipped_mean(dark_data, 3.0, axis=2, mask=dark_mask, inplace=True)

        # Memory cleanup
        del dark_data
        del dark_mask

        master_bpm = np.array(master_dark == 0.0, dtype=np.uint8)

        # Save the master dark image with all of the combined images in the header
        master_dark_header = fits_utils.create_master_calibration_header(images)
        master_dark_image = Image(self.pipeline_context, data=master_dark,
                                  header=master_dark_header)
        master_dark_image.filename = master_dark_filename
        master_dark_image.bpm = master_bpm

        logs.pop_tag(logging_tags, 'master_dark')
        logs.add_tag(logging_tags, 'filename', os.path.basename(master_dark_image.filename))
        self.logger.info('Created master dark', extra=logging_tags)

        return [master_dark_image]
