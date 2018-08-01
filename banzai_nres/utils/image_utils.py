from banzai.munge import munge as banzai_munge
from banzai import logs, dbs
from banzai.utils import image_utils

from banzai.images import Image


logger = logs.get_logger(__name__)


def munge(image, pipeline_context):
    image.header['CRPIX1'] = 0
    image.header['CRPIX2'] = 0
    # adding headers which get called inside of _trim_image of trim.py in Banzai.
    banzai_munge(image, pipeline_context)


def read_images(image_list, pipeline_context):
    """
    This is a copy of banzai.images.read_images
    which includes the above munge patch.
    """
    images = []
    for filename in image_list:
        try:
            logger.info('above Image call')
            image = Image(pipeline_context, filename=filename)
            logger.info('created the Image object')
            if image.telescope is None:
                logger.info('None telescope')
                error_message = 'Telescope is not in the database: {site}/{instrument}'
                error_message = error_message.format(site=image.site, instrument=image.instrument)
                raise dbs.TelescopeMissingException(error_message)
            munge(image, pipeline_context)
            if image.bpm is None:
                image_utils.load_bpm(image, pipeline_context)
            images.append(image)
        except Exception as e:
            logger.error('Error loading image: {error}'.format(error=e), extra={'tags': {'filename': filename}})
            continue

    return images
