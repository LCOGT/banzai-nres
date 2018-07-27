from banzai.utils.image_utils import get_bpm

from banzai.munge import munge

from banzai_nres.images import Image
from banzai.utils import fits_utils

from banzai import logs


logger = logs.get_logger(__name__)

def read_images(image_list, pipeline_context):
    """
    This is a copy of banzai.images.read_images
    which will properly handle images which already have a Bad Pixel Mask (BPM).
    Once this is fixed upstream in banzai, this should be deleted.
    The Image class will not work in upstream banzai for NRES frames (as is) because
    it fails to find the instrument in the DB.
    """

    images = []
    for filename in image_list:
        try:
            image = Image(pipeline_context, filename=filename)
            munge(image, pipeline_context)
            if image.bpm is None:
                logger.info('tele id and ccdsum: ' + str(image.telescope_id) + ' ,ccdsum:' + str(image.ccdsum))
                logger.info('instrument name and site:' + str(image.instrument) + str(image.site))
                bpm = get_bpm(image, pipeline_context)
                logger.info(str(bpm.shape))
                if bpm is None:
                    logger.error('No BPM file exists for this image.',
                                 extra={'tags': {'filename': image.filename}})
                else:
                    image.bpm = bpm
                    logger.info('images good')
                    images.append(image)
                    logger.info(str(image.bpm.shape))
            else:
                images.append(image)
                logger.info(str(image.bpm.shape))
        except Exception as e:
            logger.error('Error loading {0}'.format(filename))
            logger.error(e)
            continue
    return images
