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
    """

    images = []
    for filename in image_list:
        try:
            logger.info('in banzai_nres read', extra={'tags': {'raw_path': pipeline_context.raw_path}})
            fits_utils.open_image(filename)  #
            logger.info('file opened', extra={'tags': {'raw_path': pipeline_context.raw_path}})
            image = Image(pipeline_context, filename=filename)
            logger.info('built image', extra={'tags': {'inst,site,teleid': image.instrument + image.site + str(image.telescope_id)}})
            munge(image, pipeline_context)
            if image.bpm is None:
                bpm = get_bpm(image, pipeline_context)
                if bpm is None:
                    logger.error('No BPM file exists for this image.',
                                 extra={'tags': {'filename': image.filename}})
                else:
                    image.bpm = bpm
                    images.append(image)
            else:
                images.append(image)
        except Exception as e:
            logger.error('Error loading {0}'.format(filename))
            logger.error(e)
            continue
    return images
