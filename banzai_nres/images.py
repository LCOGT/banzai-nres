from banzai.images import Image as banzaiImage
from banzai import dbs

from banzai import logs


logger = logs.get_logger(__name__)


class Image(banzaiImage):

    def __init__(self, pipeline_context, per_pixel_variance=None, filename=None, data=None, header={},
                 extension_headers=[], bpm=None):
        super(Image, self).__init__(pipeline_context, filename=filename, data=data, header=header,
                                    extension_headers=extension_headers, bpm=bpm)
        self.per_pixel_variance = per_pixel_variance
        self.instrument = header.get('TELESCOP')

        if self.site is not None and self.instrument is not None:
            self.telescope_id = dbs.get_telescope_id(self.site, self.instrument,
                                                     db_address=pipeline_context.db_address)
        else:
            self.telescope_id = None

        logger.info('should be nres01:' + header.get('TELESCOP'), extra={'tags': {'raw_path': pipeline_context.raw_path}})
        logger.info('telescope id:' + self.telescope_id, extra={'tags': {'raw_path': pipeline_context.raw_path}})