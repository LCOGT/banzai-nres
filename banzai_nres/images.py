from banzai.images import Image as banzaiImage

from banzai import logs


logger = logs.get_logger(__name__)


class Image(banzaiImage):

    def __init__(self, pipeline_context, per_pix_vars=None, filename=None, data=None, header={},
                 extension_headers=[], bpm=None):
        super(Image, self).__init__(pipeline_context, filename=filename, data=data, header=header,
                                    extension_headers=extension_headers, bpm=bpm)
        self.per_pix_vars = per_pix_vars
        self.instrument = header.get('TELESCOP')