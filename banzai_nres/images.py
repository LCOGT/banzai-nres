from banzai.images import Image as BanzaiImage


class Image(BanzaiImage):

    def __init__(self, pipeline_context, filename=None, data=None, header=None,
                 extension_headers=None, bpm=None):
        super(Image, self).__init__(pipeline_context, filename=filename, data=data, header=header,
                                    extension_headers=extension_headers, bpm=bpm)

        self.trace = None
