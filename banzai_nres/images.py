from banzai.images import Image as BanzaiImage
from banzai_nres import traces


class Image(BanzaiImage):

    def __init__(self, pipeline_context, filename=None, data=None, header=None,
                 extension_headers=None, bpm=None):
        super(Image, self).__init__(pipeline_context, filename=filename, data=data, header=header,
                                    extension_headers=extension_headers, bpm=bpm)

        self.trace_fit_coefficients = None
        self.fiber_order = None
