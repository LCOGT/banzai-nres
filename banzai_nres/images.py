from banzai.images import Image as BanzaiImage
from banzai_nres.dbs import get_trace_coefficients


class Image(BanzaiImage):

    def __init__(self, pipeline_context, filename=None, data=None, header=None,
                 extension_headers=None, bpm=None):
        super(Image, self).__init__(pipeline_context, filename=filename, data=data, header=header,
                                    extension_headers=extension_headers, bpm=bpm)

        if self.header is not None:
            self.trace_fit_coefficients, self.fiber_order = get_trace_coefficients(self)
