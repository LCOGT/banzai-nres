from banzai.images import Image
from banzai_nres.fibers import FibersState


class NRESImage(Image):
    def __init__(self, pipeline_context, filename=None, data=None, data_tables=None,
                 header=None, extension_headers=None, bpm=None):
        super(NRESImage, self).__init__(pipeline_context, filename=filename, data=data, data_tables=data_tables,
                                        header=header, extension_headers=extension_headers, bpm=bpm)
        self.fibers_state = FibersState.from_header(self.header)
