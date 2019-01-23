from banzai.images import Image
from banzai_nres.fibers import fiber_states_from_header


class NRESImage(Image):
    def __init__(self, pipeline_context, filename=None, data=None, data_tables=None,
                 header=None, extension_headers=None, bpm=None):
        super(NRESImage, self).__init__(pipeline_context, filename=filename, data=data, data_tables=data_tables,
                                        header=header, extension_headers=extension_headers, bpm=bpm)
        self.trace = None
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = fiber_states_from_header(self.header)
