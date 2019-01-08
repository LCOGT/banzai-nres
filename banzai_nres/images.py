from banzai.images import Image as BanzaiImage


class NRESImage(BanzaiImage):
    def __init__(self, pipeline_context, filename=None, data=None, data_tables=None,
                 header=None, extension_headers=None, bpm=None):
        super(NRESImage, self).__init__(pipeline_context, filename=filename, data=data, data_tables=data_tables,
                                        header=header, extension_headers=extension_headers, bpm=bpm)
        self.trace = None
