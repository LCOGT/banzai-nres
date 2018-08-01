from banzai.stages import Stage


class TraceFitOrderbyOrder(Stage):
    def __init__(self, pipeline_context):
        super(TraceFitOrderbyOrder, self).__init__(pipeline_context)
