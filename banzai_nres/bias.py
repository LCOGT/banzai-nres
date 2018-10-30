from banzai.bias import BiasMaker as BanzaiBiasMaker
import logging

logger = logging.getLogger(__name__)


class BiasMaker(BanzaiBiasMaker):

    def __init__(self, pipeline_context):
        super(BiasMaker, self).__init__(pipeline_context)

    @property
    def min_images(self):
        return 5
