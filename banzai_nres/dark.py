from banzai import logs
from banzai.dark import DarkMaker as BanzaiDarkMaker


logger = logs.get_logger(__name__)


class DarkMaker(BanzaiDarkMaker):

    def __init__(self, pipeline_context):
        super(DarkMaker, self).__init__(pipeline_context)

    @property
    def min_images(self):
        return 3
