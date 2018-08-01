import os.path

import numpy as np

from banzai.images import Image
from banzai import logs
from banzai.bias import BiasMaker as BanzaiBiasMaker
from banzai.utils import stats, fits_utils


logger = logs.get_logger(__name__)


class BiasMaker(BanzaiBiasMaker):

    def __init__(self, pipeline_context):
        super(BiasMaker, self).__init__(pipeline_context)

    @property
    def min_images(self):
        return 5
