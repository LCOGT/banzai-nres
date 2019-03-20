import numpy as np
import abc

from banzai_nres.utils.extract_utils import Extract

import logging

logger = logging.getLogger(__name__)


class BoxExtract(Extract):
    def __init__(self, rectified_spectrum_per_order=None):
        self.rectified_spectrum_per_order = rectified_spectrum_per_order

    @abc.abstractmethod
    def extract(self):
        extracted_spectrum_per_order = {}
        for order_id in list(self.rectified_spectrum_per_order.keys()):
            extracted_spectrum_per_order[order_id] = self.extract_order(self.rectified_spectrum_per_order[order_id])
