import numpy as np
import abc
import logging

from banzai_nres.utils.extract_utils import Extract
from banzai_nres.utils import extract_utils

from banzai.stages import Stage


logger = logging.getLogger(__name__)


class BoxExtract(Extract):
    def __init__(self, rectified_spectrum_per_order=None):
        self.rectified_spectrum_per_order = rectified_spectrum_per_order

    @abc.abstractmethod
    def extract(self):
        extracted_spectrum_per_order = {}
        # this relies on rectified spectrum being a dict, which will not be the case.
        # this class still seems clunky. Do we need it? Can we wrap this action into the BoxExtractor class?
        for order_id in list(self.rectified_spectrum_per_order.keys()):
            extracted_spectrum_per_order[order_id] = self.extract_order(self.rectified_spectrum_per_order[order_id])
        return extracted_spectrum_per_order


class BoxExtractor(Stage):
    def __init__(self, pipeline_context):
        super(BoxExtractor, self).__init__(pipeline_context)

    def do_stage(self, image):
        logger.info('Box extracting spectrum', image=image)
        if image.extracted_spectrum is None:
            image.extracted_spectrum = {}
        image.extracted_spectrum['box'] = BoxExtract(image.rectified_data).extract()
        # we need to append the extracted spectrum onto the image in a more intelligent way.
        # maybe we should append this as a table, then it can be saved through banzai in the same way as other tables,
        # and then we can clearly label wavelength etc.
        return image


class RectifyTwodSpectrum(Stage):
    def __init__(self, pipeline_context):
        super(RectifyTwodSpectrum, self).__init__(pipeline_context)
        self.max_extraction_half_window = 10

    def do_stage(self, image):
        logger.info('Rectifying the 2d spectrum', image=image)
        if image.trace is None:
            logger.error('Image has empty trace attribute', image=image)
            raise ValueError('image.trace is None')
        rectified_data = extract_utils.rectify_orders(image.data, image.trace,
                                                      half_window=self.max_extraction_half_window)
        image.rectified_data = rectified_data
        # rectified_data should probably be a class like Trace,
        # which has a get_spectrum for an order id as an attribute etc.
        # also we need to append rectified_data onto the image in a more intelligent way.
        return image
