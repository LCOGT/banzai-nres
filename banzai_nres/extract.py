import numpy as np
import abc
import logging
from astropy.table import Table

from banzai_nres.utils.extract_utils import Extract
from banzai_nres.utils import extract_utils
import banzai_nres.settings as nres_settings

from banzai.stages import Stage
from banzai.images import DataTable


logger = logging.getLogger(__name__)


class BoxExtract(Extract):
    def __init__(self, rectified_spectrum_per_order=None):
        self.rectified_spectrum_per_order = rectified_spectrum_per_order

    @abc.abstractmethod
    def extract(self):
        extracted_spectrum_per_order = {'id': [], 'flux': [], 'pixel': []}
        for order_id in list(self.rectified_spectrum_per_order.keys()):
            flux = self.extract_order(self.rectified_spectrum_per_order[order_id])
            extracted_spectrum_per_order['flux'].append(flux)
            extracted_spectrum_per_order['pixel'].append(np.arange(len(flux)))
            extracted_spectrum_per_order['id'].append(order_id)
        return Table(extracted_spectrum_per_order)


class BoxExtractor(Stage):
    def __init__(self, pipeline_context):
        super(BoxExtractor, self).__init__(pipeline_context)
        self.extraction_half_window = nres_settings.BOX_EXTRACTION_HALF_WINDOW
        self.max_extraction_half_window = nres_settings.MAX_EXTRACTION_HALF_WINDOW

    def do_stage(self, image):
        logger.info('Box extracting spectrum', image=image)
        rectified_twod_spectrum = self._trim_rectified_data(image.rectified_data)
        spectrum = BoxExtract(rectified_twod_spectrum).extract()
        image.data_tables['box_extracted_spectrum'] = DataTable(data_table=spectrum, name='SPECBOX')
        return image

    def _trim_rectified_data(self, rectified_twod_spectrum):
        """
        :param rectified_twod_spectrum: dictionary of 2d spectra, where each index order_id gives a spectrum in a window
        of size 2 * MAX_EXTRACTION_HALF_WINDOW + 1 around the trace labelled by order_id.
        :return: dictionary of 2d spectra, where each index order_id gives a spectrum in a window of size
        2 * BOX_EXTRACTION_HALF_WINDOW + 1 around the trace labelled by order_id.

        NOTE: this assumes implicitly that the center of the rectified traces lies at the array index
        MAX_EXTRACTION_HALF_WINDOW in each of the elements of rectified_twod_spectrum. I.e. the center of the trace
        is always at the center of the rectified spectrum.
        """
        trimmed_rectified_spectrum = {}
        if self.extraction_half_window >= self.max_extraction_half_window:
            # short circuit
            logger.warning('Box extraction window was chosen to be >= the max extraction window '
                           '(which was used in rectification). Defaulting to the max extraction window.')
            return rectified_twod_spectrum
        for order_id in list(rectified_twod_spectrum.keys()):
            trim = self.max_extraction_half_window - self.extraction_half_window
            trimmed_rectified_spectrum[order_id] = rectified_twod_spectrum[order_id][trim:-trim]
        return trimmed_rectified_spectrum


class RectifyTwodSpectrum(Stage):
    def __init__(self, pipeline_context):
        super(RectifyTwodSpectrum, self).__init__(pipeline_context)
        self.max_extraction_half_window = nres_settings.MAX_EXTRACTION_HALF_WINDOW

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
