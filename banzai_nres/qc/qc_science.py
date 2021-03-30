import numpy as np
from banzai.stages import Stage
import logging

logger = logging.getLogger('banzai')

# This defines the order for which we calculate the SNR value for the header.
# By default this is the order containing the Mg b lines.
SNR_ORDER = 90


class CalculateScienceFrameMetrics(Stage):

    def __init__(self, runtime_context):
        super(CalculateScienceFrameMetrics, self).__init__(runtime_context)

    def do_stage(self, image):
        snr = get_snr(image, SNR_ORDER)
        image.meta['SNR'] = snr, 'Signal-to-noise ratio at 5180 Angstroms'

        return image


def get_snr(image, order):
    snr_all = image.spectrum[image.science_fiber, order]['flux'] / image.spectrum[image.science_fiber, order]['uncertainty']
    snr = np.percentile(snr_all, 90)
    return snr
