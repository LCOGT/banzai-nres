import numpy as np
from banzai.stages import Stage
from banzai.logs import get_logger

logger = get_logger()

# This defines the order for which we calculate the SNR value for the header.
# By default this is the order containing the Mg b lines.
SNR_ORDER = 90


class CalculateScienceFrameMetrics(Stage):

    def __init__(self, runtime_context):
        super(CalculateScienceFrameMetrics, self).__init__(runtime_context)

    def do_stage(self, image):
        snr, snr_wave = get_snr(image, SNR_ORDER, self.runtime_context.PIXELS_PER_RESOLUTION_ELEMENT)

        image.meta['SNR'] = snr, 'Signal-to-noise ratio/res element @ 5180 Angstrom'

        image.meta['SCIFIBER'] = image.science_fiber, 'Fiber ID of the science spectrum'

        return image


def get_snr(image, order, pixels_per_resolution_element):
    snr_all = image.spectrum[image.science_fiber, order]['flux'] / \
              image.spectrum[image.science_fiber, order]['uncertainty']
    snr = np.percentile(snr_all, 90)
    snr_wave = np.mean(image.spectrum[image.science_fiber, order]['wavelength'][snr_all > snr])
    snr *= np.sqrt(pixels_per_resolution_element)
    return snr, snr_wave
