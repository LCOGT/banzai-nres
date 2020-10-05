import numpy as np
from scipy.ndimage import binary_dilation, percentile_filter, binary_fill_holes
from scipy.stats import median_absolute_deviation
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import logging

logger = logging.getLogger('banzai')


def mark_features(flux, sigma=3, continuum_formal_error=None, detector_resolution=4, percentile=10, binary_dilations=3):
    """
    :param flux:
    :param line_widths: the standard deviation of the gaussian profile if profile='gaussian', else the width
                        of the lorentz_profile if profile='lorentzian'. See lorentz_profile below.
    :param min_snr_at_line_center:
    :param continuum_formal_error:
    :param detector_resolution:
    :return:
    """
    # start with a mask that marks every pixel as "ignore"
    mask = np.ones_like(flux, dtype=bool)
    if continuum_formal_error is None:
        # get the noisy scatter in the continuum-only portion of the spectrum by high-pass filtering the spectrum and
        # then taking the mad times 1.4826 (the robust standard deviation)
        continuum_formal_error = 1.4826 * median_absolute_deviation(flux - gaussian_filter1d(flux, sigma=detector_resolution/2))
    # take the brightest 10% (percentile=10% by default) of pixels to trace above the tops of the continuum
    continuum_estimate = percentile_filter(flux, percentile=-percentile, size=3*detector_resolution, mode='nearest')
    # keep the pixels that are close to the continuum estimate (mark as "keep")
    mask[np.isclose(continuum_estimate, flux, atol=sigma*continuum_formal_error)] = 0
    # binary dilate to cover the wings of lines which were clipped.
    mask = binary_dilation(mask, iterations=binary_dilations)
    if np.count_nonzero(mask) == len(mask):
        logger.warning('Masking unsuccessful. Entire spectral region was masked. Aborting and returning the spectrum'
                       'as it was input, i.e. without masking any elements.')
        # if all the elements are masked, raise a warning then return a all zeros mask (no masked elements)
        # so that later stages do not crash.
        mask = np.zeros_like(mask)
    return mask


def smooth_derivative(signal, smoothing_kernal):
    # smooth the data so that numerical differentiation does not amplify the noise
    smoothed_signal = gaussian_filter1d(signal, sigma=smoothing_kernal)
    first_derivative = np.hstack([[0], smoothed_signal[1:] - smoothed_signal[:-1]])
    return first_derivative


def lorentz_profile(x, center, width):
    return 1 / np.pi * 1 / 2 * width / ((x - center) ** 2 + 1 / 4 * width ** 2)
