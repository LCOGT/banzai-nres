import numpy as np
from scipy.ndimage import binary_dilation, percentile_filter
from scipy.stats import median_absolute_deviation
from scipy.ndimage import gaussian_filter1d
import logging

logger = logging.getLogger('banzai')


def mark_features(flux, sigma=3, continuum_formal_error=None, detector_resolution=4):
    """
    :param flux:
    :param sigma:
        :param continuum_formal_error:
    :param min_feature_width:
    :param detector_resolution:
    :return:
    """
    # start with a mask that marks every pixel as "ignore"
    mask = np.ones_like(flux, dtype=bool)
    if continuum_formal_error is None:
        # get the noisy scatter in the continuum-only portion of the spectrum by high-pass filtering the spectrum and
        # then taking the mad times 1.4826 (the robust standard deviation)
        continuum_formal_error = 1.4826 * median_absolute_deviation(flux -
                                                                    gaussian_filter1d(flux,
                                                                                      sigma=detector_resolution/2))
    # take the brightest 10% of pixels to trace above the tops of the continuum
    continuum_estimate = percentile_filter(flux, percentile=-10, size=3 * detector_resolution, mode='nearest')
    # keep the pixels that are close to the continuum estimate (mark as "keep")
    mask[np.isclose(continuum_estimate, flux, atol=sigma*continuum_formal_error)] = 0
    # binary dilate to cover the wings of lines which were clipped.
    mask = binary_dilation(mask, iterations=detector_resolution)
    if np.count_nonzero(mask) == len(mask):
        logger.warning('Masking unsuccessful. Entire spectral region was masked. Aborting and returning the spectrum'
                       'as it was input, i.e. without masking any elements.')
        # if all the elements are masked, raise a warning then return a all zeros mask (no masked elements)
        # so that later stages do not crash.
        mask = np.zeros_like(mask)
    return mask
