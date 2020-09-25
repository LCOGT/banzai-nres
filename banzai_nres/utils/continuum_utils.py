import numpy as np
from scipy.ndimage import binary_dilation
from scipy.stats import median_absolute_deviation


def mark_absorption_or_emission_features(flux, line_width):
    """
    :param flux: array of flux in the spectrum
    :param line_width: the approximate widths of sharp absorption lines in the spectrum
    :return: mask: array, size of flux. 0 represents an element to keep, 1 represents an element to throw out (to mask)
    """
    mask = np.zeros_like(flux, dtype=int)
    first_derivative = np.hstack([[0], flux[1:] - flux[:-1]])
    mad = median_absolute_deviation(first_derivative)
    # 1.5 mad is roughly 1 sigma, which is very aggressive. However, empirically 1.5 mad works well on absorption line
    # spectra. See e.g. NRES absolute orders 75, 85, 95, and 110.
    features = np.abs(first_derivative) > 1.5 * mad
    # Mask features
    mask[features] = 1
    # Mask the immediate +- pixels around each absorption feature as well.
    # we binary dilate up to the width of the absorption lines so that we cover their wings as well.
    mask = binary_dilation(mask, iterations=line_width)
    return mask


def mark_pressure_broadened_features(flux, min_line_width):
    """
    :param mask: array, size of flux. 0 represents an element to keep, 1 represents an element to throw out (to mask)
    :param flux: array of flux in the spectrum
    :param min_line_width: the approximate minimum width of a pressure broadened line (in pixels)
    :return:
    """
    mask = np.zeros_like(flux, dtype=int)
    return mask
