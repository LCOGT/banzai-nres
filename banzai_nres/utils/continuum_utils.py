import numpy as np
from scipy.ndimage import binary_dilation
from scipy.stats import median_absolute_deviation
from scipy.ndimage import gaussian_filter1d


def mark_absorption_or_emission_features(flux, line_width, sigma=1):
    """
    :param flux: array of flux in the spectrum
    :param line_width: the approximate widths of absorption or emission lines in the spectrum that you want to reject.
    :param sigma: approximate signal to noise a line has to have to be rejected. Keep in mind for a feature to be rejected
    it has to be both have a signal-to-noise (s/n) > sigma AND have a width comparable to line_width. Therefore single pixel fluctuations
    with a s/n of order 1 sigma will not be rejected because their widths are too small (unless you set line_width to 1 pixel)
    :return: mask: array, size of flux. 0 represents an element to keep, 1 represents an element to throw out (to mask)
    """
    mask = np.zeros_like(flux, dtype=bool)
    first_derivative = smooth_derivative(flux, max(line_width//2, 1))
    mad = median_absolute_deviation(first_derivative)
    # 1.4826 mad is 1 sigma for a Gaussian distribution. 1 sigma is very aggressive.
    # However, empirically 1 sigma with a line width = to the detector resolution (e.g. 4 pixels)
    # works well on absorption line spectra on NRES. See e.g. NRES absolute orders 75, 85, 95, and 110.
    features = np.abs(first_derivative) > 1.4826 * sigma * mad
    # Mask features
    mask[features] = 1
    # Mask the immediate +- pixels around each absorption feature as well.
    # we binary dilate up to the width of the absorption lines so that we cover their wings as well.
    mask = binary_dilation(mask, iterations=line_width)
    return mask


def smooth_derivative(signal, smoothing_kernal):
    # smooth the data so that numerical differentiation does not amplify the noise
    smoothed_signal = gaussian_filter1d(signal, sigma=smoothing_kernal)
    first_derivative = np.hstack([[0], smoothed_signal[1:] - smoothed_signal[:-1]])
    return first_derivative
