import numpy as np
from scipy.ndimage import binary_dilation
from scipy.stats import median_absolute_deviation
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.signal.windows import gaussian


def mark_features(input_flux, line_widths, min_snr_at_line_center=3, type='absorption', profile='gaussian'):
    mask = np.zeros_like(input_flux, dtype=bool)
    # smooth away noise fluctuations
    min_feature_width = np.min(line_widths)//2
    flux = gaussian_filter1d(input_flux, sigma=min_feature_width)
    # get the noisy scatter in the continuum-only portion of the spectrum by high-pass filtering the spectrum and
    # then taking the mad times 1.4826 (the robust standard deviation)
    continuum_formal_error = 1.4826 * median_absolute_deviation(input_flux - flux)
    if type == 'absorption':
        flux *= -1
    # pressure broadened lines have less steep wings and so require more binary dilations to capture those wings
    binary_dilations = {'gaussian': 1.5, 'lorentzian': 3}[profile]
    flux -= np.median(flux)
    for width in np.array([line_widths]).flatten():
        this_mask = np.zeros_like(flux, dtype=bool)
        peaks = find_peaks(flux, distance=min_feature_width, prominence=continuum_formal_error * min_snr_at_line_center)[0]
        this_mask[peaks] = 1
        mask = np.logical_or(mask, binary_dilation(this_mask, iterations=int(binary_dilations * width)))
    return mask


def smooth_derivative(signal, smoothing_kernal):
    # smooth the data so that numerical differentiation does not amplify the noise
    smoothed_signal = gaussian_filter1d(signal, sigma=smoothing_kernal)
    first_derivative = np.hstack([[0], smoothed_signal[1:] - smoothed_signal[:-1]])
    return first_derivative


def lorentz_profile(x, center, width):
    return 1 / np.pi * 1 / 2 * width / ((x - center) ** 2 + 1 / 4 * width ** 2)
