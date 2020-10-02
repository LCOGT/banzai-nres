import numpy as np
from scipy.ndimage import binary_dilation
from scipy.stats import median_absolute_deviation
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.signal.windows import gaussian


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
    # estimate the formal error in the continuum component
    continuum_formal_error = 1.4826 * median_absolute_deviation(first_derivative)
    # 1.4826 mad is 1 sigma for a Gaussian distribution. 1 sigma is very aggressive.
    # However, empirically 1 sigma with a line width = to the detector resolution (e.g. 4 pixels)
    # works well on absorption line spectra on NRES. See e.g. NRES absolute orders 75, 85, 95, and 110.
    # Account for any offset term in the first_derivative (i.e. a constant slope in the continuum)
    # by subtracting off the median of the first_derivative
    features = np.abs(first_derivative - np.median(first_derivative)) > sigma * continuum_formal_error
    # Mask features
    mask[features] = 1
    # Mask the immediate +- pixels around each absorption feature as well.
    # we binary dilate up to the width of the absorption lines so that we cover their wings as well.
    mask = binary_dilation(mask, iterations=line_width)
    return mask


def mark_features_test(_flux, error, line_widths, kernal_choice='gaussian', min_snr_at_line_center=3, type='absorption'):
    mask = np.zeros_like(_flux, dtype=bool)
    continuum_formal_error = 1.4826 * median_absolute_deviation(_flux)
    flux = _flux - np.linspace(100, 200, len(_flux))
    for width in np.array([line_widths]).flatten():
        this_mask = np.zeros_like(flux, dtype=bool)
        if kernal_choice.lower() == 'gaussian' or kernal_choice.lower() == 'guassian':
            kernal = gaussian(2 * width * 3, width)  # go out to 3 sigma to capture 99.73% of the guassian kernal
        elif kernal_choice.lower() == 'lorentzian':
            bound = min(width*3, len(flux))
            kernal = lorentz_profile(np.arange(-bound, bound+1), center=0, width=width)
        else:
            raise ValueError(f'kernal_choice {kernal_choice} is not implemented.')
        kernal /= np.max(kernal)  # make sure the kernal is normalized at its maximum
        if type == 'absorption':
            kernal *= -1
        # gives the ccf in terms of the peak S/N at line center.
        # E.g. if flux/error is a guassian absorption/emission line with a S/N at center of 100, then ccf
        # will have a maximum at that line's position with a value of 100.
        ccf = np.correlate(flux/error, kernal, mode='same')/(kernal**2).sum()
        # remove edge effects as a consequence of the 'same' mode
        ccf[-len(kernal)//2:] = np.nan
        ccf[:len(kernal)//2] = np.nan
        import matplotlib.pyplot as plt
        plt.plot(ccf)
        plt.plot(flux)
        plt.show()
        this_mask[np.logical_and(np.abs(ccf) > min_snr_at_line_center, np.isfinite(ccf))] = 1
        mask = np.logical_or(mask, binary_dilation(this_mask, iterations=int(width/2)))
    return mask


def mark_features(input_flux, line_widths, min_snr_at_line_center=3, type='absorption', profile='gaussian'):
    mask = np.zeros_like(input_flux, dtype=bool)
    # smooth away noise fluctuations
    min_feature_width = np.min(line_widths)//2
    flux = gaussian_filter1d(input_flux, sigma=min_feature_width)
    continuum_formal_error = 1.4826 * median_absolute_deviation(input_flux - flux)
    if type == 'absorption':
        flux *= -1
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
