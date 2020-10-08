import numpy as np
from scipy.ndimage import binary_dilation, percentile_filter, binary_fill_holes
from scipy.stats import median_absolute_deviation
from scipy.ndimage import gaussian_filter1d, label
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
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
        continuum_formal_error = 1.4826 * median_absolute_deviation(flux - gaussian_filter1d(flux, sigma=detector_resolution/2))
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


def mark_broad_features(flux, flux_error, mask, broad_line_min_width, level_to_mask=1E-2):
    broad_line_mask = np.zeros_like(mask)
    # find how many widths corresponds to the percent level of the profile to mask away
    profile = lorentz_profile(np.arange(1000), 500, 10)
    num_widths = np.abs(((np.arange(1000) - 500) / 10)[np.argmin(np.abs(profile/np.max(profile) - level_to_mask))])
    # find features that have widths of at least broad_line_min_width
    broad_regions = contiguous_regions(mask, broad_line_min_width)
    # fit lorentz profiles to the pressure broadened lines
    x = np.arange(len(flux))
    for region in broad_regions:
        # remove the offset in the flux so that we do not have to fit a fourth parameter (an offset)
        flux_to_fit = flux[region] - np.percentile(flux[region], 95)
        # fit a lorentzian profile to the offset flux
        center, width, amplitude = np.mean(x[region]), broad_line_min_width, np.percentile(flux_to_fit, 10)
        popt, pcov = curve_fit(lorentz_profile, x[region], flux_to_fit, p0=[center, width, amplitude],
                               sigma=flux_error[region], absolute_sigma=True)
        center, width, amplitude = popt[0], popt[1], popt[2]
        #dof = len(flux_to_fit) - 3  # degrees of freedom
        #norm_chi2 = np.sum(((lorentz_profile(x[region], center, width, amplitude) - flux_to_fit)/flux_error[region])**2)/dof
        #import matplotlib.pyplot as plt
        #plt.plot(x[region], flux_to_fit)
        #plt.plot(x[region], lorentz_profile(x[region], center, width, amplitude))
        #plt.title(f'{width} , {norm_chi2}')
        #plt.show()
        # vet the fit:
        if center > np.max(x[region]) or center < np.min(x[region]):
            # if the fitted center is outside of the original masked region, something went wrong.
            continue
        if width < broad_line_min_width/10:
            # if the fitted width is substantially smaller than the expected width, then something went wrong.
            continue
        # mask away +-num_widths to remove the profiles up to level_to_mask
        lower = max(0, int(center - num_widths*width))
        upper = min(int(center + num_widths*width), len(mask))
        broad_line_mask[lower:upper] = 1
    return broad_line_mask


def contiguous_regions(a, min_structure_size):
    regions, num_regions = label(a)
    broad_regions = []
    for region_label in range(1, num_regions + 1):
        region_mask = regions == region_label
        if np.count_nonzero(region_mask) >= min_structure_size:
            broad_regions.append(region_mask)
    return broad_regions


def smooth_derivative(signal, smoothing_kernal):
    # smooth the data so that numerical differentiation does not amplify the noise
    smoothed_signal = gaussian_filter1d(signal, sigma=smoothing_kernal)
    first_derivative = np.hstack([[0], smoothed_signal[1:] - smoothed_signal[:-1]])
    return first_derivative


def lorentz_profile(x, center, width, amplitude=1):
    return amplitude / np.pi * 1 / 2 * width / ((x - center) ** 2 + 1 / 4 * width ** 2)
