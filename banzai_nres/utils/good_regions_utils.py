import numpy as np
from .trace_utils import get_trace_centroids_from_coefficients
from scipy import ndimage, signal
import matplotlib.pyplot as plt


def identify_high_SN_region_bounds_along_traces(image, sig_to_noise_threshold=20):
    """
    :param image: banzai_nres image object with Trace() object attached. Image.data should be a lampflat.
    :param sig_to_noise_threshold: threshold value of trace_centroid_flux/np.sqrt(variance_of_pixel_nearest_centroid)
                                    whereby values beneath that constrain the 'good' region of the detector.
    :return: None
            sets image.trace.bright_x_region_bounds , ndarray num_traces x 2 : [[low_x_bound_for_trace1, high_x_bound],[],....]
            e.g. low_x_bound_for_trace1, high_x_bound = 0, 4095 for a high S/N trace.

    """
    trace_y_values_vs_x, num_traces, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
    rounded_trace_y_values_vs_x = trace_y_values_vs_x.astype(np.int)
    trace_x_coords = np.ones_like(rounded_trace_y_values_vs_x, dtype=np.int) * x
    indices_of_trace_centers = [rounded_trace_y_values_vs_x.flatten(), trace_x_coords.flatten()]

    flux_values = ndimage.map_coordinates(image.data.astype(np.float64), indices_of_trace_centers, order=0, prefilter=False)
    ivars = ndimage.map_coordinates(image.ivar, indices_of_trace_centers, order=0, prefilter=False)

    signal_to_noise = (flux_values * np.sqrt(ivars)).reshape(trace_x_coords.shape)
    window = np.ones((1, 5))/5
    # light filtering along the rows (not columns!) to remove extraneous high S/N points in an otherwise low S/N trace.
    signal_to_noise_filt = signal.convolve2d(signal_to_noise, window, mode='same')
    good_values_mask = signal_to_noise_filt > sig_to_noise_threshold
    good_trace_x_coords = np.zeros_like(trace_x_coords, dtype=np.float64)
    good_trace_x_coords[good_values_mask] = trace_x_coords[good_values_mask]
    good_trace_x_coords[~good_values_mask] = np.nan

    bright_lower_x_bounds = np.nanmin(good_trace_x_coords, axis=1)
    bright_upper_x_bounds = np.nanmax(good_trace_x_coords, axis=1)
    # any nan values mean there are no good values along that trace. Set them to zero so that later stages
    # can identify that those traces are too dim, and act accordingly.
    bright_lower_x_bounds[np.isnan(bright_lower_x_bounds)] = 0
    bright_upper_x_bounds[np.isnan(bright_upper_x_bounds)] = 0

    good_x_bounds = (np.stack((bright_lower_x_bounds, bright_upper_x_bounds)).T).astype(np.int)

    return good_x_bounds


def flag_traces_with_insufficient_high_SN_region(image, min_region_size=500):
    """
    :param image:
    :param min_region_size:

            sets image.trace.has_sufficient_signal_to_noise, ndarray type bool, length of the number of traces.
            image.trace.has_sufficient_signal_to_noise[i] = True if the span of the high SN (bright region) in x
            is larger than min_region_size in pixels
    """
    bright_upper_bounds = image.trace.high_signal_to_noise_region_bounds[:, 1]
    bright_lower_bounds = image.trace.high_signal_to_noise_region_bounds[:, 0]
    return (bright_upper_bounds - bright_lower_bounds > min_region_size)