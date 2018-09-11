from banzai_nres.utils.trace_utils import get_trace_centroids_from_coefficients
from banzai_nres.tests.utils import gaussian, noisify_image, trim_image, FakeImage
import numpy as np


def fill_image_with_traces(image, trimmed_shape, order_width=1.25, odd_fiber_intensity=1E4, even_fiber_intensity=5E3):
    """
    :param image: Banzai image object.
    :param trimmed_shape: The trimmed shape of the image, tuple (rows, columns)
    :param order_width: the standard deviation of an unnormalized guassian e.g. sigma
           from np.exp(-(x - b) ** 2 / (2 * sigma ** 2))
    :param odd_fiber_intensity: peak intensity of the 'odd' fiber (e.g. all fibers starting at index 1, stepping by 2)
    :param even_fiber_intensity: peak intensity of the 'even' fiber (e.g. all fibers starting at index 0, stepping by 2)
    Appends the bright streaks across the detector onto image.data according to image.trace.coefficients.
    """
    image.data = np.zeros_like(image.data)
    even_fiber = np.zeros(trimmed_shape)
    odd_fiber = np.zeros(trimmed_shape)

    trace_values_versus_xpixel, num_traces, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
    vgauss = np.vectorize(gaussian)  # prepare guassian for evaluation along a slice centered at each trace point.
    # these are realistic intensity values according to a 120 sec LSC exposure.
    for x_pixel in range(even_fiber.shape[1]):
        for i in range(num_traces):
            centroid = trace_values_versus_xpixel[i, x_pixel]
            low, high = max(0, int(centroid - 5 * order_width)), min(even_fiber.shape[1]-1, int(centroid + 5 * order_width)) + 1
            evalwindow = np.arange(low, high, 1)
            if len(evalwindow) > 0:
                if i % 2 == 0:
                    even_fiber[low:high, x_pixel] += vgauss(evalwindow, 1, centroid, order_width)
                else:
                    odd_fiber[low:high, x_pixel] += vgauss(evalwindow, 1, centroid, order_width)

    # add traces
    image.data[:trimmed_shape[0], :trimmed_shape[1]] += (odd_fiber_intensity* odd_fiber + even_fiber_intensity * even_fiber)


def trim_coefficients_to_fit_image(image, trace_fit_coefficients_no_indices):
    min_y, max_y = 0, image.data.shape[0]
    order_indices = np.array([i for i in range(0, trace_fit_coefficients_no_indices.shape[0])])
    trace_fit_coefficients = np.insert(trace_fit_coefficients_no_indices, obj=0, values=order_indices, axis=1)
    trace_values_versus_xpixel, num_traces, x = get_trace_centroids_from_coefficients(trace_fit_coefficients, image)
    good_indices = []
    for i in range(trace_values_versus_xpixel.shape[0]):
        if 1.1*np.mean(trace_values_versus_xpixel[i, :]) < max_y and (trace_values_versus_xpixel[i, :] > min_y).all():
            good_indices += [i]
    trimmed_trace_fit_coefficients_and_indices = trace_fit_coefficients[good_indices]
    assert (np.array(good_indices) - np.min(np.array(good_indices)) == np.array(list(range(len(good_indices))))).all()
    # replacing order indices with proper indicators
    trimmed_trace_fit_coefficients_and_indices[:, 0] = np.array(list(range(len(good_indices))))
    return trimmed_trace_fit_coefficients_and_indices


def generate_image_with_two_flat_traces(nx=1000, ny=50, readnoise=10, order_width=1.25, fiber_1_intensity=1E4, fiber_2_intensity=5E3, normalized_traces=False, add_noise=True):
    """
    :param nx:
    :param ny:
    :param readnoise:
    :param order_width:
    :param fiber_1_intensity: peak intensity of trace. Default 1E4 mimics a standard 120 sec NRES exposure.
    :param fiber_2_intensity: peak intensity of trace. Default 5E3 mimics a standard 120 sec NRES exposure.
    :param normalized_traces:
    :param add_noise:
    :return:
    """
    overscan_size = 2
    trimmed_shape = (ny, nx)
    image = FakeImage(nx=nx+overscan_size, ny=ny, overscan_size=overscan_size)
    trace_coefficients_no_indices = np.array([[image.data.shape[0]*1/3, 0, 0],
                                              [image.data.shape[0]*2/3, 0, 0]])

    image.trace.coefficients = trim_coefficients_to_fit_image(image, trace_coefficients_no_indices)
    if normalized_traces:
        gaussian_norm_factor = 1/np.sqrt(2 * np.pi * order_width ** 2)
        fill_image_with_traces(image, trimmed_shape=trimmed_shape, order_width=order_width,
                               odd_fiber_intensity=gaussian_norm_factor, even_fiber_intensity=gaussian_norm_factor)
    if not normalized_traces:
        fill_image_with_traces(image, trimmed_shape=trimmed_shape, order_width=order_width,
                               odd_fiber_intensity=fiber_2_intensity, even_fiber_intensity=fiber_1_intensity)
    image.readnoise = readnoise
    if add_noise:
        noisify_image(image, trimmed_shape=trimmed_shape)
    trim_image(image, trimmed_shape=trimmed_shape)
    return image
