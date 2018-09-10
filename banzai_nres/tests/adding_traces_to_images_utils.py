from banzai_nres.utils.trace_utils import get_trace_centroids_from_coefficients
from banzai_nres.tests.utils import gaussian
import numpy as np


def fill_image_with_traces(image, trimmed_shape, order_width=1.25, odd_fiber_intensity=1E4, even_fiber_intensity=5E3):
    """
    :param image: Banzai_nres FakeImage object which is square after trimming.
    :param order_width : the standard deviation of an unnormalized guassian e.g. sigma
           from np.exp(-(x - b) ** 2 / (2 * sigma ** 2))
    fills the image.data with guassian profiled trace coefficients.
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
