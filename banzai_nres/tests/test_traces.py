# # Code for simulating traces on an iamge:
# import matplotlib.pyplot as plt
# from xwavecal.tests.utils import fill_image_with_traces
# import numpy as np
# image = type('', (), {'data':np.zeros((500, 500))})
# image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image)
# plt.imshow(image.data)
# plt.show()

import numpy as np
from banzai_nres.traces import find_y_center


def gaussian(x, a, b, sigma):
    return a * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def fill_image_with_traces(image, poly_order_of_traces=4, order_width=1.5, fiber_intensity=1E4, max_num_traces=1000):
    """
    :param image: Banzai image object, where image.data is a 2d array of the image data
    :param poly_order_of_traces: the max order of the polynomial describing the traces. Maximum of 4.
    :param order_width: width of the traces in pixels
    :param fiber_intensity: peak intensity of the orders
    :param max_num_traces: max number of traces to try and fit onto the image
    :return: An image populated with semi-realistic traces which bend like parabolas, but are of degree
    min(4, poly_order_of_traces).
    """

    num_fake_traces = min(int((image.data.shape[1] - 60)/20), max_num_traces)
    coefficients = np.zeros((num_fake_traces, 4+1))
    coefficients[:, 0] = np.linspace(30, image.data.shape[0] - 30, num=num_fake_traces)
    if poly_order_of_traces >= 2:
        coefficients[:, 2] = np.linspace(30, 40, num=num_fake_traces)
    if poly_order_of_traces >= 3:
        coefficients[:, 3] = np.linspace(1, 3, num=num_fake_traces)
    if poly_order_of_traces >= 4:
        coefficients[:, 4] = np.linspace(5, 10, num=num_fake_traces)

    trace_centers = trace_fitter._centers_from_coefficients(coefficients)
    trace_overlay = np.zeros_like(image.data).astype(np.float64)
    vectorized_gaussian = np.vectorize(gaussian)
    for x_pixel in range(trace_centers.shape[1]):
        for i in range(num_fake_traces):
            centroid = trace_centers[i, x_pixel]
            low, high = max(0, int(centroid - 5 * order_width)), min(trace_centers.shape[1] - 1,
                                                                     int(centroid + 5 * order_width)) + 1
            evalwindow = np.arange(low, high, 1)
            if len(evalwindow) > 0:
                trace_overlay[low: high, x_pixel] += vectorized_gaussian(evalwindow, 1, centroid, order_width)
    image.data += trace_overlay * fiber_intensity
    second_order_coefficient_guess = np.mean(coefficients[:, 2])
    return image, trace_centers, second_order_coefficient_guess


def test_centroid_int_weights():
    nx, ny = 103, 107
    fake_data = np.zeros((ny, nx))
    fake_data[49] = 1
    fake_data[50] = 1
    fake_data[51] = 1
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    center_region = np.logical_and(Y > 10, Y < (ny - 10))
    center = find_y_center(Y, center_region, fake_data)
    np.testing.assert_allclose(center, 50.0, atol=0.1)


def test_centroid_flux_weights():
    nx, ny = 103, 107
    fake_data = np.zeros((ny, nx))
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    center_region = np.logical_and(Y > 10, Y < (ny - 10))
    fake_data[center_region] += gaussian(Y[center_region], 10.0, 45, 3.0)
    center = find_y_center(Y, center_region, fake_data)
    np.testing.assert_allclose(center, 45.0, atol=0.1)


def test_blind_solve():
    pass


def test_blind_solve_realistic_data():
    pass


def test_refining_on_noisy_data():
    pass


def test_refine_traces_blind_solve():
    pass


def test_refine_traces_offset_centroid():
    pass

