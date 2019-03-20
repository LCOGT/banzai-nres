import numpy as np

from banzai_nres.tests.test_traces import FakeTraceImage
from banzai_nres.tests.utils import fill_image_with_traces

from banzai_nres.utils.trace_utils import Trace
from banzai_nres.utils import extract_utils

import matplotlib.pyplot as plt


def test_rectify_orders():
    image = FakeTraceImage()
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                  max_num_traces=1)
    trace = Trace(data={'id': np.arange(trace_centers.shape[0]), 'centers': trace_centers})

    assert True


def test_rectify_curved_order_maps_all_values():
    image = FakeTraceImage()
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                  max_num_traces=1)
    x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
    image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
    single_order_centers = trace_centers[0]
    rectified_order, image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                              single_order_centers, half_window=10,
                                                              null_mapped_values=True)
    assert np.allclose(image_data, 0)


def test_rectify_flat_order():
    hw = 10
    image = FakeTraceImage()
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=1,
                                                                                  max_num_traces=1)
    x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
    image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
    single_order_centers = trace_centers[0]
    rectified_order, image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                              single_order_centers, half_window=hw,
                                                              null_mapped_values=False)
    trace_y_value = int(trace_centers[0][0])
    assert np.allclose(rectified_order, image_data[trace_y_value - hw: trace_y_value + hw + 1, :])
