import numpy as np

from banzai_nres.tests.test_traces import FakeTraceImage
from banzai_nres.tests.utils import fill_image_with_traces

from banzai_nres.utils.trace_utils import Trace
from banzai_nres.utils import extract_utils

from banzai_nres.extract import BoxExtract

import matplotlib.pyplot as plt


def test_rectify_orders():
    image = FakeTraceImage()
    hw = 10
    peak_intensity = 1E4
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                  max_num_traces=10,
                                                                                  fiber_intensity=peak_intensity)
    trace = Trace(data={'id': np.arange(trace_centers.shape[0]), 'centers': trace_centers})
    rectified_orders, zeroed_image_data = extract_utils.rectify_orders(image.data, trace,
                                                                       half_window=hw,
                                                                       debug=True)
    assert np.allclose(zeroed_image_data, 0)
    assert not np.allclose(image.data, 0)
    for key, item in rectified_orders.items():
        assert np.isclose(np.median(item[hw]), peak_intensity, rtol=0.02)


def test_rectify_curved_order_maps_all_values():
    image = FakeTraceImage()
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                  max_num_traces=1)
    x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
    image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
    single_order_centers = trace_centers[0]
    rectified_order, zeroed_image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                                     single_order_centers, half_window=10,
                                                                     nullify_mapped_values=True)
    assert np.allclose(zeroed_image_data, 0)


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
                                                              nullify_mapped_values=False)
    trace_y_value = int(trace_centers[0][0])
    assert np.allclose(rectified_order, image_data[trace_y_value - hw: trace_y_value + hw + 1, :])


def test_box_extract_accuracy():
    image = FakeTraceImage()
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                  max_num_traces=1)
    spectrum = BoxExtract().extract_order(image.data)
    assert np.allclose(spectrum / np.median(spectrum), 1)


def test_rectification_does_not_change_box_extract():
    image = FakeTraceImage()
    image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                  max_num_traces=1)
    x_coordinates, y_coordinates = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
    image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
    single_order_centers = trace_centers[0]
    rectified_order, image_data = extract_utils.rectify_order(image.data, image_coordinates,
                                                              single_order_centers, half_window=10,
                                                              nullify_mapped_values=False)
    rectified_spectrum = BoxExtract().extract_order(rectified_order)
    spectrum = BoxExtract().extract_order(image_data)
    assert np.allclose(spectrum / np.median(spectrum), 1)
    assert np.allclose(rectified_spectrum, spectrum)
