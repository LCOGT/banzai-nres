import numpy as np

from banzai_nres import settings as nres_settings
from banzai_nres.tests.test_traces import FakeTraceImage
from banzai_nres.tests.utils import fill_image_with_traces
from banzai_nres.utils.trace_utils import Trace
from banzai_nres.utils import extract_utils
from banzai_nres.extract import BoxExtract, BoxExtractor, RectifyTwodSpectrum

from banzai.tests.utils import FakeContext


class TestRectify:
    def test_rectify_orders(self):
        image = FakeTraceImage()
        hw = 10
        peak_intensity = 1E4
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=3,
                                                                                      fiber_intensity=peak_intensity)
        trace = Trace(data={'id': np.arange(trace_centers.shape[0]), 'centers': trace_centers})
        rectified_orders, zeroed_image_data = extract_utils.rectify_orders(image.data, trace,
                                                                           half_window=hw,
                                                                           debug=True)
        assert np.allclose(zeroed_image_data, 0)
        assert not np.allclose(image.data, 0)
        for key, item in rectified_orders.items():
            assert np.isclose(np.median(item[hw]), peak_intensity, rtol=0.02)

    def test_rectify_curved_order_maps_all_values(self):
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

    def test_rectify_flat_order(self):
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

    def test_rectification_does_not_change_box_extract(self):
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


class TestBoxExtract:
    def test_box_extract_accuracy(self):
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=1)
        spectrum = BoxExtract().extract_order(image.data)
        assert np.allclose(spectrum / np.median(spectrum), 1)

    def test_box_extract_trims_rectified_data(self):
        fake_context = FakeContext()
        max_extract_window = 10
        for half_window in [2, 6, 10, 15]:
            fake_spectrum = np.zeros((2 * max_extract_window + 1, 5))
            fake_spectrum[max_extract_window] = 1
            extractor = BoxExtractor(fake_context)
            extractor.extraction_half_window = half_window
            extractor.max_extraction_half_window = max_extract_window
            trimmed_spectrum = extractor._trim_rectified_data(rectified_twod_spectrum={'1': fake_spectrum})
            hw = min(half_window, max_extract_window)
            assert np.isclose(trimmed_spectrum['1'].shape[0], 2*hw + 1)
            assert np.allclose(trimmed_spectrum['1'][hw], 1)

    def test_integration_box_extract(self):
        fake_context = FakeContext()
        image = FakeTraceImage()
        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image, poly_order_of_traces=4,
                                                                                      max_num_traces=2,
                                                                                      fiber_intensity=1E4)
        image.trace = Trace(data={'id': np.arange(trace_centers.shape[0]), 'centers': trace_centers})
        image = RectifyTwodSpectrum(fake_context).do_stage(image)
        image = BoxExtractor(fake_context).do_stage(image)
        for spectrum in image.data_tables['box_extracted_spectrum']['flux']:
            assert np.median(spectrum) > 1E4
