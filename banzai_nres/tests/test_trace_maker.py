import pytest
import numpy as np
from scipy import ndimage
from unittest import mock
from astropy.io import fits

from banzai_nres.traces import InitialTraceFit, TraceMaker, LoadTrace
from banzai_nres.utils.trace_utils import get_coefficients_from_meta, generate_legendre_array, Trace
from banzai_nres.tests.utils import FakeImage, noisify_image, trim_image,\
                                    gaussian, generate_sample_astropy_nres_values_table
from banzai_nres.tests.adding_traces_to_images_utils import generate_image_with_two_flat_traces
from banzai_nres.tests.adding_traces_to_images_utils import trim_coefficients_to_fit_image, fill_image_with_traces
from banzai_nres.utils import trace_utils
from banzai_nres.images import NRESImage as Image
import banzai_nres.settings

from banzai.tests.utils import FakeContext
from banzai.images import DataTable


import logging

logger = logging.getLogger(__name__)


class FakeTraceImage(FakeImage):
    """
    Image of 500x500 is recommended. Drastic changes to that dimension may break trace testing.
    """
    def __init__(self, nx=500, ny=502, *args, **kwargs):
        super(FakeTraceImage, self).__init__(*args, **kwargs)
        self.caltype = 'TRACE'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'
        self.header['OBJECTS'] = 'tung&tung&none'
        self.nx = nx
        self.ny = ny
        self.bpm = np.zeros((self.ny, self.nx), dtype=np.uint8)
        self.data = np.zeros((self.ny, self.nx))


class TinyFakeImageWithTraces(object):
    """
    A tiny image with one fake trace running down the center (just a single row of bright pixels)
    """
    def __init__(self):
        size_of_test_image = 3
        self.data = np.zeros((size_of_test_image, size_of_test_image))
        self.data[1] = np.ones(size_of_test_image)

        self.x_pixel_coords = np.arange(size_of_test_image)
        self.legendre_polynomial_array = np.ones((1, size_of_test_image))
        self.legendre_polynomial_coefficients = np.array([1])
        self.negative_expected_sum_along_trace = -1 * np.sum(self.data[1])
        self.image_filt = ndimage.spline_filter(self.data)
        self.trace = Trace()


def munge_coefficients(even_coefficients, odd_coefficients):
    if even_coefficients.shape[0] != odd_coefficients.shape[0]:
        min_shape = min(odd_coefficients.shape[0], even_coefficients.shape[0]) - 1
        return even_coefficients[:min_shape], odd_coefficients[:min_shape]
    else:
        return even_coefficients, odd_coefficients


def make_random_yet_realistic_trace_coefficients(image, order_of_poly_fit=4):
    """
    :param image: Banzai_nres Image object
    Adds realistic coefficients for traces which fit entirely in the frame, saving into image.trace.coefficients
    and an arbitrary fiber_order onto image.trace.fiber_order
    """
    meta_coefficients_even = np.zeros((order_of_poly_fit + 1, 6))
    # NOTE: This 5 is poly_order + 1 We need polyorder as a global variable, or just never change it from 4.
    meta_coefficients_even[0] = [0, image.ny*20, 400, 5.90806434e+01, 4.94386504e+00, 1.37890482e+00]
    meta_coefficients_even[1] = [-3.64547386e+01, -5.01236304e+01, -1.65331378e+01, -3.31442330e+00, -7.46833391e-01, -8.14690916e-02]
    meta_coefficients_even[2] = [8.99994878e+01,  1.69576526e+00,  2.98550032e-01, -1.12856233e-01, 5.53028580e-03, -1.45839718e-01]
    meta_coefficients_even[3] = [-2.35116489e-01, -2.64776632e-01, -1.17453609e-01, -6.87267618e-02, -6.73017389e-02, -6.30241483e-02]
    meta_coefficients_even[4] = [5.80449884e-01, -3.74220905e-01,  1.84019236e-01,  5.20675962e-02, 1.23541898e-02,  2.08741426e-01]
    meta_coefficients_odd = np.copy(meta_coefficients_even)
    meta_coefficients_odd[0, 0] = 15

    for i in range(1, meta_coefficients_even.shape[0]):
        noise_scale = 1/100
        noise = np.random.normal(loc=0, scale=noise_scale, size=meta_coefficients_even.shape[1])
        meta_coefficients_even[i] += noise
        meta_coefficients_odd[i] += noise

    meta_legendre_array, x, xnorm = generate_legendre_array(image.data.shape[1], order_of_poly_fits=meta_coefficients_even.shape[1]-1)
    trace_coefficients_odd = get_coefficients_from_meta(meta_coefficients_odd, meta_legendre_array)
    trace_coefficients_even = get_coefficients_from_meta(meta_coefficients_even, meta_legendre_array)

    trace_coefficients_odd_and_indices = trim_coefficients_to_fit_image(image, trace_coefficients_odd)
    trace_coefficients_even_and_indices = trim_coefficients_to_fit_image(image, trace_coefficients_even)

    trace_coefficients_even_and_indices, trace_coefficients_odd_and_indices = munge_coefficients(
                                                trace_coefficients_even_and_indices, trace_coefficients_odd_and_indices)

    image.trace.fiber_order = (0, 1)
    image.trace.coefficients = np.vstack((trace_coefficients_even_and_indices, trace_coefficients_odd_and_indices))


def differences_between_found_and_generated_trace_vals(found_coefficients, image):
    """
    :param found_coefficients: Trace coefficients computed
    :param image: banzai image object with trace_Fit_coefficients not None
    :return: ndarray, the difference between the y values of the found traces and the old traces at each x value.
    """
    trace_values_1, num_traces_1, x = image.trace.get_trace_centroids_from_coefficients(image.data.shape[1],
                                                                                        coefficients_and_indices=found_coefficients)
    trace_values_2, num_traces_2, x = image.trace.get_trace_centroids_from_coefficients(image.data.shape[1])
    assert num_traces_1 == num_traces_2
    return trace_values_2 - trace_values_1


def array_with_two_peaks():
    """
    :return: generates a fake line cut down the ccd with traces present.
    """
    centroids = (15, 30)
    x = np.linspace(0, 50, num=100)
    vectorized_gaussian = np.vectorize(gaussian)
    y = vectorized_gaussian(x, 1, centroids[0], 1.5) + vectorized_gaussian(x, 1, centroids[1], 1.5)
    return y, x, centroids


def test_finding_first_statistically_significant_maxima():
    """
    test type: Unit Test.
    info: tests the function which generates an approximate initial guess for the location of the next trace for
    the blind trace maker. Also tests its ability to determine that no statistally significant maximum exists.
    """
    intensity_line_cut_across_two_traces, x_coords, centroids = array_with_two_peaks()
    approximate_maximum_flux = np.max(intensity_line_cut_across_two_traces)
    index_of_first_maximum, maximum_exists = trace_utils.maxima(intensity_line_cut_across_two_traces, 5, 1 / 20,
                                                                approximate_maximum_flux)
    assert np.isclose(x_coords[index_of_first_maximum], centroids[0], atol=3, rtol=0)
    index_of_first_maximum, maximum_exists = trace_utils.maxima(np.ones_like(x_coords), 5, 1 / 20,
                                                                approximate_maximum_flux)
    assert not maximum_exists


class TestUnitBlindFitAlgorithms:
    """
    Unit tests for various algorithms involved with blind-fitting via the order-by-order technique.
    """
    def test_generating_initial_guess_for_next_blind_fit_given_no_previous_fit(self):
        fake_image_data = np.zeros((10, 10))
        expected_coefficient_guess = [int(fake_image_data.shape[0]/3), 0, 90]
        coeff_guess, max_exists, refflux = trace_utils.generate_initial_guess_for_trace_polynomial(image_data=fake_image_data, x=None,
                                                                                                   evaluated_legendre_polynomials=None, order=2,
                                                                                                   second_order_coefficient_guess=expected_coefficient_guess[2],
                                                                                                   lastcoef=None, direction='up')
        assert max_exists
        assert coeff_guess == expected_coefficient_guess

    def test_generating_initial_guess_for_next_blind_fit(self):
        image = generate_image_with_two_flat_traces(nx=100, ny=100)

        start_locations = [1/2, 1/2, 1/3]
        expected_0th_order_next_guesses = [1/3, 2/3, 1/3]
        directions = ['down', 'up', 'inplace']
        for start_location, expected_guess, direction in zip(start_locations, expected_0th_order_next_guesses, directions):
            coefficients_of_last_fit = [int(image.data.shape[0] * start_location), 0]
            evaluated_legendre_polynomials, x_coords, xnorm = trace_utils.generate_legendre_array(image.data.shape[1],
                                                                                                  order_of_poly_fits=1)
            coeff_guess, max_exists, refflux = trace_utils.generate_initial_guess_for_trace_polynomial(image_data=image.data, x=x_coords,
                                                                                                       evaluated_legendre_polynomials=evaluated_legendre_polynomials,
                                                                                                       order=1, lastcoef=coefficients_of_last_fit,
                                                                                                       direction=direction)
            assert max_exists
            assert np.isclose(coeff_guess[0], expected_guess * image.data.shape[0], atol=5, rtol=0)

    def test_trimming_bad_fits_on_way_up(self):
        length = 12
        all_coefficients_and_indices_init = [[1, 3], [2, 8], [3, 11]]
        loop_counter = len(all_coefficients_and_indices_init) + 1

        coef = [length+1]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[4] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=True, done=False,
                                                                                                direction='up')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init
        assert num_orders_found == len(trimmed_all_coefficients)

        coef = [all_coefficients_and_indices_init[-1][1] - 3]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[4] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=True, done=False,
                                                                                                direction='up')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init
        assert num_orders_found == len(trimmed_all_coefficients)

        coef = [all_coefficients_and_indices_init[-1][1] - 3]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[4] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=False, done=False,
                                                                                                direction='up')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init
        assert num_orders_found == len(trimmed_all_coefficients)

    def test_trimming_bad_fits_on_way_down(self):
        length = 12
        all_coefficients_and_indices_init = [[1, 11], [2, 8], [3, 3]]
        loop_counter = len(all_coefficients_and_indices_init) + 1

        coef = [-2]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[4] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=True, done=False,
                                                                                                direction='down')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init
        assert num_orders_found == len(trimmed_all_coefficients)

        coef = [all_coefficients_and_indices_init[-1][1] + 3]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[4] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=True, done=False,
                                                                                                direction='down')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init
        assert num_orders_found == len(trimmed_all_coefficients)
        coef = [all_coefficients_and_indices_init[-1][1] - 3]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[4] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=False, done=False,
                                                                                                direction='down')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init
        assert num_orders_found == len(trimmed_all_coefficients)

    def test_removing_repeated_fits(self):
        length = 12
        all_coefficients_and_indices_init = [[1, 11], [2, 8], [3, 3], [4, 3]]
        loop_counter = len(all_coefficients_and_indices_init) + 1

        coef = [all_coefficients_and_indices_init[-1][1]]
        all_coefficients_and_indices = all_coefficients_and_indices_init + [[5] + coef]
        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=True, done=False,
                                                                                                direction='down')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init[:-1]
        assert num_orders_found == len(trimmed_all_coefficients)

        num_orders_found, trimmed_all_coefficients, done = trace_utils.validate_fit_and_trim_erroneous_fits(coef, all_coefficients_and_indices,
                                                                                                loop_counter, length,
                                                                                                maximum_exists=True, done=False,
                                                                                                direction='up')
        assert done
        assert trimmed_all_coefficients == all_coefficients_and_indices_init[:-1]
        assert num_orders_found == len(trimmed_all_coefficients)


def test_generating_blank_evaluated_legendre_array():
    tiny_image = TinyFakeImageWithTraces()
    x_coords = np.arange(tiny_image.data.shape[1])
    xnorm = x_coords * 2. / x_coords[-1] - 1
    legendre_polynomial_array, x_coords_2, xnorm_2 = trace_utils.generate_legendre_array(tiny_image.data.shape[1],
                                                                                         order_of_poly_fits=0)
    assert (legendre_polynomial_array[0] == 1).all()
    assert np.array_equal(x_coords_2, x_coords)
    assert np.array_equal(xnorm, xnorm_2)


class TestTraceClassMethods:
    def test_getting_trace_centroids_from_coefficients(self):
        tiny_image = TinyFakeImageWithTraces()
        tiny_image.trace.coefficients = np.array([[0, 1]])
        trace_values_versus_xpixel, num_traces, x_coord_array = \
            tiny_image.trace.get_trace_centroids_from_coefficients(tiny_image.data.shape[1])
        assert np.array_equal(x_coord_array, np.arange(tiny_image.data.shape[1]))
        assert num_traces == 1
        assert np.array_equal(trace_values_versus_xpixel, np.array([np.ones_like(x_coord_array)]))

    def test_converting_coefficients_array_to_astropy_table(self):
        fiber_orders_to_try = [None, (1, 2), (2, 1)]
        for fiber_order in fiber_orders_to_try:
            outputs = generate_sample_astropy_nres_values_table(fiber_order=fiber_order)
            coefficients_table = outputs[-1]
            assert coefficients_table['order'][0] == 0

    def test_converting_trace_y_values_array_array_to_astropy_table(self):
        fiber_orders_to_try = [None, (1, 2), (2, 1)]
        for fiber_order in fiber_orders_to_try:
            test_trace = Trace()
            indices = np.array([np.concatenate((np.arange(2), np.arange(2)))])
            coefficients = np.ones((4, 4))
            trace_centroids = np.ones((4, 10))
            coefficients_and_indices = np.hstack((indices.T, coefficients))
            trace_centroids_table = test_trace.convert_numpy_array_trace_centroids_to_astropy_table(num_lit_fibers=2,
                                                                                                    trace_centroids=trace_centroids,
                                                                                                    coefficients=coefficients_and_indices,
                                                                                                    fiber_order=fiber_order)
            assert trace_centroids_table['order'][0] == 0

    def test_converting_astropy_table_coefficients_to_array_and_fiber_order(self):
        fiber_orders_to_try = [None, (1, 2), (2, 1)]
        for fiber_order in fiber_orders_to_try:
            test_trace, coefficients_and_indices, coefficients_table = generate_sample_astropy_nres_values_table(
                                                                            fiber_order=fiber_order)
            load_coefficients, load_fiber_order = test_trace.convert_astropy_table_coefficients_to_numpy_array(
                                                                            coefficients_table)
            assert np.array_equal(load_coefficients, coefficients_and_indices)
            assert load_fiber_order == fiber_order

    def test_converting_astropy_table_y_values_array_and_fiber_order(self):
        fiber_orders_to_try = [None, (1, 2), (2, 1)]
        test_trace = Trace()
        for fiber_order in fiber_orders_to_try:
            test_trace, y_values_with_indices, y_values_table = generate_sample_astropy_nres_values_table(
                                                                            fiber_order=fiber_order,
                                                                            table_name=test_trace.trace_center_table_name)
            load_y_values, load_fiber_order = test_trace.convert_astropy_table_trace_y_values_to_numpy_array(
                                                                            y_values_table)
            assert np.array_equal(load_y_values, y_values_with_indices[:, 1:])
            assert load_fiber_order == fiber_order

    def test_converting_any_astropy_table_to_values_and_fiber_order(self):
        fiber_orders_to_try = [None, (1, 2), (2, 1)]
        for fiber_order in fiber_orders_to_try:
            test_trace, values_and_indices, astropy_table = generate_sample_astropy_nres_values_table(
                                                                            fiber_order=fiber_order)
            loaded_values, load_fiber_order = test_trace.recombine_values_from_table_into_nd_array_with_order_indices(
                                                                            astropy_table,
                                                                            name_of_values='coefficients')
            assert np.array_equal(loaded_values, values_and_indices)
            assert load_fiber_order == fiber_order


class TestLoadTrace:
    @mock.patch('banzai_nres.traces.LoadTrace.get_trace_coefficients')
    def test_load_trace_stage(self, mock_get_coefficients):
        mock_get_coefficients.return_value = 0
        images = [FakeTraceImage()]
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        trace_load_stage = LoadTrace(fake_context)
        images = trace_load_stage.do_stage(images)
        assert images[0].trace.coefficients == 0

    @mock.patch('banzai_nres.traces.os.path.exists')
    @mock.patch('banzai_nres.traces.fits.open')
    @mock.patch('banzai_nres.traces.dbs.get_master_calibration_image')
    def test_loading_coefficients_from_file(self, mock_cal, mock_fits_open, mock_os):
        """
        Tests that add_data_tables_to_hdu_list and regenerate_data_table_from_fits_hdu_list
        create fits.HDUList objects correctly from astropy tables with single element entries
        and for astropy tables with columns where each element is a list.
        """
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        fake_context.db_address = ''
        test_image = Image(fake_context, filename=None)
        test_image.filename = 'test.fits'
        trace_load_stage = LoadTrace(fake_context)
        trace_class = Trace()

        mock_os.return_value = True
        mock_cal.return_value = 'fake_cal.fits'

        table_name = trace_class.coefficients_table_name
        nn, coefficients_and_indices, coefficients_table = generate_sample_astropy_nres_values_table(fiber_order=None,
                                                                    table_name=trace_class.coefficients_table_name)
        test_image.data_tables[table_name] = DataTable(data_table=coefficients_table, name=table_name)
        hdu_list = []
        hdu_list = test_image._add_data_tables_to_hdu_list(hdu_list)
        fits_hdu_list = fits.HDUList(hdu_list)

        mock_fits_open.return_value = fits_hdu_list
        loaded_coefficients = trace_load_stage.get_trace_coefficients(test_image)
        assert (loaded_coefficients == coefficients_and_indices).all()

    @mock.patch('banzai_nres.traces.os.path.exists')
    @mock.patch('banzai_nres.traces.dbs.get_master_calibration_image')
    def test_loading_coefficient_failure(self, mock_cal, mock_os):
        """
        Tests that add_data_tables_to_hdu_list and regenerate_data_table_from_fits_hdu_list
        create fits.HDUList objects correctly from astropy tables with single element entries
        and for astropy tables with columns where each element is a list.
        """
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        fake_context.db_address = ''
        test_image = Image(fake_context, filename=None)
        test_image.filename = 'test.fits'

        trace_load_stage = LoadTrace(fake_context)
        mock_os.return_value = False
        mock_cal.return_value = 'fake_cal.fits'
        trace_load_stage.get_trace_coefficients(test_image)
        assert True


class TestFindingTotalFluxAcrossTraces:
    """
    test type: Unit Test.
    info: tests the three closely tied functions which evaluate
    the negative sum of the fluxes across a trace of given coefficients across the image.
    """
    def test_finding_flux_across_single_trace(self):
        tiny_image = TinyFakeImageWithTraces()
        found_value = trace_utils.negative_flux_across_trace(tiny_image.legendre_polynomial_coefficients, tiny_image.image_filt,
                                                             tiny_image.x_pixel_coords, tiny_image.legendre_polynomial_array)
        assert np.isclose(found_value, tiny_image.negative_expected_sum_along_trace)

    def test_finding_flux_across_single_trace_at_many_points(self):
        tiny_image = TinyFakeImageWithTraces()
        testpoints = np.array([1])
        found_value = (-1) * trace_utils.flux_across_trace_up_detector(testpoints, [], tiny_image.image_filt, tiny_image.x_pixel_coords,
                                                                       tiny_image.legendre_polynomial_array)[0]
        assert np.isclose(found_value, tiny_image.negative_expected_sum_along_trace)

    def test_finding_total_flux_without_filtered_image_data(self):
        tiny_image = TinyFakeImageWithTraces()
        coefficients_and_indices = np.array([[0, 1]])
        flux = trace_utils.totalflux_all_traces(coefficients_and_indices, tiny_image)
        assert np.isclose(-1 * flux, tiny_image.negative_expected_sum_along_trace)


def test_excluding_traces_which_are_cut_in_half_by_detector():
    fake_image = FakeImage(nx=12, ny=10)
    coefficients = np.array([[0, -5],
                             [1, 3]])
    trimmed_coeffs = trace_utils.exclude_traces_which_jet_off_detector(coefficients, fake_image.data)
    assert np.array_equal(trimmed_coeffs, coefficients[1:])


def test_first_sorting_of_coefficients_from_blind_fit():
    fake_coefficients_and_indices = np.array([[0, 1],
                                              [1, 3],
                                              [2, 2],
                                              [3, 4]])
    coefficients_and_indices = trace_utils.split_and_sort_coefficients_for_each_fiber(fake_coefficients_and_indices,
                                                                                      num_lit_fibers=2)
    expected_coefficients_and_indices = np.array([[0, 1],
                                                  [1, 2],
                                                  [0, 3],
                                                  [1, 4]])
    assert np.array_equal(coefficients_and_indices, expected_coefficients_and_indices)


def test_first_sorting_of_coefficients_from_blind_fit_with_odd_num_traces():
    fake_coefficients_and_indices = np.array([[0, 1],
                                              [1, 3],
                                              [2, 2],
                                              [3, 4],
                                              [4, 5]])
    coefficients_and_indices = trace_utils.split_and_sort_coefficients_for_each_fiber(fake_coefficients_and_indices,
                                                                                      num_lit_fibers=2)
    expected_coefficients_and_indices = np.array([[0, 1],
                                                  [1, 2],
                                                  [0, 3],
                                                  [1, 4]])
    assert np.array_equal(coefficients_and_indices, expected_coefficients_and_indices)


class TestTraceMaker:
    """
    Unit tests for TraceMaker class
    """
    def test_trace_maker_properties(self):
        trace_maker = TraceMaker(FakeContext(settings=banzai_nres.settings.NRESSettings()))
        assert trace_maker.calibration_type == 'TRACE'

    @mock.patch('banzai_nres.traces.NRESImage')
    @mock.patch('banzai_nres.traces.dbs.get_master_calibration_image')
    def test_trace_maker_does_not_crash_on_blank_frame(self, mock_cal, mock_images):
        readnoise = 11.0
        order_of_poly_fit = 4
        images = [FakeTraceImage(nx=100, ny=100)]
        images[0].readnoise = readnoise
        noisify_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        trim_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        fake_context.db_address = ''
        blind_trace_maker = InitialTraceFit(fake_context)
        blind_trace_maker.always_generate_traces_from_scratch = True
        mock_cal.return_value = None
        blind_trace_maker.order_of_poly_fit = order_of_poly_fit
        images = blind_trace_maker.do_stage(images)
        assert True

    @mock.patch('banzai_nres.traces.NRESImage')
    @mock.patch('banzai_nres.traces.dbs.get_master_calibration_image')
    def test_trace_maker(self, mock_cal, mock_images):
        """
        test type: Integration Test.
        info: This tests trace making via a blind fit.
        WARNING: Because trace fitting is defined with polynomials which are normalized from -1 to 1, if one squeezes
        the x axis of the image further, then the traces bend more drastically. Thus it is recommended you do not change the
        size of the FakeTraceImage.
        """
        readnoise = 11.0
        order_of_poly_fit = 4

        image = FakeTraceImage()
        image.readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(image, order_of_poly_fit=order_of_poly_fit)
        fill_image_with_traces(image, trimmed_shape=tuple([min(image.data.shape)] * 2))
        noisify_image(image, trimmed_shape=tuple([min(image.data.shape)] * 2))
        trim_image(image, trimmed_shape=tuple([min(image.data.shape)] * 2))

        original_coefficents = image.trace.coefficients

        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        fake_context.db_address = ''
        images = [image, image]

        blind_trace_maker = InitialTraceFit(fake_context)
        blind_trace_maker.max_number_of_images_to_fit = 1
        blind_trace_maker.order_of_poly_fit = order_of_poly_fit

        for force_traces_from_scratch, value in zip([True, False], [None, '']):
            blind_trace_maker.always_generate_traces_from_scratch = force_traces_from_scratch
            mock_cal.return_value = value
            images = blind_trace_maker.do_stage(images)

            master_cal_maker = TraceMaker(fake_context)
            master_cal_maker.do_stage(images)

            args, kwargs = mock_images.call_args
            trace_coefficients_data_table_name = Trace().coefficients_table_name
            master_trace_table = kwargs['data_tables'][trace_coefficients_data_table_name]._data_table

            coefficients_array, fiber_order = Trace().convert_astropy_table_coefficients_to_numpy_array(master_trace_table)
            logger.debug(coefficients_array.shape)
            images[0].trace.coefficients = coefficients_array

            difference = differences_between_found_and_generated_trace_vals(original_coefficents, images[0])
            logger.debug('median absolute deviation in unit-test trace fitting is {0} pixels'
                         .format(np.median(np.abs(difference - np.median(difference)))))
            logger.debug('standard deviation in unit-test trace fitting is {0} pixels'
                         .format(np.std(difference)))
            logger.debug('worst error (max of true minus found) in unit-test trace fitting is {0} pixels'
                         .format(np.max(np.abs(difference))))
            logger.debug('median error (median of true minus found) in unit-test trace fitting is {0} pixels'
                         .format(np.abs(np.median(difference))))

            assert np.median(np.abs(difference - np.median(difference))) < 2/100
            assert np.abs(np.median(difference)) < 2/100
