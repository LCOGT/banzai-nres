import pytest
import mock
from banzai_nres.traces import GenerateInitialGuessForTraceFitFromScratch, TraceRefine, TraceSaver
from banzai.tests.utils import FakeContext

from scipy import ndimage
import numpy as np
from banzai import logs
from banzai_nres.utils.trace_utils import get_coefficients_from_meta, generate_legendre_array, get_trace_centroids_from_coefficients
from banzai_nres.tests.utils import FakeImage, noisify_image, trim_image, gaussian
from banzai_nres.tests.adding_traces_to_images_utils import generate_image_with_two_flat_traces
from banzai_nres.tests.adding_traces_to_images_utils import trim_coefficients_to_fit_image, fill_image_with_traces

from banzai_nres.utils import trace_utils
from astropy.io import fits


logger = logs.get_logger(__name__)


class FakeTraceImage(FakeImage):
    """
    Image must be square and at least 400x400 for trace making integration test to behave. 500x500 is recommended.
    """
    def __init__(self, *args, **kwargs):
        super(FakeTraceImage, self).__init__(*args, **kwargs)
        self.caltype = 'trace'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'
        self.header['OBJECTS'] = 'tung&tung&none'
        self.nx = 500
        self.ny = 502
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

    meta_legendre_array, x, xnorm = generate_legendre_array(image, order_of_poly_fits=meta_coefficients_even.shape[1]-1)
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
    trace_values_1, num_traces_1, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
    trace_values_2, num_traces_2, x = get_trace_centroids_from_coefficients(found_coefficients, image)
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
        coeff_guess, max_exists, refflux = trace_utils.generate_initial_guess_for_trace_polynomial(image=None,
                                                                imfilt=fake_image_data, x=None,
                                                                evaluated_legendre_polynomials=None, order=2,
                                                                second_order_coefficient_guess=expected_coefficient_guess[2],
                                                                lastcoef=None, direction='up')
        assert max_exists
        assert coeff_guess == expected_coefficient_guess

    def test_generating_initial_guess_for_next_blind_fit(self):
        image = generate_image_with_two_flat_traces(nx=100, ny=100)

        image_data_filtered = ndimage.spline_filter(image.data)

        start_locations = [1/2, 1/2, 1/3]
        expected_0th_order_next_guesses = [1/3, 2/3, 1/3]
        directions = ['down', 'up', 'inplace']
        for start_location, expected_guess, direction in zip(start_locations, expected_0th_order_next_guesses, directions):
            coefficients_of_last_fit = [int(image.data.shape[0] * start_location), 0]
            evaluated_legendre_polynomials, x_coords, xnorm = trace_utils.generate_legendre_array(image, order_of_poly_fits=1)
            coeff_guess, max_exists, refflux = trace_utils.generate_initial_guess_for_trace_polynomial(image=image,
                                                                    imfilt=image_data_filtered, x=x_coords,
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
    legendre_polynomial_array, x_coords_2, xnorm_2 = trace_utils.generate_legendre_array(tiny_image,
                                                                                         order_of_poly_fits=0)
    assert (legendre_polynomial_array[0] == 1).all()
    assert np.array_equal(x_coords_2, x_coords)
    assert np.array_equal(xnorm, xnorm_2)


def test_checking_for_close_fit_between_two_fits():
    image = FakeImage(nx=102, ny=100, overscan_size=2)
    trim_image(image, trimmed_shape=(100, 100))
    indices = np.array([np.concatenate((np.arange(6), np.arange(6)))])
    coefficients = np.ones_like(indices.T)*50
    coefficients_and_indices = np.hstack((indices.T, coefficients))
    images = [image, image]
    a_close_fit = trace_utils.check_for_close_fit([coefficients_and_indices, coefficients_and_indices],
                                                  images, num_lit_fibers=2, max_pixel_error=1E-1)
    assert a_close_fit
    coefficients_and_indices_new = np.copy(coefficients_and_indices)
    coefficients_and_indices_new[:, 1] += 5
    a_close_fit = trace_utils.check_for_close_fit([coefficients_and_indices_new, coefficients_and_indices],
                                                  images, num_lit_fibers=2, max_pixel_error=1E-1)
    assert not a_close_fit


def test_checking_for_flux_change_between_two_fits():
    tiny_image = TinyFakeImageWithTraces()
    coefficients_and_indices_init = np.array([[0, 1]])
    no_flux_change = trace_utils.check_flux_change(coefficients_and_indices_init, coefficients_and_indices_init,
                                  tiny_image, relative_tolerance=1E-2)
    assert no_flux_change
    coefficients_and_indices_new = np.array([[0, 0]])
    no_flux_change = trace_utils.check_flux_change(coefficients_and_indices_new, coefficients_and_indices_init,
                                  tiny_image, relative_tolerance=1E-2)
    assert not no_flux_change


def test_getting_trace_centroids_from_coefficients():
    tiny_image = TinyFakeImageWithTraces()
    coefficients_and_indices = np.array([[0, 1]])
    trace_values_versus_xpixel, num_traces, x_coord_array = \
        trace_utils.get_trace_centroids_from_coefficients(coefficients_and_indices, tiny_image)
    assert np.array_equal(x_coord_array, np.arange(tiny_image.data.shape[1]))
    assert num_traces == 1
    assert np.array_equal(trace_values_versus_xpixel, np.array([np.ones_like(x_coord_array)]))


class TestFindingTotalFluxAcrossTraces:
    """
    test type: Unit Test.
    info: tests the three closely tied functions which evaluate
    the negative sum of the fluxes across a trace of given coefficients across the image.
    """
    def test_finding_flux_across_single_trace(self):
        tiny_image = TinyFakeImageWithTraces()
        found_value = trace_utils.crosscoef(tiny_image.legendre_polynomial_coefficients, tiny_image.image_filt,
                                            tiny_image.x_pixel_coords, tiny_image.legendre_polynomial_array)
        assert np.isclose(found_value, tiny_image.negative_expected_sum_along_trace)

    def test_finding_flux_across_single_trace_at_many_points(self):
        tiny_image = TinyFakeImageWithTraces()
        testpoints = np.array([1])
        found_value = (-1) * trace_utils.fluxvalues(testpoints, [], tiny_image.image_filt, tiny_image.x_pixel_coords,
                                                    tiny_image.legendre_polynomial_array)[0]
        assert np.isclose(found_value, tiny_image.negative_expected_sum_along_trace)

    def test_finding_total_flux_without_filtered_image_data(self):
        tiny_image = TinyFakeImageWithTraces()
        coefficients_and_indices = np.array([[0, 1]])
        flux = trace_utils.totalflux_all_traces(coefficients_and_indices, tiny_image)
        assert np.isclose(-1 * flux, tiny_image.negative_expected_sum_along_trace)


class TestImageSplinesClassMethods:
    """
    Tests the methods which generate the spline 0th, 1st and second derivatives inside the class
    trace_utils.ImageSplines
    """
    def test_generates_with_none_bpm(self):
        image = FakeImage(nx=12, ny=10, overscan_size=2)
        trim_image(image, trimmed_shape=(10, 10))
        image.bpm = None
        image_splines = trace_utils.ImageSplines()
        image_splines.calculate_spline_derivatives_and_populate_attributes(image, image.bpm)
        for attribute in ['spline', 'first_derivative', 'second_derivative']:
            assert getattr(image_splines, attribute) is not None

    def test_generates_with_zeros_bpm(self):
        image = FakeImage(nx=12, ny=10, overscan_size=2)
        trim_image(image, trimmed_shape=(10, 10))
        image_splines = trace_utils.ImageSplines()
        image_splines.calculate_spline_derivatives_and_populate_attributes(image, image.bpm)
        for attribute in ['spline', 'first_derivative', 'second_derivative']:
            assert getattr(image_splines, attribute) is not None

    def test_spline_skips_bad_pixel(self):
        image1 = FakeImage(nx=12, ny=10, overscan_size=2)
        image2 = FakeImage(nx=12, ny=10, overscan_size=2)
        trim_image(image1, trimmed_shape=(10, 10))
        trim_image(image2, trimmed_shape=(10, 10))
        image_splines1 = trace_utils.ImageSplines()
        image_splines2 = trace_utils.ImageSplines()

        # creating a bad pixel in one of the clone images
        image1.data[4, 4] = 10000
        image1.bpm[4, 4] = 1

        image_splines1.calculate_spline_derivatives_and_populate_attributes(image1, image1.bpm)
        image_splines2.calculate_spline_derivatives_and_populate_attributes(image2, image2.bpm)
        for attribute in ['spline', 'first_derivative', 'second_derivative']:
            assert np.isclose(getattr(image_splines1, attribute)[4](4), getattr(image_splines2, attribute)[4](4))


class TestTransformationsToAndFromMetaCoeffstoTraceCoeffs:
    """
    Tests initial fits of trace coeffs which generates the meta coefficients and the transformation from meta
    coefficients back to trace coefficients.
    """
    def test_generating_meta_coefficients(self):
        order_indices = np.array([0, 1])
        trace_coefficients = np.array([[1],
                                       [1]])
        meta_fit_coefficients = trace_utils.fit_trace_coeffs_to_generate_meta_coeffs(order_indices, trace_coefficients,
                                                                                     metapolyorder=0)
        assert np.allclose(meta_fit_coefficients, np.array([[1]]))

    def test_transforming_from_meta_to_trace_coeffs(self):
        trace_coeffs = trace_utils.get_coefficients_from_meta(allmetacoeffs=np.array([[1]]),
                                                              stpolyarr=np.array([[1, 1]]))
        expected_trace_coeffs = np.array([[1],
                                          [1]])
        assert np.allclose(trace_coeffs, expected_trace_coeffs)

    def test_generating_legendre_polynomial_for_generating_meta_coefficients(self):
        order_indices = np.array([0, 1])
        order_norm = order_indices * 2. / order_indices[-1] - 1
        meta_coefficients = (1,)
        trace_coefficients_this_order = trace_utils.legpolynomial(order_norm, *meta_coefficients)
        assert np.allclose(trace_coefficients_this_order, np.array([1,1]))


class TestMetaHessianandMetaGradientandTotalFluxEvaluation:
    """
    Testing Hessian and gradient creation for the meta fit and reorganization of the hessian into a matrix.
    This class tests all of the main functions needed to perform the global meta fit on traces.
    This tests whether the elements which fill the hessian have the correct ordering, these are the p,j,q,k elements
    ordered in the same way as here: https://v2.overleaf.com/read/jtckthqsdttj
    The same file can be found in docs/algorithm_docs/Newton_s_method_and_meta_fits.pdf
    Note: np.equal(a,b) is depriciated for element wise string comparison.
    """
    def constant_function(self, x):
        return np.ones_like(x)

    def dummy_meta_hessian_element(self, p, q, j, k, *extraargs):
        return '{0}, {1}, {2}, {3}'.format(p, j, q, k)

    def dummy_meta_gradient_element(self, p, k, *extraargs):
        return '{0}, {1}'.format(p, k)

    def test_generating_meta_hessian(self):

        meta_hessian_elements = trace_utils.evaluate_list_of_elements_of_hessian(stpolyarr=None,
                                                                                 array_of_individual_hessians=None,
                                                                                 tracepolyorder=1, metapolyorder=0,
                                                                                 element_generating_function=
                                                                                 self.dummy_meta_hessian_element)
        meta_hessian = trace_utils.reshape_hessian_elements_into_twod_matrix(meta_hessian_elements,
                                                                             tracepolyorder=1, metapolyorder=0)
        correct_hessian = np.array([['0, 0, 0, 0', '0, 0, 0, 1'],
                                    ['0, 1, 0, 0', '0, 1, 0, 1']])
        equally_correct_hessian = correct_hessian.T
        assert (meta_hessian == correct_hessian).all() or (meta_hessian == equally_correct_hessian).all()

    def test_generating_meta_hessian_elements(self):
        list_of_hessians = np.ones((2, 2, 2))
        evaluated_polynomials_for_meta = np.ones((2, 2))
        assert np.isclose(trace_utils.p_q_j_k_element_of_meta_hessian(0, 0, 0, 0, evaluated_polynomials_for_meta,
                                                                      list_of_hessians), 2)

    def test_generating_gradient(self):
        meta_gradient = trace_utils.evaluate_meta_gradient(stpolyarr=None, array_of_individual_gradients=None,
                                                           tracepolyorder=1, metapolyorder=0,
                                                           element_generating_function=self.dummy_meta_gradient_element)
        correct_gradient = np.array(['0, 0', '0, 1'])
        assert (correct_gradient == meta_gradient).all()

    def test_generating_meta_gradient_elements(self):
        array_of_gradients = np.ones((2, 2))
        evaluated_polynomials_for_meta = np.ones((2, 2))
        assert np.isclose(trace_utils.p_k_element_of_meta_gradient(0, 0, evaluated_polynomials_for_meta, array_of_gradients), 2)

    def test_end_to_end_building_of_meta_hessian_and_meta_gradient(self):
        meta_fit_coefficients = np.array([[1],
                                          [1]])
        coeffs_vector = meta_fit_coefficients.reshape(meta_fit_coefficients.size)

        image = FakeImage(nx=12, ny=10, overscan_size=2)
        trim_image(image, trimmed_shape=(10, 10))
        constant_spline_functions = [self.constant_function] * image.data.shape[1]
        image_splines = trace_utils.ImageSplines(spline=constant_spline_functions,
                                                 first_derivative=constant_spline_functions,
                                                 second_derivative=constant_spline_functions)
        x_coord_array = np.arange(image.data.shape[1])
        evaluated_polynomials_for_meta = np.array([[1, 1]])
        evaluated_polynomials_for_traces = np.array([np.ones(image.data.shape[1]),
                                                     np.ones(image.data.shape[1])])

        meta_fit_extra_args = (image_splines, evaluated_polynomials_for_meta,
                               evaluated_polynomials_for_traces, x_coord_array)

        meta_hessian = (-1) * trace_utils.NegativeHessian(coeffs_vector, *meta_fit_extra_args)
        assert np.isclose(meta_hessian, image.data.shape[0] * evaluated_polynomials_for_meta.shape[1]).all()

        meta_gradient = (-1) * trace_utils.NegativeGradient(coeffs_vector, *meta_fit_extra_args)
        assert np.isclose(meta_gradient, image.data.shape[0] * evaluated_polynomials_for_meta.shape[1]).all()

    def test_finding_total_flux_from_meta_coefficients(self):
        meta_fit_coefficients = np.array([[1]])
        coeffs_vector = meta_fit_coefficients.reshape(meta_fit_coefficients.size)
        image = FakeImage(nx=12, ny=10, overscan_size=2)
        trim_image(image, trimmed_shape=(10, 10))
        image_splines = trace_utils.ImageSplines()
        image_splines.calculate_spline_derivatives_and_populate_attributes(image, image.bpm)
        evaluated_polynomials_for_meta = np.array([[1, 1]])
        evaluated_polynomials_for_traces = np.array([np.ones(image.data.shape[1])])
        x_coord_array = np.arange(image.data.shape[1])

        total_positive_flux = -1 * trace_utils.neg_totalflux_for_scipy(coeffs_vector, image_splines,
                                                                       evaluated_polynomials_for_meta,
                                                                       evaluated_polynomials_for_traces, x_coord_array)
        assert np.isclose(2 * np.sum(evaluated_polynomials_for_traces), total_positive_flux)


class TestMakingPairsofLampflatstoFit:
    """
    test type: Unit Test.
    info: tests the cross correlate images function for meta fits.
    """
    def test_making_list_of_pairs_of_lampflats_to_fit(self):
        list_of_images = [0, 1, 2]
        indices_to_try, try_combinations_of_images = trace_utils.cross_correlate_image_indices(list_of_images,
                                                                                               cross_correlate_num=2)
        assert try_combinations_of_images
        for pair in [(0, 1), (1, 2), (0, 2)]:
            assert pair in indices_to_try
        assert len(indices_to_try) == 3

    def test_zero_cross_correlate_returns_entry_list(self):
        list_of_images = [0, 1, 2]
        indices_to_try, try_combinations_of_images = trace_utils.cross_correlate_image_indices(list_of_images,
                                                                                               cross_correlate_num=0)
        assert indices_to_try == list_of_images
        assert not try_combinations_of_images

    def test_empty_list_returns_empty_list(self):
        indices_to_try, try_combinations_of_images = trace_utils.cross_correlate_image_indices(images=[],
                                                                                               cross_correlate_num=2)
        assert indices_to_try == []
        assert not try_combinations_of_images


def test_excluding_traces_which_are_cut_in_half_by_detector():
    fake_image = FakeImage(nx=12, ny=10)
    coefficients = np.array([[0, -5],
                             [1, 3]])
    trimmed_coeffs = trace_utils.exclude_traces_which_jet_off_detector(coefficients, fake_image)
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


def test_splitting_coefficients_per_fiber():
    fake_coefficients_and_indices = np.array([[0, 1],
                                              [1, 2],
                                              [0, 3],
                                              [1, 4]])
    per_fiber_coeffs = trace_utils.split_already_sorted_coefficients_into_each_fiber(fake_coefficients_and_indices,
                                                                                     num_lit_fibers=2)
    assert np.array_equal(per_fiber_coeffs[0], fake_coefficients_and_indices[:2])
    assert np.array_equal(per_fiber_coeffs[1], fake_coefficients_and_indices[2:])


def test_getting_number_of_lit_fibers():
    image = FakeImage(nx=7, ny=5, overscan_size=2)
    real_num_lit_fibers_list = [0, 1, 2, 3]
    objects_list = ['none&none&none', 'none&tung&none', 'tung&tung&none', 'tung&tung&tung']
    for objects_string, real_num_lit_fibers in zip(objects_list, real_num_lit_fibers_list):
        image.header = {'OBJECTS': objects_string}
        num_lit_fibers = trace_utils.get_number_of_lit_fibers(image)
        assert num_lit_fibers == real_num_lit_fibers


class TestTraceRefine:
    """
    Unit tests for trace refine
    """
    @mock.patch('banzai_nres.traces.Image')
    def test_trace_refine_converges(self, mock_images):
        """
        test type: Unit Test. tests that TraceRefine stage converges to true traces.
        """
        readnoise = 11.0
        order_of_meta_fit = 6
        order_of_poly_fit = 4

        images = [FakeTraceImage()]
        images[0].readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(images[0], order_of_poly_fit=order_of_poly_fit)
        fill_image_with_traces(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        noisify_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        trim_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))

        original_coefficents = images[0].trace.coefficients

        trace_refiner = TraceRefine(FakeContext())
        trace_refiner.order_of_meta_fit = order_of_meta_fit
        images = trace_refiner.do_stage(images)

        difference = differences_between_found_and_generated_trace_vals(original_coefficents, images[0])
        logger.debug('error in unit-test trace fitting is less than {0} of a pixel'
                     .format(np.median(np.abs(difference - np.median(difference)))))
        logger.debug('worst error in unit-test trace fitting is {0} pixels'.format(np.max(np.abs(difference))))
        logger.debug('systematic error (median difference) in unit-test trace fitting is less than {0} of a pixel'
                     .format(np.abs(np.median(difference))))

        assert np.median(np.abs(difference - np.median(difference))) < 2/100
        assert np.abs(np.median(difference)) < 2/100

    def test_trace_refine_refits_if_traces_shifted_too_much(self):
        readnoise = 11.0
        order_of_meta_fit = 6
        order_of_poly_fit = 4

        image = FakeTraceImage()
        image.readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(image, order_of_poly_fit=order_of_poly_fit)
        fill_image_with_traces(image, trimmed_shape=tuple([min(image.data.shape)] * 2))
        noisify_image(image, trimmed_shape=tuple([min(image.data.shape)] * 2))
        trim_image(image, trimmed_shape=tuple([min(image.data.shape)] * 2))

        shifted_coefficents = np.copy(image.trace.coefficients)
        shifted_coefficents[:, 1] -= 3/4  # shifting 'found traces' by 3/4 of a pixel.

        trace_refiner = TraceRefine(FakeContext())
        trace_refiner.order_of_meta_fit = order_of_meta_fit
        fiber_order, refit_coefficients = trace_refiner.refit_if_necessary(image, num_lit_fibers=2,
                                                refined_trace_coefficients=shifted_coefficents,
                                                absolute_pixel_tolerance=1/2)

        difference = differences_between_found_and_generated_trace_vals(refit_coefficients, image)

        assert fiber_order is None
        assert np.median(np.abs(difference - np.median(difference))) < 2/100
        assert np.abs(np.median(difference)) < 2/100


@mock.patch('banzai_nres.traces.Image')
def test_blind_trace_maker(mock_images):
    """
    test type: Integration Test.
    info: This tests blind trace making (which involves blind trace making and then refining via meta-fit).
    Note: The fake images made here must have enough orders N such that N > order_of_meta_fit + 1.
    Currently it creates ~20 traces so 10 orders. A bigger image would have more traces, but expanding the
    image size may cause the created traces to be unrealistic. My suggestion is to never change the fake_image size
    and to keep order_of_meta_fit less than 8. (6 works well so why would you want it larger?)
    WARNING: Because trace fitting is defined with polynomials which are normalized from -1 to 1, if one squeezes
    the x axis of the image further, then the traces bend more drastically. Thus it is recommended you do not change the
    size of the FakeTraceImage.
    """
    num_trials = 2
    readnoise = 11.0
    order_of_poly_fit = 4
    order_of_meta_fit = 6

    for x in range(num_trials):
        images = [FakeTraceImage()]
        images[0].readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(images[0], order_of_poly_fit=order_of_poly_fit)
        fill_image_with_traces(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        noisify_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        trim_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))

        original_coefficents = images[0].trace.coefficients

        blind_trace_maker = GenerateInitialGuessForTraceFitFromScratch(FakeContext())
        blind_trace_maker.order_of_poly_fit = order_of_poly_fit
        images = blind_trace_maker.do_stage(images)

        trace_refiner = TraceRefine(FakeContext())
        trace_refiner.order_of_meta_fit = order_of_meta_fit
        images = trace_refiner.do_stage(images)
        master_cal_maker = TraceSaver(FakeContext())
        master_cal_maker.do_stage(images)

        args, kwargs = mock_images.call_args
        master_trace = kwargs['data']
        logger.debug(master_trace.shape)
        images[0].trace.coefficients = master_trace

        difference = differences_between_found_and_generated_trace_vals(original_coefficents, images[0])
        logger.debug('error in unit-test trace fitting is less than {0} of a pixel'
                     .format(np.median(np.abs(difference - np.median(difference)))))
        logger.debug('worst error in unit-test trace fitting is {0} pixels'.format(np.max(np.abs(difference))))
        logger.debug('systematic error (median difference) in unit-test trace fitting is less than {0} of a pixel'
                     .format(np.abs(np.median(difference))))

        assert np.median(np.abs(difference - np.median(difference))) < 2/100
        assert np.abs(np.median(difference)) < 2/100
