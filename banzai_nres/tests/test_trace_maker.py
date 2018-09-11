import pytest
import mock
from banzai_nres.traces import BlindTraceMaker
from banzai.tests.utils import FakeContext

from scipy import ndimage
import numpy as np
from banzai import logs
from banzai_nres.utils.trace_utils import get_coefficients_from_meta, generate_legendre_array, get_trace_centroids_from_coefficients
from banzai_nres.tests.utils import FakeImage, noisify_image, trim_image, gaussian
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
        self.nx = 500
        self.ny = 502
        self.bpm = np.zeros((self.ny, self.nx), dtype=np.uint8)
        self.data = np.zeros((self.ny, self.nx))


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
    the blind trace maker.
    """
    intensity_line_cut_across_two_traces, x_coords, centroids = array_with_two_peaks()
    approximate_maximum_flux = np.max(intensity_line_cut_across_two_traces)
    index_of_first_maximum, maximum_exists = trace_utils.maxima(intensity_line_cut_across_two_traces, 5, 1 / 20,
                                                                approximate_maximum_flux)[0]
    assert np.isclose(x_coords[index_of_first_maximum], centroids[0], atol=3, rtol=0)


def test_finding_total_flux_across_a_trace():
    """
    test type: Unit Test.
    info: tests the two closely tied functions which evaluate
    the negative sum of the fluxes across a trace of given coefficients across the image.
    """
    size_of_test_image = 3
    image_data = np.zeros((size_of_test_image, size_of_test_image))
    image_data[1] = np.ones(size_of_test_image)

    x_pixel_coords = np.arange(size_of_test_image)
    legendre_polynomial_array = np.ones((1, size_of_test_image))
    legendre_polynomial_coefficients = np.array([1])
    expected_value = -1 * np.sum(image_data[1])

    image_filt = ndimage.spline_filter(image_data)
    found_value = trace_utils.crosscoef(legendre_polynomial_coefficients, image_filt,
                                        x_pixel_coords, legendre_polynomial_array)
    assert np.isclose(found_value, expected_value)

    testpoints = np.array([1])
    found_value = (-1) * trace_utils.fluxvalues(testpoints, [], image_filt, x_pixel_coords, legendre_polynomial_array)[0]
    assert np.isclose(found_value, expected_value)


def test_first_sorting_of_coefficients_from_blind_fit():



def test_splitting_coefficients_per_fiber():
    fake_coefficients_and_indices = np.array([[0, 1],
                                              [1, 2],
                                              [0, 3],
                                              [1, 4]])
    per_fiber_coeffs = trace_utils.split_already_sorted_coefficients_into_each_fiber(fake_coefficients_and_indices,
                                                                  num_lit_fibers=2)
    assert np.equal(per_fiber_coeffs[0], fake_coefficients_and_indices[:2]).all()
    assert np.equal(per_fiber_coeffs[1], fake_coefficients_and_indices[2:]).all()


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

    for x in range(num_trials):
        images = [FakeTraceImage()]
        images[0].readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(images[0], order_of_poly_fit=order_of_poly_fit)
        fill_image_with_traces(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        noisify_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))
        trim_image(images[0], trimmed_shape=tuple([min(images[0].data.shape)] * 2))

        maker = BlindTraceMaker(FakeContext())
        maker.order_of_poly_fit = order_of_poly_fit
        maker.do_stage(images)

        args, kwargs = mock_images.call_args
        master_trace = kwargs['data']
        logger.debug(master_trace.shape)

        difference = differences_between_found_and_generated_trace_vals(master_trace, images[0])
        logger.debug('error in unit-test trace fitting is less than %s of a pixel' %
                    np.median(np.abs(difference - np.median(difference))))
        logger.debug('worst error in unit-test trace fitting is %s pixels'%np.max(np.abs(difference)))
        logger.debug('systematic error (median difference) in unit-test trace fitting is less than %s of a pixel' %
                    np.abs(np.median(difference)))

        assert np.median(np.abs(difference - np.median(difference))) < 2/100
        assert np.abs(np.median(difference)) < 2/100
