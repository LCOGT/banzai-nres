import pytest
import mock
from banzai_nres.traces import BlindTraceMaker
from banzai.tests.utils import FakeContext
import numpy as np
from banzai import logs
from banzai_nres.utils.trace_utils import get_coefficients_from_meta, generate_legendre_array, get_trace_centroids_from_coefficients
from banzai_nres.tests.utils import FakeImage, noisify_image, trim_image
from astropy.io import fits


logger = logs.get_logger(__name__)


class FakeTraceImage(FakeImage):
    def __init__(self, *args, **kwargs):
        super(FakeTraceImage, self).__init__(*args, **kwargs)
        self.caltype = 'trace'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'
        # Note: Image must be at least 400x400 for enough traces to populate to test 'global-meta' procedure
        self.nx = 500
        self.ny = 502
        self.bpm = np.zeros((self.ny, self.nx), dtype=np.uint8)
        self.data = np.zeros((self.ny, self.nx))


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


def gaussian(x, A, b, sigma):
    return A * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def fill_image_with_traces(image, order_width=1.25):
    """
    :param image: Banzai_nres FakeImage object which is square after trimming.
    :param order_width : the standard deviation of an unnormalized guassian e.g. sigma
           from np.exp(-(x - b) ** 2 / (2 * sigma ** 2))
    fills the image.data with guassian profiled trace coefficients.
    """
    image.data = np.zeros_like(image.data)
    trimmed_shape = tuple([min(image.data.shape)]*2)
    even_fiber = np.zeros(trimmed_shape)
    odd_fiber = np.zeros(trimmed_shape)

    trace_values_versus_xpixel, num_traces, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
    vgauss = np.vectorize(gaussian)  # prepare guassian for evaluation along a slice centered at each trace point.
    odd_fiber_intensity, even_fiber_intensity = 1E4, 5E3
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


@mock.patch('banzai_nres.traces.Image')
def test_blind_trace_maker(mock_images):
    """
    This tests blind trace making (which involves blind trace making and then refining via meta-fit).
    Note: The fake images made here must have enough orders N such that N > order_of_meta_fit + 1.
    Currently it creates ~20 traces so 10 orders. A bigger image would have more traces, but expanding the
    image size may cause the created traces to be unrealistic. My suggestion is to never change the fake_image size
    and to keep order_of_meta_fit less than 8. (6 works well so why would you want it larger?)
    """
    num_trials = 2
    readnoise = 11.0
    order_of_poly_fit = 4

    for x in range(num_trials):
        images = [FakeTraceImage()]
        images[0].readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(images[0], order_of_poly_fit=order_of_poly_fit)
        fill_image_with_traces(images[0])
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

        assert np.median(np.abs(difference - np.median(difference))) < 1/100
        assert np.abs(np.median(difference)) < 1/100