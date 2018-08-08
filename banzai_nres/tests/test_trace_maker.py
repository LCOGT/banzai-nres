import pytest
import mock
from banzai_nres.traces import BlindTraceMaker
from banzai.tests.utils import FakeContext
import numpy as np
from banzai import logs
from banzai.utils import stats
from banzai_nres.utils.trace_utils import get_coefficients_from_meta, generate_legendre_array, get_trace_centroids_from_coefficients
from banzai_nres.tests.utils import FakeImage
from astropy.io import fits


logger = logs.get_logger(__name__)


class FakeTraceImage(FakeImage):
    def __init__(self, *args, **kwargs):
        super(FakeTraceImage, self).__init__(*args, **kwargs)
        self.caltype = 'trace'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'


def trim_coefficients_to_fit_image(image, trace_fit_coefficients_no_indices):
    min_y, max_y = 0, image.data.shape[0]
    order_indices = np.array([i for i in range(0, trace_fit_coefficients_no_indices.shape[0])])
    trace_fit_coefficients = np.insert(trace_fit_coefficients_no_indices, obj=0, values=order_indices, axis=1)
    trace_values_versus_xpixel, num_traces, x = get_trace_centroids_from_coefficients(trace_fit_coefficients, image)
    good_indices = []
    for i in range(trace_values_versus_xpixel.shape[0]):
        if (trace_values_versus_xpixel[i, :] < max_y).all() and (trace_values_versus_xpixel[i, :] > min_y).all():
            good_indices += [i]
    trimmed_trace_fit_coefficients_and_indices = trace_fit_coefficients[good_indices]
    assert np.array(good_indices) - np.min(np.array(good_indices)) == np.array(list(range(len(good_indices))))
    # replacing order indices with proper indicators
    trimmed_trace_fit_coefficients_and_indices[:, 0] = np.array(list(range(len(good_indices))))
    return trimmed_trace_fit_coefficients_and_indices


def munge_coefficients(even_coefficients, odd_coefficients):
    if even_coefficients.shape[0] != odd_coefficients.shape[0]:
        min_shape = min(odd_coefficients.shape[0], even_coefficients.shape[0]) - 1
    return even_coefficients[:min_shape], odd_coefficients[:min_shape]


def make_random_yet_realistic_trace_coefficients(image):
    """
    :param image: Banzai_nres Image object
    Adds realistic coefficients for traces which fit entirely in the frame, saving into image.trace_fit_coefficients
    and an arbitrary fiber_order onto image.fiber_order
    """
    meta_coefficients_even = np.zeros((5, 6))
    meta_coefficients_even[0] = [1664.8, 1957, 394, 61, 5, 1.62]
    meta_coefficients_even[1] = [-36.4, -50, -16.5, -3.31, -0.74, 0]
    meta_coefficients_even[2] = [89.9, 1.69, 0.294, -0.112, 0, 1.62]
    meta_coefficients_even[3] = [-0.234, -0.264, -0.114, 0, 0, 0]
    meta_coefficients_even[4] = [0.581, -0.374, 0.188, 0, 0.0207, 0]
    meta_coefficients_odd = np.copy(meta_coefficients_even)
    meta_coefficients_odd[0] = [1668.1, 1940, 386, 59, 4.9, 1.37]
    for i in range(1, meta_coefficients_even.shape[0]):
        noise_scale = np.abs(np.median(meta_coefficients_even[i, 0]))/100
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

    image.fiber_order = (0, 1)
    image.trace_fit_coefficients = np.vstack((trace_coefficients_even_and_indices, trace_coefficients_odd_and_indices))


def gaussian(x, A, b, sigma):
    return A * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def fill_image_with_traces(image):
    """
    :param image: Banzai_nres FakeImage object which is square after trimming.
    fills the image.data with guassian profiled trace coefficients.
    """
    image.data = np.zeros_like(image.data)
    trimmed_shape = tuple([min(image.data.shape)]*2)
    even_fiber = np.zeros(trimmed_shape)
    odd_fiber = np.zeros(trimmed_shape)

    trace_values_versus_xpixel, num_traces, x = get_trace_centroids_from_coefficients(image.trace_fit_coefficients, image)
    vgauss = np.vectorize(gaussian)  # prepare guassian for evaluation along a slice centered at each trace point.
    order_width, odd_fiber_intensity, even_fiber_intensity = 1.8, 1E4, 5E3
    # these are realistic intensity values according to a 120 sec LSC exposure.
    for x_pixel in range(even_fiber.shape[1]):
        for i in range(num_traces):
            centroid = trace_values_versus_xpixel[i, x_pixel]
            low, high = max(0, int(centroid - 4 * order_width)), min(even_fiber.shape[1]-1, int(centroid + 4 * order_width)) + 1
            evalwindow = np.arange(low, high, 1)
            if len(evalwindow) > 0:
                if i % 2 == 0:
                    even_fiber[low:high, x_pixel] += vgauss(evalwindow, 1, centroid, order_width)
                else:
                    odd_fiber[low:high, x_pixel] += vgauss(evalwindow, 1, centroid, order_width)

    # add traces
    image.data[:trimmed_shape[0], :trimmed_shape[1]] += (odd_fiber_intensity* odd_fiber + even_fiber_intensity * even_fiber)


def add_overscan_and_noisify_image(image):
    """
    :param image: Banzai_nres FakeImage object.
    This adds poisson, overscan and readnoise to an image with traces already on it, in that order.
    """
    trimmed_shape = tuple([min(image.data.shape)] * 2)
    # poisson noise
    poissonnoise_mask = np.random.poisson(image.data[:trimmed_shape[0], :trimmed_shape[1]])
    image.data[:trimmed_shape[0], :trimmed_shape[1]] += poissonnoise_mask
    # overscan
    image.data += 800
    # read noise
    image.data += np.random.normal(0, image.readnoise, image.data.shape)


def differences_between_found_and_generated_trace_vals(found_coefficients, image):
    trace_values_1, num_traces_1, x = get_trace_centroids_from_coefficients(image.trace_fit_coefficients, image)
    trace_values_2, num_traces_2, x = get_trace_centroids_from_coefficients(found_coefficients, image)
    assert num_traces_1 == num_traces_2
    return trace_values_2 - trace_values_1


@mock.patch('banzai_nres.traces.Image')
def test_blind_trace_maker(mock_images):
    num_trials = 2
    readnoise = 11.0

    for x in range(num_trials):
        images = [FakeTraceImage()]
        images[0].readnoise = readnoise

        make_random_yet_realistic_trace_coefficients(images[0])
        fill_image_with_traces(images[0])
        add_overscan_and_noisify_image(images[0])

        maker = BlindTraceMaker(FakeContext())
        maker.do_stage(images)

        args, kwargs = mock_images.call_args
        master_trace = kwargs['data']
        logger.debug(master_trace.shape)

        difference = differences_between_found_and_generated_trace_vals(master_trace, images[0])
        logger.info('error in trace fitting is less than %s of a pixel' % stats.absolute_deviation(np.abs(difference)))
        logger.info('systematic error (median difference) in trace fitting is less than %s of a pixel' %
                    np.abs(np.median(difference)))

        assert stats.absolute_deviation(np.abs(difference)) < 1/10
        assert np.abs(np.median(difference)) < 1/100
        assert False
