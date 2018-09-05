from datetime import datetime
import numpy as np

from banzai_nres.traces import Trace
from banzai_nres.images import Image
from banzai_nres.coordinate_transform import Coordinates
from banzai_nres.fiber_profile import FiberProfile
from banzai_nres.extraction import Spectra

from banzai_nres.utils.trace_utils import get_trace_centroids_from_coefficients

"""
General test utils:
"""


class FakeImage(Image):
    def __init__(self, nx=102, ny=100, overscan_size=2, ccdsum='2 2', epoch='20180807'):
        self.nx = nx
        self.ny = ny
        self.telescope_id = -1
        self.site = 'lsc'
        self.instrument = 'nres01'
        self.ccdsum = ccdsum
        self.epoch = epoch
        self.overscan_size = overscan_size
        self.data = np.ones((ny, nx))
        self.filename = 'test.fits'
        self.filter = 'U'
        self.dateobs = datetime(2018, 8, 7)
        self.header = {}
        self.caltype = ''
        self.bpm = np.zeros((ny, nx-overscan_size), dtype=np.uint8)
        self.request_number = '0000331403'
        self.readnoise = 11
        self.block_id = '254478983'
        self.molecule_id = '544562351'
        self.exptime = 30.0
        self.obstype = 'TEST'

        self.trace = Trace()
        self.coordinates = Coordinates()
        self.fiber_profile = FiberProfile()
        self.spectra = Spectra()
        self.ivar = None


def gaussian(x, A, b, sigma):
    return A * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def noisify_image(image, trimmed_shape):
    """
    :param image: Banzai_nres FakeImage object.
    This adds poisson, readnoise to an image with traces already on it, in that order.
    """
    # poisson noise
    poissonnoise_mask = np.random.poisson(image.data[:trimmed_shape[0], :trimmed_shape[1]])
    image.data[:trimmed_shape[0], :trimmed_shape[1]] += poissonnoise_mask
    # read noise
    image.data += np.random.normal(0, image.readnoise, image.data.shape)


def trim_image(image, trimmed_shape):
    """
    :param image:
    Squares up the image, thus fake images which are not square are a bad idea.
    this trim_image may be unneccessary.
    """
    image.data = image.data[:trimmed_shape[0], :trimmed_shape[1]]
    image.ny, image.nx = trimmed_shape


def append_good_region_info(image):
    image.trace.has_sufficient_signal_to_noise = np.ones(len(image.trace.coefficients)).astype(bool)
    image.trace.high_signal_to_noise_region_bounds = np.ones((len(image.trace.coefficients), 2))
    image.trace.high_signal_to_noise_region_bounds[:, 0] = 0
    image.trace.high_signal_to_noise_region_bounds[:, 1] = image.data.shape[1]-1


def fill_with_simple_inverse_variances(image):
    read_var = image.readnoise ** 2
    shot_var = image.data.astype(np.float32)
    vars = shot_var + read_var
    vars[vars < read_var] = read_var
    image.ivar = np.reciprocal(vars)


"""
Trace related test utils:
"""


def generate_image_with_two_flat_traces(readnoise=10, order_width=1.25, normalized_traces=False, add_noise=True):
    overscan_size = 2
    nx = 1000 + overscan_size
    ny = 50
    trimmed_shape = (ny, nx - overscan_size)
    image = FakeImage(nx=nx, ny=(ny+2))
    trace_coefficients_no_indices = np.array([[image.data.shape[0]*1/3, 0, 0],
                                              [image.data.shape[0]*2/3, 0, 0]])

    image.trace.coefficients = trim_coefficients_to_fit_image(image, trace_coefficients_no_indices)
    if normalized_traces:
        gaussian_norm_factor = 1/np.sqrt(2 * np.pi * order_width ** 2)
        fill_image_with_traces(image, trimmed_shape=trimmed_shape, order_width=order_width,
                               odd_fiber_intensity=gaussian_norm_factor, even_fiber_intensity=gaussian_norm_factor)
    if not normalized_traces:
        # adopt standard intensities which mimic a 120 sec frame
        fill_image_with_traces(image, trimmed_shape=trimmed_shape, order_width=order_width)
    image.readnoise = readnoise
    if add_noise:
        noisify_image(image, trimmed_shape=trimmed_shape)
    trim_image(image, trimmed_shape=trimmed_shape)
    return image


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