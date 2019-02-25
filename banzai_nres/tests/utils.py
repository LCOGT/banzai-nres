from datetime import datetime
import numpy as np
from astropy.io import fits

from banzai_nres.utils.trace_utils import Trace, SingleTraceFitter


class FakeImage(object):
    def __init__(self, nx=102, ny=100, overscan_size=2, ccdsum='2 2', epoch='20180807'):
        self.nx = nx
        self.ny = ny
        self.telescope_id = -1
        self.site = 'lsc'
        self.instrument = 'nres01'
        self.camera = None
        self.ccdsum = ccdsum
        self.epoch = epoch
        self.overscan_size = overscan_size
        self.data = np.ones((ny, nx))
        self.filename = 'test.fits'
        self.filter = 'U'
        self.dateobs = datetime(2018, 8, 7)
        self.header = fits.Header({'RDNOISE': 11})
        self.caltype = ''
        self.bpm = np.zeros((ny, nx-overscan_size), dtype=np.uint8)
        self.request_number = '0000331403'
        self.block_id = '254478983'
        self.molecule_id = '544562351'
        self.exptime = 30.0
        self.obstype = 'TEST'
        self.data_tables = None
        self.trace = Trace(num_centers_per_trace=nx)
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = False, True, True


def gaussian(x, A, b, sigma):
    return A * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def noisify_image(image):
    """
    :param image: Banzai_nres FakeImage object.
    This adds poisson and read noise to an image with traces already on it, in that order.
    """
    image.data = np.random.poisson(image.data) + np.random.normal(0, image.header['RDNOISE'], image.data.shape)


def array_with_peaks(x, centroids, amplitudes, stds):
    vectorized_gaussian = np.vectorize(gaussian)
    y = np.zeros_like(x)
    for centroid, amplitude, std in zip(centroids, amplitudes, stds):
        y += vectorized_gaussian(x, amplitude, centroid, std)
    return y


def fill_image_with_traces(image, poly_fit_order=4, order_width=1.5, fiber_intensity=1E4):
    if image.data.shape[1] < 200:
        raise ValueError
    trace_fitter = SingleTraceFitter(image_data=image.data,
                                     second_order_coefficient_guess=0,
                                     poly_fit_order=poly_fit_order)
    num_fake_traces = int((image.data.shape[1] - 60)/20)
    coefficients = np.ones((num_fake_traces, poly_fit_order+1))
    coefficients[:, 0] = np.linspace(30, image.data.shape[0] - 30, num=num_fake_traces)
    coefficients[:, 2] = np.linspace(30, 40, num=num_fake_traces)
    coefficients[:, 3] = np.linspace(1, 3, num=num_fake_traces)
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
    image.data += trace_overlay*fiber_intensity
    second_order_coefficient_guess = np.mean(coefficients[:, 2])
    return image, trace_centers, second_order_coefficient_guess
