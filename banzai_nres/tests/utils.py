from datetime import datetime
import numpy as np

from banzai_nres.utils.trace_utils import Trace


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
        self.header = {}
        self.caltype = ''
        self.bpm = np.zeros((ny, nx-overscan_size), dtype=np.uint8)
        self.request_number = '0000331403'
        self.readnoise = 11
        self.block_id = '254478983'
        self.molecule_id = '544562351'
        self.exptime = 30.0
        self.obstype = 'TEST'
        self.data_tables = None
        self.trace = Trace()
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = False, True, True


def gaussian(x, A, b, sigma):
    return A * np.exp(-(x - b) ** 2 / (2 * sigma ** 2))


def noisify_image(image):
    """
    :param image: Banzai_nres FakeImage object.
    This adds poisson and read noise to an image with traces already on it, in that order.
    """
    image.data = np.random.poisson(image.data) + np.random.normal(0, image.readnoise, image.data.shape)


def array_with_two_peaks():
    """
    :return: generates a fake signal with two peaks with height 1.
    """
    centroids = (15, 30)
    x = np.linspace(0, 50, num=100)
    vectorized_gaussian = np.vectorize(gaussian)
    y = vectorized_gaussian(x, 1, centroids[0], 1.5) + vectorized_gaussian(x, 1, centroids[1], 1.5)
    return y, centroids, x


def fill_image_with_traces(image):
    image.data[int(image.data.shape[0]/2), :] = 1E5
