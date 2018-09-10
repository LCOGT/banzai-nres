from datetime import datetime
import numpy as np

from banzai_nres.traces import Trace
from banzai_nres.images import Image


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
