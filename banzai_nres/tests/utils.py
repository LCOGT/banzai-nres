from banzai_nres.images import Image
import datetime
from astropy.io import fits
import numpy as np
from banzai_nres.utils.trace_utils import get_coefficients_from_meta, generate_legendre_array


class FakeImage(Image):
    def __init__(self, nx=256, ny=259, ccdsum='2 2', epoch='20180807', n_amps=1):
        self.nx = nx
        self.ny = ny
        self.telescope_id = -1
        self.site = 'lsc'
        self.instrument = 'nres01'
        self.ccdsum = ccdsum
        self.epoch = epoch
        self.data = np.ones((ny, nx))
        self.filename = 'test.fits'
        self.filter = 'U'
        self.dateobs = datetime(2018, 8, 7)
        self.header = {}
        self.caltype = ''
        self.bpm = np.zeros((ny, nx), dtype=np.uint8)
        self.request_number = '0000331403'
        self.readnoise = 11.0
        self.block_id = '254478983'
        self.molecule_id = '544562351'
        self.exptime = 30.0
        self.obstype = 'TEST'
        self.trace_fit_coefficients = None
        self.fiber_order = None


class FakeTraceImage(FakeImage):
    def __init__(self, *args, **kwargs):
        super(FakeImage, self).__init__(*args, **kwargs)
        self.caltype = 'trace'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'


def trim_coefficients_to_fit_image(image, trace_fit_coefficients_no_indices):
    return 0

def random_yet_realistic_trace_coefficients(image):
    """
    :param image: Banzai_nres Image object
    :return: realistic coefficients for traces which fit entirely in the frame.
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
        noise_scale = meta_coefficients_even[i, 0]/100
        noise = np.random.normal(loc=0,scale=noise_scale,size=meta_coefficients_even.shape[1])
        meta_coefficients_even[i] += noise
        meta_coefficients_odd[i] += noise

    meta_legendre_array, x, xnorm = generate_legendre_array(image, meta_coefficients_even.shape[1]-1)
    trace_coefficients_odd = get_coefficients_from_meta(meta_coefficients_odd, meta_legendre_array)
    trace_coefficients_even = get_coefficients_from_meta(meta_coefficients_even, meta_legendre_array)