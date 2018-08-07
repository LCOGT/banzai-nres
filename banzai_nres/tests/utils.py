from banzai_nres.images import Image
import datetime
from astropy.io import fits
import numpy as np


class FakeImage(Image):
    def __init__(self, nx=256, ny=259, image_multiplier=1.0,
                 ccdsum='2 2', epoch='20180807', n_amps=1):
        self.nx = nx
        self.ny = ny
        self.telescope_id = -1
        self.site = 'lsc'
        self.instrument = 'nres01'
        self.ccdsum = ccdsum
        self.epoch = epoch
        self.data = image_multiplier * np.ones((ny, nx), dtype=np.float32)
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


class FakeTraceImage(FakeImage):
    def __init__(self, *args, **kwargs):
        super(FakeImage, self).__init__(*args, **kwargs)
        self.caltype = 'trace'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'


def random_yet_realistic_trace_coefficients(nx, ny):
    characteristic_meta_coefficients_even = np.zeros((5, 6))
    characteristic_meta_coefficients_odd = np.zeros((5, 6))
    characteristic_meta_coefficients_even[0] = [1664, 1957, 394, 61, 5, 1.62]
    characteristic_meta_coefficients_even[1] = [-36.4, -50, -16.5, -3.31, -0.74, 0]
    characteristic_meta_coefficients_even[2] = [89.9, 1.69, 0.294, -0.112, 0, 1.62]
    characteristic_meta_coefficients_even[3] = [1664, 1957, 394, 61, 5, 1.62]
    characteristic_meta_coefficients_even[4] = [1664, 1957, 394, 61, 5, 1.62]
