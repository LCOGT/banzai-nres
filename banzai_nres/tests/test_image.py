from banzai_nres.frames import NRESObservationFrame, Spectrum1D
from banzai.data import HeaderOnly, CCDData, ArrayData, DataTable
from banzai import context
import numpy as np


def test_get_num_lit_fibers():
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    assert image.num_lit_fibers() == 2
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'none&tung&none'})], 'foo.fits')
    assert image.num_lit_fibers() == 1


def test_to_fits():
    data = np.ones((2, 2))
    spec = Spectrum1D({'fiber': [0], 'order': [1], 'flux': [np.arange(10)], 'wavelength': [np.arange(10)]})
    image = NRESObservationFrame([CCDData(data=data, uncertainty=2 * data, name='SCI',
                                          meta={'OBJECTS': 'tung&tung&none', 'EXTNAME': ''})], 'foo.fits')
    image['TRACES'] = ArrayData(3 * data, name='TRACES')
    image['WEIGHTS'] = ArrayData(4 * data, name='WEIGHTS')
    image['BACKGROUND'] = ArrayData(5 * data, name='BACKGROUND')
    image['SPECTRUM'] = DataTable(spec.table, name='SPECTRUM')
    hdulist = image.to_fits(context.Context({'EXTENSION_NAMES_TO_CONDENSE': ['SCI'],
                                             'REDUCED_DATA_EXTENSION_TYPES': {},
                                             'fpack': True,
                                             'LOSSLESS_EXTENSIONS': []}))
    hdu_ext_names = [hdulist[i].header.get('EXTNAME') for i in range(len(hdulist))
                     if hdulist[i].header.get('EXTNAME') is not None]
    for name in ['SCI', 'BPM', 'ERR', 'TRACES', 'WEIGHTS', 'BACKGROUND', 'SPECTRUM']:
        assert name in hdu_ext_names
