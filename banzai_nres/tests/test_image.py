from banzai_nres.frames import NRESObservationFrame, EchelleSpectralCCDData, Spectrum1D
from banzai.data import HeaderOnly
from banzai import context

import numpy as np


def test_get_num_lit_fibers():
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    assert image.num_lit_fibers() == 2
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'none&tung&none'})], 'foo.fits')
    assert image.num_lit_fibers() == 1


def test_to_fits():
    data = np.ones((2, 2))
    spec = Spectrum1D({'fiber': 0, 'order': 1, 'flux': np.arange(10)})
    image = EchelleSpectralCCDData(data=data, uncertainty=2 * data,
                                   traces=3 * data, weights=4 * data,
                                   background=5 * data, spectrum=spec,
                                   meta={'OBJECTS': 'tung&tung&none', 'EXTNAME': ''})
    hdulist = image.to_fits(context.Context({'EXTENSION_NAMES_TO_CONDENSE': []}))
    attribute = {'TRACES': 'traces', 'BACKGROUND': 'background', '1DSPEC': 'spectrum',
                 'WEIGHTS': 'weights'}
    hdu_ext_names = [hdulist[i].header['EXTNAME'] for i in range(len(hdulist))]
    for name in attribute.keys():
        i = np.where([n == name for n in hdu_ext_names])[0][0]
        if name != '1DSPEC':
            assert np.allclose(hdulist[i].data, getattr(image, attribute[name]))
        else:
            assert np.allclose(hdulist[i].data[0]['flux'], image.spectrum[0, 1]['flux'])
