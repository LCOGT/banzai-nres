from banzai_nres.images import NRESObservationFrame
from banzai.images import HeaderOnly


def test_get_num_lit_fibers():
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    assert image.num_lit_fibers() == 2
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'none&tung&none'})], 'foo.fits')
    assert image.num_lit_fibers() == 1
