from banzai_nres.utils import runtime_utils
from banzai_nres.images import NRESObservationFrame
from banzai.images import HeaderOnly


def test_get_telescope_filename():
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'tung&tung&none', 'TELESCOP': 'nres01'})], 'foo.fits')
    assert runtime_utils.get_telescope_filename(image) == 'nrs01'
