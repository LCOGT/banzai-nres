from banzai_nres.utils import runtime_utils
from banzai_nres.frames import NRESObservationFrame
from banzai.data import HeaderOnly


def test_get_telescope_filename():
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'tung&tung&none', 'TELESCOP': 'nres01'}, name='test')],
                                 'foo.fits')
    assert runtime_utils.get_telescope_filename(image) == 'nrs01'
