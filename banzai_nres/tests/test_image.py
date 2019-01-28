from banzai_nres.images import NRESImage
import banzai_nres.settings
from banzai.tests.utils import FakeContext
import mock


@mock.patch('banzai.images.Image._init_instrument_info')
def test_image_class_loads(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.trace is None

@mock.patch('banzai.images.Image._init_instrument_info')
def test_get_num_lit_fibers(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()),
                      header={'OBJECTS': 'tung&tung&none'})
    image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = False, False, False
    assert image.num_lit_fibers() == 0
