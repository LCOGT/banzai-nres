from banzai_nres.images import NRESImage
from banzai.tests.utils import FakeContext

import mock


@mock.patch('banzai.images.Image._init_instrument_info')
def test_image_class_loads(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.trace is None
    assert image.rectified_2d_spectrum is None


@mock.patch('banzai.images.Image._init_instrument_info')
def test_get_num_lit_fibers(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert image.num_lit_fibers() == 2
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'none&tung&none'})
    assert image.num_lit_fibers() == 1
