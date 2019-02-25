import mock

from banzai_nres.fibers import fiber_states_from_header, fibers_state_to_filename
from banzai_nres.images import NRESImage
import banzai_nres.settings
from banzai.tests.utils import FakeContext


def test_creation_from_header():
    header = {'OBJECTS': 'targ&ThAr&none'}
    assert (True, True, False) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&targ'}
    assert (False, True, True) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&none'}
    assert (False, True, False) == fiber_states_from_header(header)


@mock.patch('banzai.images.Image._init_instrument_info')
def test_fiber_state_to_filename(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()),
                      header={'OBJECTS': 'tung&tung&none'})
    assert fibers_state_to_filename(image) == '110'
