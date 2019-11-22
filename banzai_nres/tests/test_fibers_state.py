import mock

from banzai_nres.fibers import fiber_states_from_header, fibers_state_to_filename
from banzai_nres.images import NRESCalibrationFrame
from banzai.tests.utils import FakeContext


def test_creation_from_header():
    header = {'OBJECTS': 'targ&ThAr&none'}
    assert (1, 1, 0) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&targ'}
    assert (0, 1, 1) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&none'}
    assert (0, 1, 0) == fiber_states_from_header(header)


@mock.patch('banzai.images.Image._init_instrument_info')
def test_fiber_state_to_filename(mock_instrument):
    mock_instrument.return_value = None, None, None
    image = NRESImage(runtime_context=FakeContext(),
                      header={'OBJECTS': 'tung&tung&none'})
    assert fibers_state_to_filename(image) == '110'
