from banzai_nres.fibers import fiber_states_from_header, fibers_state_to_filename
from banzai_nres.frames import NRESObservationFrame
from banzai.data import HeaderOnly


def test_creation_from_header():
    header = {'OBJECTS': 'targ&ThAr&none'}
    assert (1, 1, 0) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&targ'}
    assert (0, 1, 1) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&none'}
    assert (0, 1, 0) == fiber_states_from_header(header)


def test_fiber_state_to_filename():
    image = NRESObservationFrame([HeaderOnly(meta={'OBJECTS': 'tung&tung&none'}, name='fibertest')], 'foo.fits')
    assert fibers_state_to_filename(image) == '110'
