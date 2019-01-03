from banzai_nres.fibers import fiber_states_from_header


def test_creation_from_header():
    header = {'OBJECTS': 'targ&ThAr&none'}
    assert (True, True, False) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&targ'}
    assert (False, True, True) == fiber_states_from_header(header)

    header = {'OBJECTS': 'none&ThAr&none'}
    assert (False, True, False) == fiber_states_from_header(header)
