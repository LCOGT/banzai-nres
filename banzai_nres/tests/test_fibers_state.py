from banzai_nres.fibers import FibersState


def test_fiber_states_equal():
    f1 = FibersState(True, False, False)
    f2 = FibersState(True, False, False)
    assert f1 == f2
    assert f1 >= f2
    assert f1 <= f2
    f1 = FibersState(True, True, False)
    f2 = FibersState(True, True, False)
    assert f1 == f2
    assert f1 >= f2
    assert f1 <= f2
    f1 = FibersState(True, False, True)
    f2 = FibersState(True, False, True)
    assert f1 == f2
    assert f1 >= f2
    assert f1 <= f2
    f1 = FibersState(False, False, True)
    f2 = FibersState(False, False, True)
    assert f1 == f2
    assert f1 >= f2
    assert f1 <= f2


def test_fiber_states_not_equal_and_inequalities():
    f1 = FibersState(True, True, False)
    f2 = FibersState(True, False, False)
    assert f1 != f2
    assert f1 > f2
    f1 = FibersState(False, True, False)
    f2 = FibersState(True, True, False)
    assert f1 != f2
    assert f1 < f2
    f1 = FibersState(True, False, True)
    f2 = FibersState(True, False, False)
    assert f1 != f2
    assert f1 > f2
    f1 = FibersState(False, False, True)
    f2 = FibersState(False, False, False)
    assert f1 != f2
    assert f1 > f2


def test_on_off_off_str():
    fibers_state = FibersState(True, False, False)
    assert str(fibers_state) == '100'


def test_off_on_off_str():
    fibers_state = FibersState(False, True, False)
    assert str(fibers_state) == '010'


def test_off_off_on_str():
    fibers_state = FibersState(False, False, True)
    assert str(fibers_state) == '001'


def test_creation_from_header():
    header = {'OBJECTS': 'targ&ThAr&none'}
    assert FibersState(True, True, False) == FibersState.from_header(header)

    header = {'OBJECTS': 'none&ThAr&targ'}
    assert FibersState(False, True, True) == FibersState.from_header(header)

    header = {'OBJECTS': 'none&ThAr&none'}
    assert FibersState(False, True, False) == FibersState.from_header(header)
