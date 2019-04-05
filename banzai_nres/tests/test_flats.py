from banzai.tests.utils import FakeContext
from banzai_nres.flats import FlatStacker


def test_flatstacker():
    assert FlatStacker(FakeContext()).calibration_type == 'LAMPFLAT'
