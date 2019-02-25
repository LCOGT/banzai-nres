from banzai_nres.utils import munge_utils
from banzai_nres.tests.utils import FakeImage


def test_get_telescope_filename():
    image = FakeImage()
    image.header['TELESCOP'] = 'nres01'
    assert munge_utils.get_telescope_filename(image) == 'nrs01'
