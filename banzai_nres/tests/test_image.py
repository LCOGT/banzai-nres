from banzai_nres.images import NRESImage
import banzai_nres.settings
from banzai.tests.utils import FakeContext


class AClass(object):
    def __init__(self):
        self.an_attribute = None


def test_image_class_loads():
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()))
    assert image.trace is None


def test_get_num_lit_fibers():
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()))
    image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = False, False, False
    assert image.num_lit_fibers() == 0
