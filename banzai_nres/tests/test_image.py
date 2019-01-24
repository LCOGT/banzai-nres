from banzai_nres.images import NRESImage
import banzai_nres.settings
from banzai.tests.utils import FakeContext


class AClass(object):
    def __init__(self):
        self.an_attribute = None


def test_image_class_loads():
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()))
    assert image.trace is None
