from banzai_nres.images import NRESImage
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute
import banzai_nres.settings
from banzai.tests.utils import FakeContext


class AClass(object):
    def __init__(self):
        self.an_attribute = None


def test_image_class_loads():
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()))
    assert image.trace is None


def test_adding_class_as_attribute():
    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()))
    setattr(image, 'a_class', None)
    images = [image]
    add_class_as_attribute(images, 'a_class', AClass)
    assert images[0].a_class.an_attribute is None

    image = NRESImage(pipeline_context=FakeContext(settings=banzai_nres.settings.NRESSettings()))
    images = [image]
    add_class_as_attribute(images, 'a_class', AClass)
    assert images[0].a_class.an_attribute is None
