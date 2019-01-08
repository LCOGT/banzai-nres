from banzai_nres.images import NRESImage
from banzai.tests.utils import FakeContext
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute


class AClass(object):
    def __init__(self):
        self.an_attribute = None


def test_image_class_loads():
    image = NRESImage(pipeline_context=FakeContext())
    assert image.trace is None


def test_adding_class_as_attribute():
    image = NRESImage(pipeline_context=FakeContext())
    setattr(image, 'a_class', None)
    images = [image]
    add_class_as_attribute(images, 'a_class', AClass)
    assert images[0].a_class.an_attribute is None

    image = NRESImage(pipeline_context=FakeContext())
    images = [image]
    add_class_as_attribute(images, 'a_class', AClass)
    assert images[0].a_class.an_attribute is None
