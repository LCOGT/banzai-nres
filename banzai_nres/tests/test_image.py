from banzai_nres import images
from banzai.tests.utils import FakeContext


def test_image_class_loads():
    image = images.Image(pipeline_context=FakeContext())
    assert image.trace is None
