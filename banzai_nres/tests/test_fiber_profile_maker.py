"""
Here we should generate a fake trace image with just two flat traces e.g an image which is like 100 by 50.
, give it the coefficients and measure
how far off the fit is from the true fiber profile guassian (including noise).
"""
from banzai_nres.tests.utils import FakeImage, noisify_image
from banzai_nres.tests.test_trace_maker import fill_image_with_traces

def test_fiber_profile_maker():
    assert True