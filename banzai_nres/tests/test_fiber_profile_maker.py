"""
Here we should generate a fake trace image with just two flat traces e.g an image which is like 100 by 50.
, give it the coefficients and measure
how far off the fit is from the true fiber profile guassian (including noise).
"""
import numpy as np
import mock

from banzai.tests.utils import FakeContext, GenerateFiberProfileImage
from banzai import logs

from banzai_nres.tests.utils import FakeImage, noisify_image, trim_image
from banzai_nres.tests.test_trace_maker import fill_image_with_traces, trim_coefficients_to_fit_image
from banzai_nres.fiber_profile import SampleFiberProfileAcrossImage


logger = logs.get_logger(__name__)


def _generate_image_with_two_flat_traces(readnoise=10, order_width=1.25):
    image = FakeImage(nx=1000, ny=102)
    trace_coefficients_no_indices = np.array([[image.data.shape[0]*1/3, 0, 0],
                                              [image.data.shape[0]*2/3, 0, 0]])
    image.trace.coefficients = trim_coefficients_to_fit_image(image, trace_coefficients_no_indices)
    fill_image_with_traces(image, order_width=order_width)
    image.readnoise = readnoise
    noisify_image(image, trimmed_shape=(100, 1000))
    trim_image(image, trimmed_shape=(100, 1000))
    return image


def test_fiber_profile_maker():
    real_full_width_half_max = 1.25
    images = [_generate_image_with_two_flat_traces(order_width=real_full_width_half_max)]
    sampling_stage = SampleFiberProfileAcrossImage(FakeContext())
    sampling_stage.do_stage(images)
    fiber_profile_maker_stage = GenerateFiberProfileImage(FakeContext())
    fiber_profile_maker_stage.do_stage(images)
    fwhm_estimate = images[0].median_full_width_half_max

    fwhm_abs_error = np.abs(real_full_width_half_max -  fwhm_estimate)
    fwhm_fractional_error = fwhm_abs_error/real_full_width_half_max
    logger.info('%s = |real_fwhm - fwhm_estimate|'%fwhm_abs_error)
    logger.info('%s = |real_fwhm - fwhm_estimate|/real_fwhm'%fwhm_fractional_error)
    assert (fwhm_fractional_error < 0.05)