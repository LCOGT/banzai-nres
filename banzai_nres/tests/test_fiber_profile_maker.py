"""
Here we should generate a fake trace image with just two flat traces e.g an image which is like 100 by 50.
, give it the coefficients and measure
how far off the fit is from the true fiber profile guassian (including noise).
"""
import numpy as np

from banzai.tests.utils import FakeContext
from banzai import logs

from banzai_nres.tests.utils import FakeImage, noisify_image, trim_image, \
    append_good_region_info, fill_with_simple_inverse_variances
from banzai_nres.tests.test_trace_maker import fill_image_with_traces, trim_coefficients_to_fit_image
from banzai_nres.fiber_profile import SampleFiberProfileAcrossImage, GenerateFiberProfileImage
from banzai_nres.coordinate_transform import MakeTraceCentricCoordinates


logger = logs.get_logger(__name__)


def generate_image_with_two_traces(readnoise=10, order_width=1.25):
    overscan_size = 2
    nx = 1000 + overscan_size
    ny = 50
    trimmed_shape = (ny, nx - overscan_size)
    image = FakeImage(nx=nx, ny=(ny+2))
    trace_coefficients_no_indices = np.array([[image.data.shape[0]*1/3, 0, 0],
                                              [image.data.shape[0]*2/3, 0, 0]])

    image.trace.coefficients = trim_coefficients_to_fit_image(image, trace_coefficients_no_indices)
    fill_image_with_traces(image, trimmed_shape=trimmed_shape, order_width=order_width)
    image.readnoise = readnoise
    noisify_image(image, trimmed_shape=trimmed_shape)
    trim_image(image, trimmed_shape=trimmed_shape)
    return image


def test_fiber_profile_maker():
    """
    tests the accuracy of the profile fitter by fitting a guassian profile about known trace centers.
    Because we are fitting a gaussian (Which is normalized) the guassian is completely determined by the sigma parameter
    so if the found_sigma and the real_sigma (called here full width half maxes fwhm) are equal, then we have matched
    the true profile.
    """
    real_full_width_half_max = 1.25
    image = generate_image_with_two_traces(order_width=real_full_width_half_max)
    append_good_region_info(image)

    fill_with_simple_inverse_variances(image)
    images = [image]
    # appending coordinate info
    coordinate_stage = MakeTraceCentricCoordinates(FakeContext())
    images = coordinate_stage.do_stage(images)
    sampling_stage = SampleFiberProfileAcrossImage(FakeContext())
    # modifying fit parameters to work for such a small image.
    sampling_stage.wing_intervals = 1
    sampling_stage.middle_intervals = 1
    sampling_stage.size_of_basis = 2

    images = sampling_stage.do_stage(images)
    fiber_profile_maker_stage = GenerateFiberProfileImage(FakeContext())
    images = fiber_profile_maker_stage.do_stage(images)
    fwhm_estimate = images[0].median_full_width_half_max

    fwhm_abs_error = np.abs(real_full_width_half_max - fwhm_estimate)
    fwhm_fractional_error = fwhm_abs_error/real_full_width_half_max
    logger.info('%s = |real_fwhm - fwhm_estimate|' % fwhm_abs_error)
    logger.info('%s = |real_fwhm - fwhm_estimate|/real_fwhm' % fwhm_fractional_error)
    assert (fwhm_fractional_error < 0.01)
