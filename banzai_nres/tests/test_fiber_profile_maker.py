"""
Here we should generate a fake trace image with just two flat traces e.g an image which is like 100 by 50.
, give it the coefficients and measure
how far off the fit is from the true fiber profile guassian (including noise).
"""
import numpy as np

from banzai.tests.utils import FakeContext
from banzai import logs

from banzai_nres.tests.utils import append_good_region_info, fill_with_simple_inverse_variances,\
    generate_image_with_two_flat_traces
from banzai_nres.fiber_profile import SampleFiberProfileAcrossImage, GenerateFiberProfileImage
from banzai_nres.coordinate_transform import MakeTraceCentricCoordinates


logger = logs.get_logger(__name__)


def test_fiber_profile_maker():
    """
    tests the accuracy of the profile fitter by fitting a guassian profile about known trace centers.
    Because we are fitting a gaussian (Which is normalized) the guassian is completely determined by the sigma parameter
    so if the found_sigma and the real_sigma (called here full width half maxes fwhm) are equal, then we have matched
    the true profile.
    """
    real_full_width_half_max = 1.25
    image = generate_image_with_two_flat_traces(order_width=real_full_width_half_max)
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
