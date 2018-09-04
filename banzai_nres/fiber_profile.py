"""
Scripts for fitting the fiber intensity perpendicular to the dispersion.
Author
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import numpy as np
from scipy import ndimage


from banzai_nres.images import Image
from banzai_nres.utils import fiber_profile_utils
from banzai_nres.utils.trace_utils import get_trace_centroids_from_coefficients
from banzai_nres.utils.fiber_profile_utils import fit_vertical_fiber_intensity_functions_over_given_horizontal_ranges
from banzai_nres.utils.array_utils import fill_with_nearest_left_value_if_flagged_as_false
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute

from banzai import logs
from banzai.stages import CalibrationMaker, Stage, MasterCalibrationDoesNotExist


logger = logs.get_logger(__name__)


class FiberProfile(object):
    def __init__(self):
        self.fit_coefficients = None
        self.horizontal_ranges = None
        self.normalized_fiber_profile_image = None
        self.median_full_width_half_max = 1.5


class FiberStage(Stage):
    """
    the model with which one fits the profiles can be changed via the self.fiber_profile_model attribute.
    the default option is Shapelets. See fiber_profile_utils.py for the Shapelets Class.
    """

    def __init__(self, pipeline_context):
        super(FiberStage, self).__init__(pipeline_context)
        self.fiber_profile_model = fiber_profile_utils.Shapelets

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'fiber_profile'

    def do_stage(self, images):
        pass


class SampleFiberProfileAcrossImage(FiberStage):
    """
    This stage samples the fiber profile of across the detector, for each trace.
    """
    def __init__(self, pipeline_context):
        super(SampleFiberProfileAcrossImage, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'fiber_profile'

    def do_stage(self, images):
        add_class_as_attribute(images, 'fiber_profile', FiberProfile)
        for image in images:
            filtered_image_data = ndimage.spline_filter(image.data)
            trace_centroid_values_per_order, num_traces, x = get_trace_centroids_from_coefficients(
                image.trace.coefficients,
                image)
            fiber_profile_coeffs_per_trace = []
            horizontal_ranges_per_profile_fit_per_trace = []
            for trace_number in range(num_traces):
                if image.trace.has_sufficient_signal_to_noise[trace_number]:
                    horizontal_ranges = horizontal_windows_to_fit_fiber_profile(image, trace_number)
                    coeffs_at_each_window = fit_vertical_fiber_intensity_functions_over_given_horizontal_ranges(
                        filtered_image_data, image, trace_centroid_values_per_order, trace_number, horizontal_ranges,
                        Model=self.fiber_profile_model)

                    normalized_coeffs = fiber_profile_utils.normalize_fiber_fits(coeffs_at_each_window)
                    fiber_profile_coeffs_per_trace.append(normalized_coeffs)
                    horizontal_ranges_per_profile_fit_per_trace.append(horizontal_ranges)
                if not image.trace.has_sufficient_signal_to_noise[trace_number]:
                    fiber_profile_coeffs_per_trace.append([None])
                    horizontal_ranges_per_profile_fit_per_trace.append([None])
            # Below: filling fiber profiles for dim traces with the adjacent profile from the same fiber.
            # NOTE: if the zeroth order of the first fiber is dim, then the zeroth order of the second fiber
            # will incorrectly inherit the profile from the 66th order of the first fiber.

            fill_with_nearest_left_value_if_flagged_as_false(image.trace.has_sufficient_signal_to_noise,
                                                             horizontal_ranges_per_profile_fit_per_trace)
            fill_with_nearest_left_value_if_flagged_as_false(image.trace.has_sufficient_signal_to_noise,
                                                             fiber_profile_coeffs_per_trace)

            image.fiber_profile.fit_coefficients = np.array(fiber_profile_coeffs_per_trace)
            image.fiber_profile.horizontal_ranges = np.array(horizontal_ranges_per_profile_fit_per_trace)


class GenerateFiberProfileImage(FiberStage):
    """
    This stage interpolates sampled fiber profiles by interpolating the parameters which describe their best fit.
    This then uses the interpolated parameters to evaluate the profile across the detector (e.g. at 4096 points) and
    at each trace (e.g. 134 traces, so 134*4096 total different fiber profiles).
    """
    def __init__(self, pipeline_context):
        super(GenerateFiberProfileImage, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'fiber_profile'

    def do_stage(self, images):
        add_class_as_attribute(images, 'fiber_profile', FiberProfile)
        for image in images:
            normalized_fiber_profile_image = np.zeros_like(image.data)
            num_traces = image.trace.coefficients.shape[0]
            mean_of_horizontal_ranges = []
            full_width_half_maxes = []
            for trace_number in range(num_traces):
                horizontal_ranges = image.fiber_profile.horizontal_ranges[trace_number]
                normalized_coeffs = image.fiber_profile.fit_coefficients[trace_number]
                fwhm = fiber_profile_utils.Shapelets().full_width_half_max(normalized_coeffs[0])
                # interpolating coefficients into generating functions
                mean_of_horizontal_ranges.append(np.mean(horizontal_ranges, axis=1))
                coefficient_gen_funcs = fiber_profile_utils.interpolate_fiber_fits(normalized_coeffs,
                                                                                   mean_of_horizontal_ranges[-1])
                normalized_fiber_profile_values, indices = fiber_profile_utils.evaluate_normalized_fiber_profile_across_detector_x_y_per_trace(
                                                           image, coefficient_gen_funcs,
                                                           trace_number, window=int(7 * fwhm),
                                                           Model=self.fiber_profile_model,
                                                           renormalize=False)
                # += is important so that regions of high probability do not get accidentally overridden by the tail of
                # an adjacent profile.
                normalized_fiber_profile_image[indices] += normalized_fiber_profile_values
                full_width_half_maxes.append(fwhm)

            image.fiber_profile.normalized_fiber_profile_image = normalized_fiber_profile_image
            image.median_full_width_half_max = np.median(np.array(full_width_half_maxes))


class FiberProfileMaker(CalibrationMaker):
    def __init__(self, pipeline_context):
        super(FiberProfileMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'fiber_profile'

    @property
    def min_images(self):
        return 1

    def make_master_calibration_frame(self, images, image_config, logging_tags):
        single_image_list = [images[0]]
        logger.info('making master profile image on only the first image in the image list')
        SampleFiberProfileAcrossImage(self.pipeline_context).do_stage(single_image_list)
        GenerateFiberProfileImage(self.pipeline_context).do_stage(single_image_list)
        master_profile_image = Image(pipeline_context=self.pipeline_context,
                                     data=single_image_list[0].fiber_profile.normalized_fiber_profile_image, header=None)
        # also save all aspects of FiberProfile object as extra cards and the good regions we extracted. - but do
        # good region stuff in a different stage.

        return [master_profile_image]


class LoadFiberProfileImage(Stage):
    """
    This stage loads the fiber profile image (used for calculating optimal extraction weights) into each banzai Image.
    """
    def __init__(self, pipeline_context):
        super(LoadFiberProfileImage, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'fiber_profile'

    def do_stage(self, images):
        add_class_as_attribute(images, 'fiber_profile', FiberProfile)
        for image in images:
            normalized_fiber_profile_image = load_master_profile()
            image.fiber_profile.normalized_fiber_profile_image = normalized_fiber_profile_image


def horizontal_windows_to_fit_fiber_profile(image, trace_number):
    fitting_window_size = int(image.data.shape[1] / 26)
    # this is set so that for a 4k x 4k image, we get a window of roughly 150 points.
    smallest_good_x, largest_good_x = image.trace.high_signal_to_noise_region_bounds[trace_number]
    # the detector is too dim in some regions to fit.
    x_range = largest_good_x - smallest_good_x
    middle_intervals = 8
    wing_intervals = 4
    low_extents = np.concatenate((np.linspace(smallest_good_x, smallest_good_x + x_range / 3,
                                              num=wing_intervals, endpoint=False),
                                  np.linspace(smallest_good_x + x_range / 3,
                                              smallest_good_x + 2 * x_range / 3,
                                              num=middle_intervals, endpoint=False),
                                  np.linspace(smallest_good_x + 2 * x_range / 3,
                                              largest_good_x - fitting_window_size,
                                              num=wing_intervals)))
    low_extents = low_extents[low_extents.argsort()]
    high_extents = low_extents + fitting_window_size
    horizontal_ranges = list(zip(low_extents.astype(np.int), high_extents.astype(np.int)))
    return horizontal_ranges


if __name__ == "__main__":
    None

