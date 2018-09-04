"""
Module for identifying the regions on the CCD where the traces have a S/N above some threshold.

Author
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from banzai.stages import Stage, MasterCalibrationDoesNotExist
from banzai_nres.utils.good_regions_utils import identify_high_SN_region_bounds_along_traces,\
     flag_traces_with_insufficient_high_SN_region


class IdentifyRegionsWithGoodSignaltoNoise(Stage):
    """
    This stage identifies the x-bounds per trace where the S/N exceeds the given threshold value.
    Thus the x-bounds give the column indices per trace where the data is clean enough to fit with a fiber profile.
    """
    def __init__(self, pipeline_context):
        super(IdentifyRegionsWithGoodSignaltoNoise, self).__init__(pipeline_context)
        self.signal_to_noise_threshold = 22

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'fiber_profile'

    def do_stage(self, images):
        for image in images:
            image.trace.high_signal_to_noise_region_bounds = identify_high_SN_region_bounds_along_traces(image,
                                                             sig_to_noise_threshold=self.signal_to_noise_threshold)

            min_good_region_width = int(image.data.shape[1] / 10) # roughly 400 pixels.
            image.trace.has_sufficient_signal_to_noise = flag_traces_with_insufficient_high_SN_region(image,
                                                                                      min_good_region_width)
        return images
