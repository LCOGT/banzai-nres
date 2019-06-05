"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import Trace, AllTraceFitter
from banzai.calibrations import Stage, ApplyCalibration, create_master_calibration_header
from banzai.utils import file_utils

import banzai_nres.settings as nres_settings
import banzai.settings as banzai_settings

import sep
import os
import logging

logger = logging.getLogger(__name__)


class TraceMaker(Stage):
    def __init__(self, runtime_context):
        super(TraceMaker, self).__init__(runtime_context)
        self.runtime_context = runtime_context
        self.order_of_poly_fit = nres_settings.TRACE_FIT_POLYNOMIAL_ORDER
        self.second_order_coefficient_guess = nres_settings.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS
        self.trace_table_name = nres_settings.TRACE_TABLE_NAME
        self.xmin = nres_settings.WINDOW_FOR_TRACE_IDENTIFICATION['min']
        self.xmax = nres_settings.WINDOW_FOR_TRACE_IDENTIFICATION['max']
        self.min_peak_to_peak_spacing = nres_settings.MIN_FIBER_TO_FIBER_SPACING
        self.min_snr = nres_settings.MIN_SNR_FOR_TRACE_IDENTIFICATION

    @property
    def calibration_type(self):
        return 'TRACE'

    def do_stage(self, image):
        master_header = create_master_calibration_header(old_header=image.header, images=[image])
        master_header['OBSTYPE'] = self.calibration_type
        master_filename = banzai_settings.CALIBRATION_FILENAME_FUNCTIONS[self.calibration_type](image)
        master_filepath = self._get_filepath(self.runtime_context, image, master_filename)
        logger.info('fitting traces order by order', image=image)

        bkg_subtracted_image_data = image.data - sep.Background(image.data).back()
        fitter = AllTraceFitter(xmin=self.xmin, xmax=self.xmax,
                                min_peak_to_peak_spacing=self.min_peak_to_peak_spacing,
                                min_snr=self.min_snr)
        trace = Trace(data=None, filepath=master_filepath, header=master_header, image=image,
                      num_centers_per_trace=image.data.shape[1], table_name=self.trace_table_name,
                      obstype=self.calibration_type)
        trace = fitter.fit_traces(trace=trace, image_data=bkg_subtracted_image_data,
                                  poly_fit_order=self.order_of_poly_fit,
                                  second_order_coefficient_guess=self.second_order_coefficient_guess,
                                  image_noise_estimate=image.header['RDNOISE'])
        logger.info('Created master trace', image=image, extra_tags={'calibration_type': self.calibration_type,
                                                                     'output_path': master_filepath,
                                                                     'calibration_obstype': master_header['OBSTYPE']})
        return trace

    @staticmethod
    def _get_filepath(runtime_context, lampflat_image, master_filename):
        output_directory = file_utils.make_output_directory(runtime_context, lampflat_image)
        return os.path.join(output_directory, os.path.basename(master_filename))


class LoadTrace(ApplyCalibration):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    def __init__(self, runtime_context):
        super(LoadTrace, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'TRACE'

    def apply_master_calibration(self, image, master_calibration_path):
        image.trace = Trace.load(master_calibration_path, extension_name=nres_settings.TRACE_TABLE_NAME)
        master_trace_filename = os.path.basename(master_calibration_path)
        image.header['L1IDTRAC'] = (master_trace_filename, 'ID of trace centers file')
        logger.info('Loading trace centers', image=image,  extra_tags={'L1IDTRAC': image.header['L1IDTRAC']})
        return image

    def do_stage(self, image):
        master_calibration_path = self.get_calibration_filename(image)
        if master_calibration_path is None:
            self.on_missing_master_calibration(image)
            return image
        return self.apply_master_calibration(image, master_calibration_path)
