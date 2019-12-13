"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import Trace, AllTraceFitter
from banzai.calibrations import CalibrationMaker, CalibrationUser
from banzai.utils import file_utils

import sep
import os
import logging

logger = logging.getLogger('banzai')


class TraceMaker(CalibrationMaker):
    def __init__(self, runtime_context):
        super(TraceMaker, self).__init__(runtime_context)
        self.runtime_context = runtime_context
        self.order_of_poly_fit = runtime_context.TRACE_FIT_POLYNOMIAL_ORDER
        self.second_order_coefficient_guess = runtime_context.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS
        self.trace_table_name = runtime_context.TRACE_TABLE_NAME
        self.xmin = runtime_context.WINDOW_FOR_TRACE_IDENTIFICATION['min']
        self.xmax = runtime_context.WINDOW_FOR_TRACE_IDENTIFICATION['max']
        self.min_peak_to_peak_spacing = runtime_context.MIN_FIBER_TO_FIBER_SPACING
        self.min_snr = runtime_context.MIN_SNR_FOR_TRACE_IDENTIFICATION

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        image = images[0]
        master_header = create_master_calibration_header(old_header=image.header, images=[image])
        master_header['OBSTYPE'] = self.calibration_type
        make_calibration_name = file_utils.make_calibration_filename_function(self.calibration_type,
                                                                              self.runtime_context)

        master_filename = make_calibration_name(images[0])
        logger.info('fitting traces order by order', image=image)
        bkg_subtracted_image_data = image.data - sep.Background(image.data).back()
        fitter = AllTraceFitter(xmin=self.xmin, xmax=self.xmax,
                                min_peak_to_peak_spacing=self.min_peak_to_peak_spacing,
                                min_snr=self.min_snr)
        calibration_set_criteria = self.runtime_context.CALIBRATION_SET_CRITERIA.get(self.calibration_type, {})
        trace = Trace(data=None, filepath=master_filename, header=master_header, image=image,
                      num_centers_per_trace=image.data.shape[1], trace_table_name=self.trace_table_name,
                      obstype=self.calibration_type, calibration_set_criteria=calibration_set_criteria)
        trace = fitter.fit_traces(trace=trace, image_data=bkg_subtracted_image_data,
                                  poly_fit_order=self.order_of_poly_fit,
                                  second_order_coefficient_guess=self.second_order_coefficient_guess,
                                  image_noise_estimate=image.header['RDNOISE'])
        logger.info('Created master trace', image=image, extra_tags={'calibration_type': self.calibration_type,
                                                                     'master_filename': master_filename,
                                                                     'calibration_obstype': master_header['OBSTYPE']})
        return trace


class LoadTrace(CalibrationUser):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    def __init__(self, runtime_context):
        super(LoadTrace, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'TRACE'

    def apply_master_calibration(self, image, master_calibration_path):
        image.trace = Trace.load(master_calibration_path, trace_table_name=self.runtime_context.TRACE_TABLE_NAME)
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
