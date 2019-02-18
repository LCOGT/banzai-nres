"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import Trace, AllTraceFitter
from banzai.stages import Stage
from banzai.calibrations import CalibrationMaker
from banzai import dbs
from astropy.io import fits
import os
import logging

logger = logging.getLogger(__name__)


class TraceMaker(CalibrationMaker):
    def __init__(self, pipeline_context):
        super(TraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.master_selection_criteria = self.pipeline_context.CALIBRATION_SET_CRITERIA.get(
            self.calibration_type.upper(), [])
        self.order_of_poly_fit = 4
        self.second_order_coefficient_guess = self.pipeline_context.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS
        self.trace_table_name = self.pipeline_context.TRACE_TABLE_NAME
        self.fit_march_parameters = {'window': 100, 'step_size': 6}
        self.match_filter_parameters = {'min_peak_spacing': 5, 'neighboring_peak_flux_ratio': 5}

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        traces = []
        for image in images:
            logger.debug('fitting traces order by order', image=image)
            fitter = AllTraceFitter(march_parameters=self.fit_march_parameters)
            trace = fitter.fit_traces(cls=Trace, image=image, poly_fit_order=self.order_of_poly_fit,
                                      second_order_coefficient_guess=self.second_order_coefficient_guess,
                                      match_filter_parameters=self.match_filter_parameters)
            trace.trace_table_name = self.trace_table_name
            traces.append(trace)
        return traces


class LoadTrace(Stage):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    def __init__(self, pipeline_context):
        super(LoadTrace, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.trace_table_name = self.pipeline_context.TRACE_TABLE_NAME
        self.master_selection_criteria = self.pipeline_context.CALIBRATION_SET_CRITERIA.get(
            self.calibration_type.upper(), [])

    @property
    def calibration_type(self):
        return 'TRACE'

    def do_stage(self, images):
        images_to_remove = []
        for image in images:
            logger.info('loading trace data onto image', image=image)
            master_trace_path = dbs.get_master_calibration_image(image, self.calibration_type,
                                                                 self.master_selection_criteria,
                                                                 db_address=self.pipeline_context.db_address)

            if master_trace_path is None or not os.path.exists(master_trace_path):
                logger.error('Master trace fit file not found for {0}.'.format(image.filename))
                logger.error("Can't find trace coefficients for image, stopping reduction", image=image)
                images_to_remove.append(image)
                continue
            else:
                hdu_list = fits.open(master_trace_path)
                image.trace = Trace.load(hdu_list, trace_table_name=self.trace_table_name)
        for image in images_to_remove:
            images.remove(image)
        return images
