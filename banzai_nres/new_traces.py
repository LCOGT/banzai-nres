"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.new_trace_utils import Trace
from banzai_nres.utils import trace_utils
from banzai.stages import Stage
from banzai.calibrations import CalibrationMaker
from banzai import dbs
from banzai.images import DataTable, regenerate_data_table_from_fits_hdu_list
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
        self.fit_march_parameters = {'window': 100, 'step_size': 6}
        self.match_filter_parameters = {'min_peak_spacing': 5, 'neighboring_peak_flux_ratio': 20}

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        traces = []
        for image in images:
            logger.debug('fitting traces order by order', image=image)
            trace = Trace.fit_traces(image, poly_fit_order=self.order_of_poly_fit,
                                     second_order_coefficient_guess=self.second_order_coefficient_guess,
                                     fit_march_parameters= self.fit_march_parameters,
                                     match_filter_parameters=self.match_filter_parameters)
            traces.append(trace)
        return traces


class LoadTrace(Stage):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    def __init__(self, pipeline_context):
        super(LoadTrace, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
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
            #NOTE this will fail if master_trace_full_path does not exist or leads to a none file. is there
            # a banzai util which will return None which I can use instead of fits.open?
            hdu_list = fits.open(master_trace_path)
            image.trace = Trace.load(hdu_list)

            if image.trace is None:
                logger.error("Can't find trace coefficients for image, stopping reduction", image=image)
                images_to_remove.append(image)
                continue
        for image in images_to_remove:
            images.remove(image)
        return images

    def get_trace_coefficients(self, image):
        """
        :param image: Banzai Image
        :return: The coefficients and indices (ndarray) array. The first column are the diffraction
        order designations (e.g. 0,1,2...67)
        """
        coefficients_and_indices = None
        master_trace_full_path = dbs.get_master_calibration_image(image, self.calibration_type,
                                                                  self.master_selection_criteria,
                                                                  db_address=self.pipeline_context.db_address)

        if master_trace_full_path is None or not os.path.exists(master_trace_full_path):
            logger.error('Master trace fit file not found for {0}.'.format(image.filename))

        else:
            hdu_list = fits.open(master_trace_full_path)
            coeffs_name = Trace().coefficients_table_name
            dict_of_table = regenerate_data_table_from_fits_hdu_list(hdu_list, table_extension_name=coeffs_name)
            coefficients_and_indices_table = dict_of_table[coeffs_name]
            coefficients_and_indices, lit_fibers = trace_utils.convert_astropy_table_coefficients_to_numpy_array(
                                                                                        coefficients_and_indices_table,
                                                                                        coefficients_table_name=coeffs_name)
            logger.info('Imported master trace coefficients array with '
                        'shape {0}'.format(str(coefficients_and_indices.shape)))

        return coefficients_and_indices
