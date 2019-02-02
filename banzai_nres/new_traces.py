"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import fit_traces_order_by_order, Trace
from banzai_nres.utils import trace_utils
from banzai_nres.images import NRESImage
from banzai.stages import Stage
from banzai.calibrations import CalibrationMaker
from banzai.utils import fits_utils
from banzai import dbs
from banzai.images import DataTable, regenerate_data_table_from_fits_hdu_list
from astropy.io import fits
import numpy as np
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

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        traces = []
        for image in images:
            trace = Trace()
            logger.debug('fitting traces order by order', image=image)
            trace.coefficients = fit_traces_order_by_order(image.data,
                                                                 self.second_order_coefficient_guess,
                                                                 order_of_poly_fits=self.order_of_poly_fit,
                                                                 num_lit_fibers=image.num_lit_fibers())
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
            image.trace = Trace()
            image.trace.coefficients = self.get_trace_coefficients(image)

            if image.trace.coefficients is None:
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
