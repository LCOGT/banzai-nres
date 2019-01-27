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


class SaveTrace(CalibrationMaker):
    """
    Updates the master calibration trace file.

    Note:  if self.try_combinations_of_images, cross_correlate_num is the
    number of images you compare to see if their trace fits agree
    """
    def __init__(self, pipeline_context):
        super(SaveTrace, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.master_selection_criteria = self.pipeline_context.CALIBRATION_SET_CRITERIA.get(
            self.calibration_type.upper(), [])
        self.fiber_order_header_name = 'FIBRORDR'

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        """
        :param images: list of Banzai-NRES Image objects.
        :return: frame with trace coefficients appended as astropy tables etc.
        """
        good_frame = images[0]
        num_lit_fibers = good_frame.num_lit_fibers()

        make_calibration_name = self.pipeline_context.CALIBRATION_FILENAME_FUNCTIONS[self.calibration_type]
        master_trace_filename = make_calibration_name(good_frame)

        logger.debug('master_trace', os.path.basename(master_trace_filename))

        header = fits_utils.create_master_calibration_header(images)
        header[self.fiber_order_header_name] = good_frame.trace.fiber_order
        header['FITFLAT'] = good_frame.filename
        header['OBSTYPE'] = 'TRACE'
        header['DATE-OBS'] = good_frame.header.get('DATE-OBS')
        header['DAY-OBS'] = good_frame.header.get('DAY-OBS')
        header['INSTRUME'] = good_frame.header.get('TELESCOP')
        header['OBJECTS'] = good_frame.header.get('OBJECTS')
        logger.debug('master calibration filename in TraceMaker is {0}'.format(os.path.basename(master_trace_filename)))

        master_trace_coefficients_table = trace_utils.convert_numpy_array_coefficients_to_astropy_table(num_lit_fibers,
                                                                                                        fiber_order=good_frame.trace.fiber_order)
        trace_centroids = good_frame.trace.get_trace_centers()[0]
        master_trace_centroids_table = trace_utils.convert_numpy_array_trace_centroids_to_astropy_table(num_lit_fibers,
                                                                                                        trace_centroids,
                                                                                                        good_frame.trace.coefficients,
                                                                                                        good_frame.trace.fiber_order)
        center_name = good_frame.trace.trace_center_table_name
        coefficients_name = good_frame.trace.coefficients_table_name

        master_trace_centroids_table[center_name].description = 'The y position of the trace centers for each x pixel' \
                                                                ' from 0 to the number of columns in the image'
        master_trace_centroids_table[center_name].unit = 'pixel'
        master_trace_coefficients_table[coefficients_name].description = 'Coefficients for the ' \
        '0th, 1st,..., Nth order legendre polynomials which describe the trace center across the detector. The ' \
        'Legendre polynomials are normalized from -1 to 1 over the range 0,..., 4096=number of columns of image.'

        master_trace_coefficients_table = DataTable(data_table=master_trace_coefficients_table, name=coefficients_name)
        master_trace_centroids_table = DataTable(data_table=master_trace_centroids_table, name=center_name)
        master_cal_data_tables = {coefficients_name: master_trace_coefficients_table,
                                  center_name: master_trace_centroids_table}
        master_trace_calibration = NRESImage(pipeline_context=self.pipeline_context,
                                             data=np.zeros((2, 2)), header=header, data_tables=master_cal_data_tables)

        master_trace_calibration.filename = master_trace_filename
        return master_trace_calibration


class FitTrace(Stage):
    """
    Loads trace coefficients from file and appends them onto the image object.
    If no master file is found or self.always_generate_traces_from_scratch, then it will do a blind fit:

    Generates an initial guess for the trace global-meta fitting by fitting the traces order by order.
    :param second_order_coefficient_guess : should in no instance ever be changed unless the detector drastically
    changes. This should be the approximate (good to within \pm 30 pixels) difference between the position of the bottom
    of the trace and its position when it contacts the edge of the detector. E.g. if you were to surround a trace in
    the minimum sized box possible, this is the y-height of the box (parallel to increasing order direction). Sign matters
    the convention is that if the traces curve upwards then this is positive. Negative if they curve downwards. Sign matters
    because this is the guess for the second order coefficient of the blind order-by-order trace fit.
    """
    def __init__(self, pipeline_context):
        super(FitTrace, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.master_selection_criteria = self.pipeline_context.CALIBRATION_SET_CRITERIA.get(
            self.calibration_type.upper(), [])
        self.order_of_poly_fit = 4
        self.second_order_coefficient_guess = self.pipeline_context.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS

    @property
    def calibration_type(self):
        return 'TRACE'

    def do_stage(self, images):
        for image in images:
            image.trace = Trace()
            logger.debug('fitting order by order on {0}'.format(image.filename))
            image.trace.coefficients = fit_traces_order_by_order(image.data,
                                                                 self.second_order_coefficient_guess,
                                                                 order_of_poly_fits=self.order_of_poly_fit,
                                                                 num_lit_fibers=image.num_lit_fibers())
        return images


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
        for image in images:
            image.trace = Trace()
            image.trace.coefficients = self.get_trace_coefficients(image)
        return images

    def get_trace_coefficients(self, image):
        """
        :param image: Banzai Image
        :return: The coefficients and indices (ndarray), and fiber order tuple from the nearest master trace file.
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
                                                                                        coefficients_and_indices_table)

            assert coefficients_and_indices is not None
            logger.info('Imported master trace coefficients array with '
                        'shape {0}'.format(str(coefficients_and_indices.shape)))

        return coefficients_and_indices
