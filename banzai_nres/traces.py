"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import fit_traces_order_by_order, get_number_of_lit_fibers, Trace
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute
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
    """
    Updates the master calibration trace file.

    Note:  if self.try_combinations_of_images, cross_correlate_num is the
    number of images you compare to see if their trace fits agree
    """
    def __init__(self, pipeline_context):
        super(TraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.fiber_order_header_name = 'FIBRORDR'

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        """
        :param images: list of Banzai-NRES Image objects.
        :return: frame with trace coefficients appended as astropy tables etc.
        """
        if all(image.trace.coefficients is None for image in images):
            logger.error('Each image in images are missing trace coefficients, aborting master trace calibration saving')
            return []

        good_frame = images[0]
        num_lit_fibers = get_number_of_lit_fibers(images[0])

        logger.debug(str(self.pipeline_context.FRAME_CLASS))
        logger.debug(self.pipeline_context.CALIBRATION_FILENAME_FUNCTIONS)
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

        master_trace_coefficients_table = good_frame.trace.convert_numpy_array_coefficients_to_astropy_table(num_lit_fibers,
                                                                                                             fiber_order=good_frame.trace.fiber_order)
        trace_centroids = good_frame.trace.get_trace_centroids_from_coefficients(image_width=good_frame.data.shape[1])[0]
        master_trace_centroids_table = good_frame.trace.convert_numpy_array_trace_centroids_to_astropy_table(num_lit_fibers,
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
        return [master_trace_calibration]


class InitialTraceFit(Stage):
    """
    Loads trace coefficients from file and appends them onto the image object.
    If no master file is found or self.always_generate_traces_from_scratch, then it will do a blind fit:

    Generates an initial guess for the trace global-meta fitting by fitting the traces order by order.
    :param average_trace_vertical_extent : should in no instance ever be changed unless the detector drastically
    changes. This should be the approximate (good to within \pm 30 pixels) difference between the position of the bottom
    of the trace and its position when it contacts the edge of the detector. E.g. if you were to surround a trace in
    the minimum sized box possible, this is the y-height of the box (parallel to increasing order direction). Sign matters
    the convention is that if the traces curve upwards then this is positive. Negative if they curve downwards. Sign matters
    because this is the guess for the second order coefficient of the blind order-by-order trace fit.

    :param max_number_of_images_to_fit : the number of images from the larger list images that you wish to actually fit.
    For instance if 1, we do one fit then we adopt that fit onto all other images in the list. Set to some large number
    if you want to do a fit onto every image in the stack.
    """
    def __init__(self, pipeline_context):
        super(InitialTraceFit, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.always_generate_traces_from_scratch = False
        self.order_of_poly_fit = 4
        self.average_trace_vertical_extent = 90  # do NOT haphazardly change this.
        self.max_number_of_images_to_fit = 1

    @property
    def calibration_type(self):
        return 'TRACE'

    def do_stage(self, images):
        add_class_as_attribute(images, 'trace', Trace)
        if self.always_generate_traces_from_scratch:
            images = self.blind_fit_traces_on_images(images)
        if not self.always_generate_traces_from_scratch:
            for image in images:
                logger.debug('importing master coeffs from %s' % image.filename)
                coefficients_and_indices_initial, fiber_order = self.get_trace_coefficients(image)
                image.trace.fiber_order = fiber_order
                image.trace.coefficients = coefficients_and_indices_initial
                if image.trace.coefficients is None:
                    logger.error('No master trace file found for {0}'.format(image.filename))

        if any(image.trace.coefficients is None for image in images):
            logger.error('No master calibration found for at least one image, refitting all images')
            images = self.blind_fit_traces_on_images(images)
        return images

    def blind_fit_traces_on_images(self, images):
        for i, image in enumerate(images):
            num_lit_fibers = get_number_of_lit_fibers(image)
            if i < self.max_number_of_images_to_fit:
                logger.debug('fitting order by order on {0}'.format(image.filename))
                second_order_coefficient_guess = self.average_trace_vertical_extent
                coefficients_and_indices_initial = fit_traces_order_by_order(image, second_order_coefficient_guess,
                                                                             order_of_poly_fits=self.order_of_poly_fit,
                                                                             num_lit_fibers=num_lit_fibers)
            else:
                logger.debug('adopting last order-order fit onto {0}'.format(image.filename))
                coefficients_and_indices_initial = images[self.max_number_of_images_to_fit - 1].trace.coefficients
            image.trace.fiber_order = None
            image.trace.coefficients = coefficients_and_indices_initial
        return images

    def get_trace_coefficients(self, image):
        """
        :param image: Banzai Image
        :return: The coefficients and indices (ndarray), and fiber order tuple from the nearest master trace file.
        """
        coefficients_and_indices, fiber_order = None, None
        master_trace_full_path = dbs.get_master_calibration_image(image, self.calibration_type,
                                                                  self.group_by_keywords,
                                                                  db_address=self.pipeline_context.db_address)

        if master_trace_full_path is None or not os.path.exists(master_trace_full_path):
            logger.error('Master trace fit file not found, will '
                         'attempt a blind fit on file {0}.'.format(image.filename))

        else:
            fiber_order_header_name = TraceMaker(self.pipeline_context).fiber_order_header_name
            fiber_order = fits.getheader(master_trace_full_path).get(fiber_order_header_name)
            hdu_list = fits.open(master_trace_full_path)
            coeffs_name = Trace().coefficients_table_name
            dict_of_table = regenerate_data_table_from_fits_hdu_list(hdu_list, table_extension_name=coeffs_name)
            coefficients_and_indices_table = dict_of_table[coeffs_name]
            coefficients_and_indices, loaded_fiber_order = Trace().convert_astropy_table_coefficients_to_numpy_array(
                                                                                        coefficients_and_indices_table)

            assert coefficients_and_indices is not None
            logger.info('Imported master trace coefficients array with '
                        'shape {0}'.format(str(coefficients_and_indices.shape)))
            assert fiber_order == loaded_fiber_order

        return coefficients_and_indices, fiber_order
