"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import check_for_close_fit, cross_correlate_image_indices, \
    optimize_coeffs_entire_lampflat_frame, fit_traces_order_by_order, get_number_of_lit_fibers, Trace
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute
from banzai_nres.images import Image
from banzai.stages import CalibrationMaker, Stage, MasterCalibrationDoesNotExist
from banzai.utils import fits_utils
from banzai import logs
from banzai import dbs
from banzai.images import DataTable, regenerate_data_table_from_fits_hdu_list
from astropy.io import fits
import numpy as np
import os


logger = logs.get_logger(__name__)


class TraceMaker(CalibrationMaker):
    """
    Updates the master calibration trace file.
    """
    def __init__(self, pipeline_context):
        super(TraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.try_combinations_of_images = False
        self.cross_correlate_num = 1
        # if self.try_combinations_of_images, this is the number of images you compare to see if their trace fits agree

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'TRACE'

    @property
    def min_images(self):
        return 1

    def make_master_calibration_frame(self, images, image_config, logging_tags):
        """
        :param images:
        :param image_config:
        :param logging_tags:
        :return: frame with trace coefficients appended as the primary data object.
        This function cross correlates image.trace.coefficients between different fits (self.try_combinations_of_images)
        and saves the trace information for the set of images which agree the most.
        """
        master_trace_filename = self.get_calibration_filename(image_config)
        logs.add_tag(logging_tags, 'master_trace', os.path.basename(master_trace_filename))

        satisfactory_fit = False
        image_indices_to_try, try_combinations_of_images = cross_correlate_image_indices(images,
                                                                                         self.cross_correlate_num)
        counter = 0
        num_lit_fibers = get_number_of_lit_fibers(images[0])
        while not satisfactory_fit:
            if self.try_combinations_of_images:
                images_to_try = [images[i] for i in image_indices_to_try[counter]]
            else:
                images_to_try = [images[counter]]
            counter += 1

            coefficients_and_indices_list = [image.trace.coefficients for image in images_to_try]

            satisfactory_fit = check_for_close_fit(coefficients_and_indices_list, images_to_try, num_lit_fibers,
                                                   max_pixel_error=1E-1)

            if not satisfactory_fit:
                logger.warning(
                    "Too much disagreement between the master trace fits on this subset of lampflats"
                    ". Trying a new set of lampflats. {0} sets exist".format(len(image_indices_to_try)))

            if counter >= len(image_indices_to_try) and not satisfactory_fit:
                raise "No set of {0} master meta-fits agree enough. Master Trace fitting failed on this batch {1} ! " \
                      "Try generating order by order. This set of lamp flat frames possibly disagree." .format(
                          self.cross_correlate_num, images_to_try[0].filename)

        good_frame = images_to_try[0]

        header = fits_utils.create_master_calibration_header(images)
        header['FIBRORDR'] = good_frame.trace.fiber_order
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
        master_trace_calibration = Image(pipeline_context=self.pipeline_context,
                                         data=np.zeros((2, 2)), header=header, data_tables=master_cal_data_tables)

        master_trace_calibration.filename = master_trace_filename

        return [master_trace_calibration]


class TraceRefine(Stage):
    """
    Updates image.trace.coefficients via global-meta fitting. This procedure is robust enough that
    one can run this one arc or science frames if they so wished to refine their traces.
    :param max_number_of_images_to_refine : the number of images from the larger list images that you wish to actually fit.
    For instance if 1, we do one fit then we adopt that fit onto all other images in the list. Must be at least 1. Set to
    some large number if you want to always fit every frame in images.
    """
    def __init__(self, pipeline_context):
        super(TraceRefine, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.order_of_meta_fit = 6
        self.max_number_of_images_to_refine = 1
        self.refit_if_traces_shifted_substantially = True

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'trace'

    def do_stage(self, images):
        add_class_as_attribute(images, 'trace', Trace)
        for i, image in enumerate(images):
            num_lit_fibers = get_number_of_lit_fibers(image)

            if i < self.max_number_of_images_to_refine:
                logger.debug('refining with global meta on {0}'.format(image.filename))
                refined_trace_coefficients = optimize_coeffs_entire_lampflat_frame(
                    image.trace.coefficients, image, num_of_lit_fibers=num_lit_fibers,
                    order_of_meta_fit=self.order_of_meta_fit, bpm=image.bpm)
                fiber_order, refined_trace_coefficients = self.refit_if_necessary(image,
                                                                                  num_lit_fibers,
                                                                                  refined_trace_coefficients,
                                                                                  absolute_pixel_tolerance=1.0)
            else:
                logger.debug('adopting last global-meta fit onto {0}'.format(image.filename))
                refined_trace_coefficients = images[self.max_number_of_images_to_refine - 1].trace.coefficients
                fiber_order = images[self.max_number_of_images_to_refine - 1].trace.fiber_order

            image.trace.coefficients = refined_trace_coefficients
            image.trace.fiber_order = fiber_order
        return images

    def refit_if_necessary(self, image, num_lit_fibers, refined_trace_coefficients, absolute_pixel_tolerance=1.0):
        # consider using an independent function to see if traces have shifted. However this may be bad
        # because getting under a pixel accuracy with just a peak finder is quite hard. Might just be worth it
        # to use the algorithm which we have already developed which does so with extremely high accuracy.
        """
        If the new coefficients differ by more than absolute_pixel_tolerance on average compared to the
        as loaded coefficients, then we
        refit the traces on image order-by-order then redoes the global-meta on those new coefficients.

        Note that because this can come after the order-by-order fitting, too small of a pixel tolerance will
        just cause us to redo fits every time. A pixel tolerance of 1 is good because standard NRES shifts are on the
        order of less than a pixel.
        """
        coefficients_to_compare = [refined_trace_coefficients, image.trace.coefficients]
        traces_did_not_shift_substantially = check_for_close_fit(coefficients_to_compare,
                                                                 [image, image], num_lit_fibers,
                                                                 max_pixel_error=absolute_pixel_tolerance)

        if self.refit_if_traces_shifted_substantially and (not traces_did_not_shift_substantially):
            logger.warning('traces shifted beyond allowable amount, refitting lampflat order-by-order.')
            blind_fit_maker = GenerateInitialGuessForTraceFitFromScratch(pipeline_context=self.pipeline_context)
            image = blind_fit_maker.do_stage([image])[0]
            refined_trace_coefficients = optimize_coeffs_entire_lampflat_frame(
                image.trace.coefficients, image, num_of_lit_fibers=num_lit_fibers,
                order_of_meta_fit=self.order_of_meta_fit, bpm=image.bpm)
            fiber_order = None
        else:
            fiber_order = image.trace.fiber_order
        return fiber_order, refined_trace_coefficients


class GenerateInitialGuessForTraceFitFromScratch(Stage):
    """
    Generates an initial guess for the trace global-meta fitting by fitting the traces order by order.
    :param average_trace_vertical_extent : should in no instance ever be changed unless the detector drastically
    changes. This should be the approximate (good to within \pm 30 pixels) difference between the position of the bottom
    of the trace and its position when it contacts the edge of the detector. E.g. if you were to surround a trace in
    the minimum sized box possible, this is the y-height of the box (parallel to increasing order direction). Sign matters
    the convention is that if the traces curve upwards then this is positive. Negative if they curve downwards. Sign matters
    because this is the guess for the second order coefficient of the blind order-by-order trace fit.

    :param max_number_of_images_to_fit : the number of images from the larger list images that you wish to actually fit.
    For instance if 1, we do one fit then we adopt that fit onto all other images in the list.
    """
    def __init__(self, pipeline_context):
        super(GenerateInitialGuessForTraceFitFromScratch, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.order_of_poly_fit = 4
        self.average_trace_vertical_extent = 90  # do NOT haphazardly change this.
        self.max_number_of_images_to_fit = 1

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'trace'

    def do_stage(self, images):
        add_class_as_attribute(images, 'trace', Trace)
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


class LoadInitialGuessForTraceFit(Stage):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    def __init__(self, pipeline_context):
        super(LoadInitialGuessForTraceFit, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'trace'

    def do_stage(self, images):
        add_class_as_attribute(images, 'trace', Trace)
        for image in images:
            logger.debug('importing master coeffs from %s' % image.filename)
            coefficients_and_indices_initial, fiber_order = self.get_trace_coefficients(image)
            image.trace.fiber_order = fiber_order
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
        if image.header['OBSTYPE'] != 'TRACE' and os.path.isfile(master_trace_full_path):
            fiber_order = fits.getheader(master_trace_full_path).get('FIBRORDR')
            hdu_list = fits.open(master_trace_full_path)
            coeffs_name = Trace().coefficients_table_name
            dict_of_table = regenerate_data_table_from_fits_hdu_list(hdu_list,table_extension_name=coeffs_name)
            coefficients_and_indices_table = dict_of_table[coeffs_name]
            coefficients_and_indices, loaded_fiber_order = Trace().convert_astropy_table_coefficients_to_numpy_array(
                                                                                        coefficients_and_indices_table)

            assert coefficients_and_indices is not None
            logger.debug('Imported master trace coefficients shape: ' + str(coefficients_and_indices.shape))
            assert fiber_order == loaded_fiber_order

        if image.header['OBSTYPE'] != 'LAMPFLAT' and not os.path.isfile(master_trace_full_path):
            raise MasterCalibrationDoesNotExist

        return coefficients_and_indices, fiber_order
