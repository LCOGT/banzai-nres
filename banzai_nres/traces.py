"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import check_for_close_fit, check_flux_change, cross_correlate_image_indices, \
    optimize_coeffs_entire_lampflat_frame, fit_traces_order_by_order, get_number_of_lit_fibers
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute
from banzai_nres.images import Image
from banzai.stages import CalibrationMaker, Stage, MasterCalibrationDoesNotExist
from banzai.utils import fits_utils
from banzai import logs
from banzai import dbs
from astropy.io import fits
import os


logger = logs.get_logger(__name__)


class Trace(object):
    """
    Object for storing all the trace related attributes. This gets appended to each Image instance.
    """
    def __init__(self):
        self.coefficients = None
        self.fiber_order = None


class TraceMaker(CalibrationMaker):
    """
    Updates the master calibration trace file.
    """
    def __init__(self, pipeline_context):
        super(TraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.order_of_meta_fit = 6

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
        master_traces = make_master_traces(images, self, image_config, logging_tags,
                                           'global-meta', cross_correlate_num=1,
                                           order_of_meta_fit=self.order_of_meta_fit)

        return [master_traces]


class BlindTraceMaker(CalibrationMaker):
    """
    Fits traces order by order. Only use if you want to generate a new master
    trace file without loading any trace locations from the data-base
    :param average_trace_vertical_extent : should in no instance ever be changed unless the detector drastically
    changes. This is the approximate (good to within \pm 30 pixels) difference between the position of the bottom
    of the trace and its position when it contacts the edge of the detector. E.g. if you were to surround a trace in
    the minimum sized box possible, this is the y-height of the box (parallel to increasing order direction). Sign matters
    the convention is that if the traces curve upwards then this is positive. Negative if they curve downwards. Sign matters
    because this is exactly the guess for the second order coefficient of the blind trace-trace fit.
    """
    def __init__(self, pipeline_context):
        super(BlindTraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.order_of_meta_fit = 6
        self.order_of_poly_fit = 4
        self.average_trace_vertical_extent = 90  # do NOT haphazardly change this.

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'trace'

    @property
    def min_images(self):
        return 1

    def make_master_calibration_frame(self, images, image_config, logging_tags):
        master_traces = make_master_traces(images, self, image_config, logging_tags,
                                           'order-by-order', cross_correlate_num=1,
                                           order_of_poly_fits=self.order_of_poly_fit,
                                           order_of_meta_fit=self.order_of_meta_fit)
        return [master_traces]


class TraceUpdater(Stage):
    """
    Loads the most recent master trace file and stores it on the image under trace.coefficients and trace.fiber_order.
    Updates the trace centroid locations for a frame of any observation type.
    Will keep the as-imported master_trace locations if the reasonable criterion which indicate a good fit are not met.
    Right now, any updating of traces away from the recent master calibration file is not allowed. So technically,
    this stage only serves the purpose of loading the loading the master calibration traces. One can enable on the
    fly updating by allow_on_the_fly_updating = False, however that may pose issues downstream.
    """
    def __init__(self, pipeline_context):
        super(TraceUpdater, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.order_of_meta_fit = 6
        self.allow_on_the_fly_updating = False

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'trace'

    def do_stage(self, images):
        add_class_as_attribute(images, 'trace', Trace)
        for image in images:
            # getting coefficients from master trace file
            coefficients_and_indices_initial, fiber_order = get_trace_coefficients(image, self)
            num_lit_fibers = get_number_of_lit_fibers(coefficients_and_indices_initial)
            # once we write fiber_order correctly this should turn into len(fiber_order) or maybe assert
            # num_lit_fibers = len(fiber_order)
            image.trace.coefficients = coefficients_and_indices_initial
            image.trace.fiber_order = fiber_order

            # optimizing master traces on this frame in particular
            coefficients_and_indices_new = optimize_coeffs_entire_lampflat_frame(
                coefficients_and_indices_initial, image, num_of_lit_fibers=num_lit_fibers,
                order_of_meta_fit=self.order_of_meta_fit, bpm=None)
            logger.debug('refining trace coefficients on %s' % image.filename)

            close_fit = check_for_close_fit([coefficients_and_indices_new, coefficients_and_indices_initial],
                                            [image, image], num_lit_fibers, max_pixel_error=3)
            reasonable_flux_change = check_flux_change(coefficients_and_indices_new, coefficients_and_indices_initial,
                                                       image)
            # keeping the optimized traces only if they satisfy certain conditions
            if close_fit and reasonable_flux_change and self.allow_on_the_fly_updating:
                image.trace.coefficients = coefficients_and_indices_new
                logger.debug('New trace fit accepted on %s' % image.filename)
            if not close_fit or not reasonable_flux_change:
                logger.warning('Either 1. orders have possibly shifted drastically on %s, or it is very low S/N \n'
                               'resorting to as imported coefficients.' % image.filename)
        return images


def make_master_traces(images, maker_object, image_config, logging_tags, method,
                       cross_correlate_num=2, order_of_poly_fits=4, order_of_meta_fit=6):
    """
    :param images: List of banzai Image classes
    :param method: 'order-by-order' or 'global-meta'. Order by order should only be used when making a brand new Master
                if the current master traces have floated too far away.
    :param cross_correlate_num: number of frames to cross correlate trace solutions if there are sufficient number of
                frames. cross_correlate_num must be at least 1. If =1, no cross correlation will be done.
    :param maker_object: CalibrationMaker object.
    :return: Banzai image object where image.data are the trace coefficients. with order indices as the first column.
    """
    add_class_as_attribute(images, 'trace', Trace)
    master_trace_filename = maker_object.get_calibration_filename(image_config)
    logs.add_tag(logging_tags, method + 'master_trace', os.path.basename(master_trace_filename))

    satisfactory_fit = False
    image_indices_to_try, try_combinations_of_images = cross_correlate_image_indices(images, cross_correlate_num)
    coefficients_and_indices_list = []
    counter = 0
    while not satisfactory_fit:
        if try_combinations_of_images:
            images_to_try = [images[i] for i in image_indices_to_try[counter]]
        else:
            images_to_try = [images[counter]]
        counter += 1
        for image in images_to_try:
            if method == 'order-by-order':
                logger.debug('fitting order by order on %s' % image.filename)
                second_order_coefficient_guess = maker_object.average_trace_vertical_extent
                coefficients_and_indices_initial = fit_traces_order_by_order(image, second_order_coefficient_guess,
                                                                             order_of_poly_fits=order_of_poly_fits,
                                                                             num_lit_fibers=2)
                fiber_order = None
            if method == 'global-meta':
                logger.debug('importing master coeffs and refining fit on %s' % image.filename)
                coefficients_and_indices_initial, fiber_order = get_trace_coefficients(image, maker_object)

            num_lit_fibers = get_number_of_lit_fibers(coefficients_and_indices_initial)
            # once we write fiber_order correctly this should turn into len(fiber_order) or maybe assert
            # num_lit_fibers = len(fiber_order)

            coefficients_and_indices_list += [optimize_coeffs_entire_lampflat_frame(
                coefficients_and_indices_initial, image, num_of_lit_fibers=num_lit_fibers,
                order_of_meta_fit=order_of_meta_fit, bpm=None)]

        satisfactory_fit = check_for_close_fit(coefficients_and_indices_list, images_to_try, num_lit_fibers,
                                               max_pixel_error=1E-1)

        assert coefficients_and_indices_initial.shape == coefficients_and_indices_list[0].shape

        if not satisfactory_fit:
            logger.warning(
                "Unsatisfactory master fit. Trying a new set of lampflats. %s sets exist" % len(image_indices_to_try))

        if counter >= len(image_indices_to_try) and not satisfactory_fit:
            raise "No set of %s master meta-fits agree enough. Master Trace fitting failed on this batch%s ! " \
                  "Try generating order by order. This set of lamp flat frames possibly disagree." % (
                      cross_correlate_num, images_to_try[0].filename)

    header = fits_utils.create_master_calibration_header(images)
    header['FIBRORDR'] = fiber_order
    header['FITFLAT'] = images[0].filename
    header['OBSTYPE'] = 'TRACE'
    header['DATE-OBS'] = images[0].header.get('DATE-OBS')
    header['DAY-OBS'] = images[0].header.get('DAY-OBS')
    header['INSTRUME'] = images[0].header.get('TELESCOP')

    logger.info(os.path.basename(master_trace_filename))

    master_trace_coefficients = Image(pipeline_context=maker_object.pipeline_context,
                                      data=coefficients_and_indices_list[0], header=header)

    master_trace_coefficients.filename = master_trace_filename
    assert master_trace_coefficients.data.shape is not None

    return master_trace_coefficients


def get_trace_coefficients(image, maker_object):
    """
    :param image: Banzai Image
    :param maker_object: CalibrationMaker or Stage object - must have attribute pipeline context
    :return: The coefficients and indices (ndarray), and fiber order tuple from the nearest master trace file.
    """
    coefficients_and_indices, fiber_order = None, None
    master_trace_full_path = dbs.get_master_calibration_image(image, maker_object.calibration_type,
                                                              maker_object.group_by_keywords,
                                                              db_address=maker_object.pipeline_context.db_address)
    if image.header['OBSTYPE'] != 'TRACE' and os.path.isfile(master_trace_full_path):
        fiber_order = fits.getheader(master_trace_full_path).get('FIBRORDR')
        coefficients_and_indices = fits.getdata(master_trace_full_path)

        logger.debug('Imported master trace coefficients shape: ' + str(coefficients_and_indices.shape))
        assert coefficients_and_indices is not None
        assert fiber_order is not None

    if image.header['OBSTYPE'] != 'LAMPFLAT' and not os.path.isfile(master_trace_full_path):
        raise MasterCalibrationDoesNotExist

    return coefficients_and_indices, fiber_order
