"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import check_for_close_fit, check_flux_change, cross_correlate_image_indices, \
    optimize_coeffs_entire_lampflat_frame, fit_traces_order_by_order
from banzai_nres.dbs import get_trace_coefficients
from banzai_nres.images import Image
from banzai.stages import CalibrationMaker, Stage, MasterCalibrationDoesNotExist
from banzai.utils import fits_utils
from banzai import logs
from astropy.io import fits
import os


logger = logs.get_logger(__name__)


class TraceMaker(CalibrationMaker):
    """
    Updates the master calibration trace file.
    """
    def __init__(self, pipeline_context):
        super(TraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context

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
        master_bias_filename = self.get_calibration_filename(image_config)
        master_traces = make_master_traces(images, self, image_config, logging_tags,
                                           'global-meta', cross_correlate_num=1)

        return [master_traces]


class BlindTraceMaker(CalibrationMaker):
    """
    Fits traces order by order. Only use if you want to generate a new master
    trace file without loading any trace locations from the data-base
    """
    def __init__(self, pipeline_context):
        super(BlindTraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context

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
                                           'order-by-order', cross_correlate_num=1)
        logger.info(master_traces.data.shape)
        return [master_traces]


class TraceUpdater(Stage):
    """
    Updates the trace centroid locations for a frame
    of any observation type. Will keep the as-imported master_trace locations
    if the reasonable criterion which indicate a good fit are not met.
    """
    def __init__(self, pipeline_context):
        super(TraceUpdater, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'trace'

    def do_stage(self, images):
        for image in images:
            coefficients_and_indices_initial, fiber_order = get_trace_coefficients(image, self.pipeline_context)
            assert image.fiber_order == fiber_order
            coefficients_and_indices_new = optimize_coeffs_entire_lampflat_frame(
                coefficients_and_indices_initial, image, order_of_meta_fit=6)
            logger.info('refining trace coefficients on %s' % image.filename)

            close_fit = check_for_close_fit([coefficients_and_indices_new, coefficients_and_indices_initial],
                                            [image, image], max_pixel_error=3)
            reasonable_flux_change = check_flux_change(coefficients_and_indices_new, coefficients_and_indices_initial,
                                                       image)

            if close_fit and reasonable_flux_change:
                image.trace_fit_coefficients = coefficients_and_indices_new
                logger.info('New trace fit accepted on %s' % image.filename)
            if not close_fit or not reasonable_flux_change:
                logger.info('Either 1. orders have possibly shifted drastically on %s, or it is very low S/N \n'
                            'resorting to as imported coefficients.' % image.filename)
        return images


def make_master_traces(images, maker_object, image_config, logging_tags, method, cross_correlate_num=2):
    """
    :param images: List of banzai Image classes
    :param method: 'order-by-order' or 'global-meta'. Order by order should only be used when making a brand new Master
                if the current master traces have floated too far away.
    :param cross_correlate_num: number of frames to cross correlate trace solutions if there are sufficient number of
                frames. cross_correlate_num must be at least 1. If =1, no cross correlation will be done.
    :param maker_object: CalibrationMaker object.
    :return: Banzai image object where image.data are the trace coefficients. with order indices as the first column.
    """
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
                logger.info('fitting order by order on %s' % image.filename)
                coefficients_and_indices_initial, fiber_order = fit_traces_order_by_order(image, order_of_poly_fits=4)
            if method == 'global-meta':
                logger.info('importing master coeffs and refining fit on %s' % image.filename)
                coefficients_and_indices_initial, fiber_order = get_trace_coefficients(image)

            coefficients_and_indices_list += [optimize_coeffs_entire_lampflat_frame(
                coefficients_and_indices_initial, image, order_of_meta_fit=6)]

        satisfactory_fit = check_for_close_fit(coefficients_and_indices_list, images_to_try, max_pixel_error=1E-1)

        if not satisfactory_fit:
            logger.warning(
                "Unsatisfactory master fit. Trying a new set of lampflats. %s sets exist" % len(image_indices_to_try))

        if counter >= len(image_indices_to_try) and not satisfactory_fit:
            raise "No set of %s master meta-fits agree enough. Master Trace fitting failed on this batch%s ! " \
                  "Try generating order by order. This set of lamp flat frames possibly disagree." % (
                      cross_correlate_num, images_to_try[0].filename)

    # choice is arbitrary since they all agree.
    header = fits_utils.create_master_calibration_header(images)
    header['FIBRORDR'] = fiber_order
    header['FITFLAT'] = images[0].filename
    header['OBSTYPE'] = 'TRACE'
    header['DATE-OBS'] = images[0].header['DATE-OBS']
    header['DAY-OBS'] = images[0].header['DAY-OBS']
    header['INSTRUME'] = images[0].header['TELESCOP']

    logger.info(os.path.basename(master_trace_filename))

    master_trace_coefficients = Image(pipeline_context=maker_object.pipeline_context,
                                      data=coefficients_and_indices_list[0], header=header)

    master_trace_coefficients.filename = master_trace_filename
    logger.info('found coefficients after image call:')
    logger.info(master_trace_coefficients.data.shape)

    return master_trace_coefficients


def get_trace_coefficients(image, pipeline_context):
    coefficients_and_indices, fiber_order = None, None

    master_trace_filename = TraceMaker(pipeline_context).get_calibration_filename(image)
    master_trace_file_path = os.path.join(pipeline_context.processed_path, master_trace_filename)
    if image.header['OBSTYPE'] != 'TRACE' and os.path.isfile(master_trace_file_path):
        fiber_order = fits.getheader(master_trace_file_path).get('FIBRORDR')
        coefficients_and_indices = fits.getdata(master_trace_file_path)

        assert coefficients_and_indices is not None
        assert fiber_order is not None

    if image.header['OBSTYPE'] != 'LAMPFLAT' and not os.path.isfile(master_trace_file_path):
        raise MasterCalibrationDoesNotExist

    return coefficients_and_indices, fiber_order
