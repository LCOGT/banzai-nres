"""
traces.py: Driver scripts for finding echelle orders across a CCD.

Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai_nres.utils.trace_utils import Trace, AllTraceFitter
from banzai.stages import Stage
from banzai.calibrations import CalibrationMaker, create_master_calibration_header
from banzai.utils import file_utils
from banzai import dbs

import sep
import os
import logging

logger = logging.getLogger(__name__)


class TraceMaker(CalibrationMaker):
    def __init__(self, pipeline_context):
        super(TraceMaker, self).__init__(pipeline_context)
        self.pipeline_context = pipeline_context
        self.order_of_poly_fit = self.pipeline_context.TRACE_FIT_POLYNOMIAL_ORDER
        self.second_order_coefficient_guess = self.pipeline_context.TRACE_FIT_INITIAL_DEGREE_TWO_GUESS
        self.trace_table_name = self.pipeline_context.TRACE_TABLE_NAME
        self.xmin = self.pipeline_context.WINDOW_FOR_TRACE_IDENTIFICATION['min']
        self.xmax = self.pipeline_context.WINDOW_FOR_TRACE_IDENTIFICATION['max']
        self.min_peak_to_peak_spacing = self.pipeline_context.MIN_PEAK_TO_PEAK_SPACING
        self.min_snr = self.pipeline_context.MIN_SNR_FOR_TRACE_IDENTIFICATION

    @property
    def calibration_type(self):
        return 'TRACE'

    def make_master_calibration_frame(self, images):
        traces = []
        for image in images:
            master_header = create_master_calibration_header(old_header=image.header, images=[image])
            master_filename = self.pipeline_context.CALIBRATION_FILENAME_FUNCTIONS[self.calibration_type](image)
            master_filepath = self._get_filepath(self.pipeline_context, image, master_filename)
            bkg_subtracted_image_data = image.data - sep.Background(image.data).back()

            logger.info('fitting traces order by order', image=image)
            fitter = AllTraceFitter(xmin=self.xmin, xmax=self.xmax,
                                    min_peak_to_peak_spacing=self.min_peak_to_peak_spacing,
                                    min_snr=self.min_snr)
            trace = Trace(data=None, filepath=master_filepath, header=master_header,
                          num_centers_per_trace=image.data.shape[1], trace_table_name=self.trace_table_name)
            trace = fitter.fit_traces(trace=trace, image_data=bkg_subtracted_image_data,
                                      poly_fit_order=self.order_of_poly_fit,
                                      second_order_coefficient_guess=self.second_order_coefficient_guess,
                                      image_noise_estimate=image.header['RDNOISE'])
            traces.append(trace)
            logger.info('Created master trace', image=image, extra_tags={'calibration_type': self.calibration_type,
                                                                         'output_path': master_filepath})
        return traces

    @staticmethod
    def _get_filepath(pipeline_context, lampflat_image, master_filename):
        output_directory = file_utils.make_output_directory(pipeline_context, lampflat_image)
        return os.path.join(output_directory, os.path.basename(master_filename))

    def do_stage(self, images):
        """
        :param images: list of images to run TraceMaker on.
        :return: [first_trace, second_trace,...] etc. This is a munge function which fixes the following problem.
        TraceMaker.make_master_calibration_frame returns [first_trace, second_trace,...] and then BANZAI's do_stage
        wraps that list in another list, e.g. do_stage natively would return [[first_trace, second_trace,...]]. This
        rewrite of do_stage simply returns [[first_trace, second_trace,...]][0] = [first_trace, second_trace,...]
        """
        master_calibrations = super(TraceMaker, self).do_stage(images)
        if len(master_calibrations) > 0:
            master_calibrations = master_calibrations[0]
        return master_calibrations


class LoadTrace(Stage):
    """
    Loads trace coefficients from file and appends them onto the image object.
    """
    #TODO make this a calibration applier stage.
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

            if master_trace_path is None or not os.path.exists(master_trace_path):
                logger.error('Master trace fit file not found for {0}.'.format(image.filename))
                logger.error("Can't find trace coefficients for image, stopping reduction", image=image)
                images_to_remove.append(image)
                continue
            else:
                image.trace = Trace.load(master_trace_path, trace_table_name=self.pipeline_context.TRACE_TABLE_NAME)
        for image in images_to_remove:
            images.remove(image)
        return images
