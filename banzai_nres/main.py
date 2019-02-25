"""
main.py: Main driver script for the banzai-NRES pipeline.
    The make_master_bias_console() and make_master_dark_console() are the entry points.
Authors
    Curtis McCully (cmccully@lcogt.net)
July 2018
    G. Mirek Brandt (gmbrandt@ucsb.edu)
July 2018
"""

from banzai_nres import settings
from banzai.main import process_directory, parse_directory_args
from banzai.utils import date_utils
from banzai import dbs
from banzai import logs
from banzai.main import run_master_maker

import logging

logger = logging.getLogger(__name__)


def reduce_bias_frames(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRESSettings())
    process_directory(pipeline_context, raw_path, ['BIAS'])


def reduce_dark_frames(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRESSettings())
    process_directory(pipeline_context, raw_path, ['DARK'])


def reduce_flat_frames(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRESSettings())
    process_directory(pipeline_context, raw_path, ['LAMPFLAT'])


def process_master_maker(pipeline_context, instrument, frame_type_to_stack, min_date, max_date,
                         master_frame_type=None, use_masters=False):
    if master_frame_type is None:
        master_frame_type = frame_type_to_stack
    extra_tags = {'instrument': instrument.camera, 'master_frame_type': master_frame_type,
                  'min_date': min_date.strftime(date_utils.TIMESTAMP_FORMAT),
                  'max_date': max_date.strftime(date_utils.TIMESTAMP_FORMAT)}
    logger.info("Making master frames", extra_tags=extra_tags)
    image_path_list = dbs.get_individual_calibration_images(instrument, frame_type_to_stack, min_date, max_date,
                                                            use_masters=use_masters,
                                                            db_address=pipeline_context.db_address)
    if len(image_path_list) == 0:
        logger.info("No calibration frames found to stack", extra_tags=extra_tags)

    try:
        run_master_maker(image_path_list, pipeline_context, master_frame_type)
    except Exception:
        logger.error(logs.format_exception())


def stack_calibrations(pipeline_context=None):
    nres_settings = settings.NRESSettings()
    extra_console_arguments = [{'args': ['--site'],
                                'kwargs': {'dest': 'site', 'help': 'Site code (e.g. ogg)', 'required': True}},
                               {'args': ['--camera'],
                                'kwargs': {'dest': 'camera', 'help': 'Camera (e.g. kb95)', 'required': True}},
                               {'args': ['--frame-type'],
                                'kwargs': {'dest': 'frame_type', 'help': 'Type of frames to process',
                                           'choices': nres_settings.CALIBRATION_STACKER_STAGE.keys(), 'required': True}},
                               {'args': ['--min-date'],
                                'kwargs': {'dest': 'min_date', 'required': True, 'type': date_utils.valid_date,
                                           'help': 'Earliest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}},
                               {'args': ['--max-date'],
                                'kwargs': {'dest': 'max_date', 'required': True, 'type': date_utils.valid_date,
                                           'help': 'Latest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}}]

    pipeline_context, raw_path = parse_directory_args(pipeline_context, None, nres_settings,
                                                      extra_console_arguments=extra_console_arguments)
    instrument = dbs.query_for_instrument(pipeline_context.db_address, pipeline_context.site, pipeline_context.camera)
    if pipeline_context.frame_type == 'TRACE':
        frame_type_to_stack = 'LAMPFLAT'
        use_masters = True
        master_frame_type = 'TRACE'
    else:
        frame_type_to_stack = pipeline_context.frame_type
        use_masters = False
        master_frame_type = None
    process_master_maker(pipeline_context, instrument,  frame_type_to_stack.upper(),
                         pipeline_context.min_date, pipeline_context.max_date,
                         master_frame_type=master_frame_type, use_masters=use_masters)
