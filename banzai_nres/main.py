"""
main.py: Main driver script for the banzai-NRES pipeline.
    reduce_night() is the console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)

    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import datetime

from banzai_nres.settings import NRESSettings
from banzai.main import process_directory, parse_directory_args
from banzai.utils import date_utils
from banzai import dbs
from banzai import logs
from banzai.main import run_master_maker
import os

import logging

logger = logging.getLogger(__name__)


class ReductionCriterion(object):
    def __init__(self, raw_path=None, pipeline_context=None, settings=NRESSettings()):
        max_date = getattr(pipeline_context, 'max_date', None)
        min_date = getattr(pipeline_context, 'min_date', None)
        if max_date is None:
            max_date = datetime.datetime.utcnow()

        if min_date is None:
            min_date = max_date - datetime.timedelta(hours=24)

        if min_date > max_date:
            logger.error('The start cannot be after the end. Aborting reduction!')
            raise ValueError('min_date > max_date.')
        self.min_date, self.max_date = min_date, max_date

        if getattr(pipeline_context, 'frame_type', None) is None:
            self.frame_types = settings.REDUCE_NIGHT_FRAME_TYPES
        else:
            self.frame_types = [pipeline_context.frame_type]
        self.raw_path = raw_path
        if raw_path is not None:
            timezone = dbs.get_timezone(pipeline_context.site, db_address=pipeline_context.db_address)
            dayobs = date_utils.get_dayobs(timezone=timezone)
            self.raw_path = os.path.join(raw_path, pipeline_context.site, pipeline_context.instrument_name, dayobs, 'raw')


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


def reduce_night(pipeline_context=None):
    nres_settings = NRESSettings()
    extra_console_arguments = [{'args': ['--site'],
                                'kwargs': {'dest': 'site', 'help': 'Site code (e.g. ogg)', 'required': True}},
                               {'args': ['--camera'],
                                'kwargs': {'dest': 'camera', 'help': 'Camera (e.g. kb95)', 'required': True}},
                               {'args': ['--instrument-name'],
                                'kwargs': {'dest': 'instrument_name', 'help': 'Instrument (e.g. nres04)', 'required': True}},
                               {'args': ['--enclosure'],
                                'kwargs': {'dest': 'enclosure', 'help': 'Enclosure code (e.g. clma)',
                                           'required': False}},
                               {'args': ['--telescope'],
                                'kwargs': {'dest': 'telescope', 'help': 'Telescope code (e.g. 0m4a)',
                                           'required': False}},
                               {'args': ['--frame-type'],
                                'kwargs': {'dest': 'frame_type', 'help': 'Type of frames to process',
                                           'choices': nres_settings.CALIBRATION_STACKER_STAGE.keys(), 'required': False}},
                               {'args': ['--min-date'],
                                'kwargs': {'dest': 'min_date', 'required': False, 'type': date_utils.valid_date,
                                           'help': 'Earliest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}},
                               {'args': ['--max-date'],
                                'kwargs': {'dest': 'max_date', 'required': False, 'type': date_utils.valid_date,
                                           'help': 'Latest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}}]

    pipeline_context, raw_path = parse_directory_args(pipeline_context, None, settings=nres_settings,
                                                      extra_console_arguments=extra_console_arguments)
    instrument = dbs.query_for_instrument(pipeline_context.db_address, pipeline_context.site,
                                          camera=pipeline_context.camera, name=pipeline_context.instrument_name,
                                          enclosure=None, telescope=None)
    reduction_criterion = ReductionCriterion(raw_path, pipeline_context, settings=nres_settings)

    for frame_type in reduction_criterion.frame_types:
        if frame_type == 'TRACE':
            frame_type_to_stack = 'LAMPFLAT'
            use_masters = True
            master_frame_type = 'TRACE'
        else:
            frame_type_to_stack = frame_type
            use_masters = False
            master_frame_type = None
            # must reduce frames before making the master calibration, unless we are making a master trace.
            process_directory(pipeline_context, reduction_criterion.raw_path, [frame_type_to_stack])

        process_master_maker(pipeline_context, instrument,  frame_type_to_stack.upper(),
                             min_date=reduction_criterion.min_date, max_date=reduction_criterion.max_date,
                             master_frame_type=master_frame_type, use_masters=use_masters)
