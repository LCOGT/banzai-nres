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
from banzai.main import process_directory, parse_directory_args, parse_args
from banzai import dbs, logs
from banzai.utils import date_utils
import os

import logging

logger = logging.getLogger(__name__)


def make_master_bias(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRESSettings())
    process_directory(pipeline_context, raw_path, ['BIAS'], last_stage=pipeline_context.LAST_STAGE['BIAS'],
                      extra_stages=pipeline_context.EXTRA_STAGES['BIAS'],
                      log_message='Making Master BIAS', calibration_maker=True)


def make_master_dark(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRESSettings())
    process_directory(pipeline_context, raw_path, ['DARK'], last_stage=pipeline_context.LAST_STAGE['DARK'],
                      extra_stages=pipeline_context.EXTRA_STAGES['DARK'],
                      log_message='Making Master Dark', calibration_maker=True)


def make_master_flat(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRESSettings())
    process_directory(pipeline_context, raw_path, ['LAMPFLAT'], last_stage=pipeline_context.LAST_STAGE['LAMPFLAT'],
                      extra_stages=pipeline_context.EXTRA_STAGES['LAMPFLAT'],
                      log_message='Making Master Flat', calibration_maker=True)


def reduce_night():
    extra_console_arguments = [{'args': ['--site'], 'kwargs': {'dest': 'site', 'help': 'Site code (e.g. ogg)'}},
                               {'args': ['--dayobs'], 'kwargs': {'dest': 'dayobs', 'default': None,
                                                               'help': 'Day-Obs to reduce (e.g. 20160201)'}},
                               {'args': ['--raw-path-root'],
                                'kwargs': {'dest': 'rawpath_root', 'default': '/archive/engineering',
                                           'help': 'Top level directory with raw data.'}}]

    pipeline_context = parse_args(settings.NRESSettings(), extra_console_arguments=extra_console_arguments,
                                  parser_description='Reduce all the data from a site at the end of a night.')

    # Ping the configdb to get instruments
    try:
        dbs.populate_instrument_tables(db_address=pipeline_context.db_address)
    except Exception:
        logger.error('Could not connect to the configdb: {error}'.format(error=logs.format_exception()))

    try:
        timezone = dbs.get_timezone(pipeline_context.site, db_address=pipeline_context.db_address)
    except dbs.SiteMissingException:
        logger.error("Site {site} not found in database {db}, exiting.".format(site=pipeline_context.site,
                                                                               db=pipeline_context.db_address),
                     extra_tags={'site': pipeline_context.site})
        return

    instruments = dbs.get_instruments_at_site(pipeline_context.site,
                                              db_address=pipeline_context.db_address,
                                              ignore_schedulability=pipeline_context.ignore_schedulability)

    # If no dayobs is given, calculate it.
    if pipeline_context.dayobs is None:
        dayobs = date_utils.get_dayobs(timezone=timezone)
    else:
        dayobs = pipeline_context.dayobs

    # For each instrument at the given site
    for instrument in instruments:
        raw_path = os.path.join(pipeline_context.rawpath_root, pipeline_context.site,
                                instrument.camera, dayobs, 'raw')

        # Run the reductions on the given dayobs
        try:
            make_master_bias(pipeline_context=pipeline_context, raw_path=raw_path)
        except Exception:
            logger.error(logs.format_exception())
        try:
            make_master_dark(pipeline_context=pipeline_context, raw_path=raw_path)
        except Exception:
            logger.error(logs.format_exception())
        try:
            make_master_flat(pipeline_context=pipeline_context, raw_path=raw_path)
        except Exception:
            logger.error(logs.format_exception())
