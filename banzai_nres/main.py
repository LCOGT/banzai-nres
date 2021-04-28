"""
main.py: Main driver script for the banzai-NRES pipeline.
    reduce_night() is the console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)

    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import banzai_nres.settings
from banzai.main import start_listener, parse_args, add_settings_to_context
from celery.schedules import crontab
from banzai.celery import app, schedule_calibration_stacking
import celery
import celery.bin.beat
from banzai.utils import date_utils, import_utils
from banzai import calibrations, logs
from banzai_nres import dbs
import banzai.dbs
import os
from banzai.data import DataProduct

import logging
import argparse
import requests

logger = logging.getLogger(__name__)


def nres_run_realtime_pipeline():
    extra_console_arguments = [{'args': ['--n-processes'],
                                'kwargs': {'dest': 'n_processes', 'default': 12,
                                           'help': 'Number of listener processes to spawn.', 'type': int}},
                               {'args': ['--queue-name'],
                                'kwargs': {'dest': 'queue_name', 'default': 'banzai_nres_pipeline',
                                           'help': 'Name of the queue to listen to from the fits exchange.'}}]

    runtime_context = parse_args(banzai_nres.settings, extra_console_arguments=extra_console_arguments)

    start_listener(runtime_context)


def nres_start_stacking_scheduler():
    logger.info('Entered entrypoint to celery beat scheduling')
    runtime_context = parse_args(banzai_nres.settings)
    for site, entry in runtime_context.SCHEDULE_STACKING_CRON_ENTRIES.items():
        app.add_periodic_task(crontab(minute=entry['minute'], hour=entry['hour']),
                              schedule_calibration_stacking.s(site=site, runtime_context=vars(runtime_context)))

    beat = celery.bin.beat.beat(app=app)
    logger.info('Starting celery beat')
    beat.run()


def nres_make_master_calibrations():
    extra_console_arguments = [{'args': ['--site'],
                                'kwargs': {'dest': 'site', 'help': 'Site code (e.g. ogg)', 'required': True}},
                               {'args': ['--camera'],
                                'kwargs': {'dest': 'camera', 'help': 'Camera (e.g. kb95)', 'required': True}},
                               {'args': ['--frame-type'],
                                'kwargs': {'dest': 'frame_type', 'help': 'Type of frames to process',
                                           'choices': ['bias', 'dark', 'lampflat'], 'required': True}},
                               {'args': ['--min-date'],
                                'kwargs': {'dest': 'min_date', 'required': True, 'type': date_utils.validate_date,
                                           'help': 'Earliest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}},
                               {'args': ['--max-date'],
                                'kwargs': {'dest': 'max_date', 'required': True, 'type': date_utils.validate_date,
                                           'help': 'Latest observation time of the individual calibration frames. '
                                                   'Must be in the format "YYYY-MM-DDThh:mm:ss".'}}]

    runtime_context = parse_args(banzai_nres.settings, extra_console_arguments=extra_console_arguments)
    instrument = banzai.dbs.query_for_instrument(runtime_context.db_address, runtime_context.site,
                                                 runtime_context.camera)
    calibrations.make_master_calibrations(instrument,  runtime_context.frame_type.upper(),
                                          runtime_context.min_date, runtime_context.max_date, runtime_context)


def add_bpm():
    parser = argparse.ArgumentParser(description="Add a bad pixel mask to the db.")
    parser.add_argument('--filename', help='Full path to Bad Pixel Mask file')
    parser.add_argument("--log-level", default='debug', choices=['debug', 'info', 'warning',
                                                                 'critical', 'fatal', 'error'])
    parser.add_argument('--db-address', dest='db_address',
                        default='mysql://cmccully:password@localhost/test',
                        help='Database address: Should be in SQLAlchemy form')
    args = parser.parse_args()
    add_settings_to_context(args, banzai_nres.settings)
    logs.set_log_level(args.log_level)
    frame_factory = import_utils.import_attribute(banzai_nres.settings.FRAME_FACTORY)()
    bpm_image = frame_factory.open({'path': args.filename}, args)
    bpm_image.is_master = True
    banzai.dbs.save_calibration_info(bpm_image.to_db_record(DataProduct(None, filename=os.path.basename(args.filename),
                                                                        filepath=os.path.dirname(args.filename))),
                                     args.db_address)


def add_bpms_from_archive():
    parser = argparse.ArgumentParser(description="Add bad pixel mask from a given archive api")
    parser.add_argument('--db-address', dest='db_address',
                        default='mysql://cmccully:password@localhost/test',
                        help='Database address: Should be in SQLAlchemy form')
    args = parser.parse_args()
    add_settings_to_context(args, banzai_nres.settings)
    # Query the archive for all bpm files
    url = f'{banzai_nres.settings.ARCHIVE_FRAME_URL}/?OBSTYPE=BPM'
    archive_auth_header = banzai_nres.settings.ARCHIVE_AUTH_HEADER
    response = requests.get(url, headers=archive_auth_header)
    response.raise_for_status()
    results = response.json()['results']

    # Load each one, saving the calibration info for each
    frame_factory = import_utils.import_attribute(banzai_nres.settings.FRAME_FACTORY)()
    for frame in results:
        frame['frameid'] = frame['id']
        bpm_image = frame_factory.open(frame, args)
        if bpm_image is not None:
            bpm_image.is_master = True
            banzai.dbs.save_calibration_info(bpm_image.to_db_record(DataProduct(None, filename=bpm_image.filename,
                                                                                filepath=None)), args.db_address)


def create_db():
    """
    Create the database structure.

    This only needs to be run once on initialization of the database.
    """
    parser = argparse.ArgumentParser("Create the database.\n\n"
                                     "This only needs to be run once on initialization of the database.")

    parser.add_argument("--log-level", default='debug', choices=['debug', 'info', 'warning',
                                                                 'critical', 'fatal', 'error'])
    parser.add_argument('--db-address', dest='db_address',
                        default='sqlite3:///test.db',
                        help='Database address: Should be in SQLAlchemy form')
    args = parser.parse_args()
    logs.set_log_level(args.log_level)

    dbs.create_db(args.db_address)


def populate_phoenix_models():
    parser = argparse.ArgumentParser("Populate the database with the Phoenix models.\n\n"
                                     "This only needs to be run once on initialization of the database.")
    parser.add_argument('--model-location', dest='model_location',
                        help='Location of the phoenix models. \
                        This should either be s3://bucket-name or an absolute directory path.')
    parser.add_argument("--log-level", default='debug', choices=['debug', 'info', 'warning',
                                                                 'critical', 'fatal', 'error'])
    parser.add_argument('--db-address', dest='db_address',
                        default='sqlite3:///test.db',
                        help='Database address: Should be in SQLAlchemy form')
    args = parser.parse_args()
    logs.set_log_level(args.log_level)

    dbs.populate_phoenix_models(args.model_location, args.db_address)
