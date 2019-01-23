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
from banzai.main import process_directory, parse_directory_args, run_end_of_night
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
    run_end_of_night(settings.NRESSettings(), [make_master_bias, make_master_dark, make_master_flat])
