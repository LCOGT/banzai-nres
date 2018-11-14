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
import logging

logger = logging.getLogger(__name__)


def make_master_bias(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRES_CRITERIA,
                                                      settings.IMAGE_CLASS)
    process_directory(pipeline_context, raw_path, ['BIAS'], last_stage=settings.BIAS_LAST_STAGE,
                      extra_stages=settings.BIAS_EXTRA_STAGES,
                      log_message='Making Master BIAS', calibration_maker=True)


def make_master_dark(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRES_CRITERIA,
                                                      settings.IMAGE_CLASS)
    process_directory(pipeline_context, raw_path, ['DARK'], last_stage=settings.DARK_LAST_STAGE,
                      extra_stages=settings.DARK_EXTRA_STAGES,
                      log_message='Making Master Dark', calibration_maker=True)


def make_master_flat(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, settings.NRES_CRITERIA,
                                                      settings.IMAGE_CLASS)
    process_directory(pipeline_context, raw_path, ['LAMPFLAT'], last_stage=settings.FLAT_LAST_STAGE,
                      extra_stages=settings.FLAT_EXTRA_STAGES,
                      log_message='Making Master Dark', calibration_maker=True)
