"""
main.py: Main driver script for the banzai-NRES pipeline.
    reduce_night() is the console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)

    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
import banzai_nres.settings as nres_settings  # import to override banzai settings.
from banzai.main import run_realtime_pipeline

import logging

logger = logging.getLogger(__name__)


def nres_run_realtime_pipeline():
    run_realtime_pipeline()
