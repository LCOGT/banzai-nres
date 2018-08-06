"""
main.py: Main driver script for the banzai-NRES pipeline.
    The make_master_bias_console() and make_master_dark_console() etc.. are the entry points.
Authors
    Curtis McCully (cmccully@lcogt.net)
July 2018
    G. Mirek Brandt (gmbrandt@ucsb.edu)
July 2018
"""

from banzai_nres.bias import BiasMaker as nres_BiasMaker
from banzai_nres.dark import DarkMaker as nres_DarkMaker
from banzai_nres import traces
import argparse

from banzai import bias, trim, dark, gain
from banzai import logs, qc
from banzai.main import get_stages_todo, run_end_of_night_from_console, run
from banzai import main as banzai_main

logger = logs.get_logger(__name__)


banzai_main.ordered_stages = [qc.HeaderSanity,
                              qc.ThousandsTest,
                              qc.SaturationTest,
                              bias.OverscanSubtractor,
                              gain.GainNormalizer,
                              trim.Trimmer,
                              bias.BiasSubtractor,
                              dark.DarkSubtractor]


def make_master_bias_console():
    run_end_of_night_from_console([make_master_bias])


def make_master_dark_console():
    run_end_of_night_from_console([make_master_dark])


def make_master_trace_console():
    run_end_of_night_from_console([make_master_trace])


def make_master_trace_blind_console():
    run_end_of_night_from_console([make_master_trace_blind])


def make_master_bias(pipeline_context):
    stages_to_do = get_stages_todo(trim.Trimmer, extra_stages=[nres_BiasMaker])
    run(stages_to_do, pipeline_context, image_types=['BIAS'], calibration_maker=True,
        log_message='Making Master BIAS')


def make_master_dark(pipeline_context):
    stages_to_do = get_stages_todo(bias.BiasSubtractor, extra_stages=[dark.DarkNormalizer, nres_DarkMaker])
    run(stages_to_do, pipeline_context, image_types=['DARK'], calibration_maker=True,
        log_message='Making Master Dark')


def make_master_trace(pipeline_context):
    stages_to_do = get_stages_todo(dark.DarkSubtractor, extra_stages=[traces.TraceFitOrderbyOrder])
    run(stages_to_do, pipeline_context, image_types=['LAMPFLAT'], calibration_maker=True,
        log_message='Making Master Trace by Updating Previous Master with global-meta Technique')


def make_master_trace_blind(pipeline_context):
    stages_to_do = get_stages_todo(dark.DarkSubtractor, extra_stages=[traces.TraceFitOrderbyOrder])
    run(stages_to_do, pipeline_context, image_types=['LAMPFLAT'], calibration_maker=True,
        log_message='Making Master Trace via order-by-order Blind Technique!')