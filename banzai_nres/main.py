"""
main.py: Main driver script for the banzai-NRES pipeline.
    The make_master_bias_console() and make_master_dark_console() are the entry points.
Authors
    Curtis McCully (cmccully@lcogt.net)
July 2018
    G. Mirek Brandt (gmbrandt@ucsb.edu)
July 2018
"""

from banzai_nres.bias import BiasMaker as nres_BiasMaker
from banzai_nres.dark import DarkMaker as nres_DarkMaker

from banzai import bias, trim, dark, gain
from banzai import qc
from banzai.main import get_stages_todo, process_directory, parse_directory_args
from banzai import main as banzai_main
from banzai.context import TelescopeCriterion
import operator
import logging

logger = logging.getLogger(__name__)

NRES_CRITERIA = [TelescopeCriterion('camera_type', operator.contains, 'NRES'),
                 TelescopeCriterion('schedulable', operator.eq, True)]

banzai_main.ORDERED_STAGES = [qc.HeaderSanity,
                              qc.ThousandsTest,
                              qc.SaturationTest,
                              bias.OverscanSubtractor,
                              gain.GainNormalizer,
                              trim.Trimmer,
                              bias.BiasSubtractor,
                              dark.DarkSubtractor]


def make_master_bias(pipeline_context=None, raw_path=None):
    pipeline_context, raw_path = parse_directory_args(pipeline_context, raw_path, NRES_CRITERIA)
    process_directory(pipeline_context, raw_path, ['BIAS'], last_stage=trim.Trimmer,
                      extra_stages=[bias.BiasMasterLevelSubtractor, nres_BiasMaker],
                      log_message='Making Master BIAS', calibration_maker=True)


def make_master_dark(pipeline_context=None, raw_path=None):
    pipeline_context = parse_directory_args(pipeline_context, raw_path, NRES_CRITERIA,
                                            parser_description='Reduce LCO NRES data.')
    process_directory(pipeline_context, raw_path, ['DARK'], last_stage=bias.BiasSubtractor,
                      extra_stages=[dark.DarkNormalizer, nres_DarkMaker],
                      log_message='Making Master Dark', calibration_maker=True)
