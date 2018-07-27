"""
main.py: Main driver script for the banzai-NRES pipeline.
    The main() function is a console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)
July 2018
    G. Mirek Brandt (gmbrandt@ucsb.edu)
July 2018
"""

import os
from banzai_nres.utils.image_utils import read_images
from banzai_nres.bias import BiasMaker as nres_BiasMaker
from banzai_nres.tests.utils import TestContext

from banzai import bias, trim, dark, gain
from banzai import logs
from banzai.utils import image_utils
from banzai.main import get_stages_todo
from banzai import main as banzai_main

logger = logs.get_logger(__name__)


banzai_main.ordered_stages = [bias.OverscanSubtractor,
                              gain.GainNormalizer,
                              trim.Trimmer,
                              bias.BiasSubtractor]


def parse_end_of_night_command_line_arguments():
    """
    Initiates the TestContext which points to NRES frames in the test directory.
    """
    return TestContext()


def run_end_of_night_from_console(scripts_to_run):
    pipeline_context = parse_end_of_night_command_line_arguments()
    for script in scripts_to_run:
        script(pipeline_context)


def make_master_bias_console():
    run_end_of_night_from_console([make_master_bias])


def make_master_bias(pipeline_context):
    stages_to_do = get_stages_todo(trim.Trimmer, extra_stages=[nres_BiasMaker])
    run(stages_to_do, pipeline_context, image_types=['BIAS'], calibration_maker=True,
        log_message='Making Master BIAS')


def run(stages_to_do, pipeline_context, image_types=[], calibration_maker=False, log_message=''):
    """
    Main driver script for banzai-NRES
    """
    if len(log_message) > 0:
        logger.info(log_message, extra={'tags': {'raw_path': pipeline_context.raw_path}})

    image_list = image_utils.make_image_list(pipeline_context)

    image_list = image_utils.select_images(image_list, image_types)

    images = read_images(image_list, pipeline_context)

    for stage in stages_to_do:
        stage_to_run = stage(pipeline_context)
        logger.info('running' + str(stage_to_run))
        images = stage_to_run.run(images)

    output_files = image_utils.save_images(pipeline_context, images,
                                           master_calibration=calibration_maker)

    logger.info(output_files, extra={'tags': {'raw_path': pipeline_context.raw_path}})
