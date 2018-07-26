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

from banzai import bias, trim, dark, gain
from banzai.qc import header_checker
from banzai import logs
from banzai.utils import image_utils
from banzai.main import make_master_dark, make_master_bias

logger = logs.get_logger(__name__)

ordered_stages = [header_checker.HeaderSanity,
                  bias.OverscanSubtractor,
                  gain.GainNormalizer,
                  trim.Trimmer,
                  bias.BiasSubtractor,
                  dark.DarkSubtractor]


class TestContext(object):
    """
    Picks out a frame or a set of frames to test. Provides the appropriate
    PipelineContext to pass an NRES test frame to the banzai modules.
    Parameters
    ----------
    filename: None if you want just the path to be included in the context (only frames with OBSTYPE = 'BIAS' are used
    Returns
    -------
    stages_todo: list of banzai.stages.Stage
                 The stages that need to be done
    """
    def __init__(self, filename=None, raw_path='/archive/engineering/lsc/nres01/20180228/raw'):
        self.processed_path = '/tmp'
        self.raw_path = raw_path
        self.filename = filename
        self.post_to_archive = False
        self.db_address = os.environ['DB_URL']
        self.preview_mode = False
        self.rlevel = 0
        self.fpack = True


def parse_end_of_night_command_line_arguments():
    """
    :return: Directory where test NRES frames live. Eventually this would be hooked up to the
    pipeline, instead of giving a fixed directory.
    """
    return TestContext()


def run_end_of_night_from_console(scripts_to_run):
    pipeline_context = parse_end_of_night_command_line_arguments()
    for script in scripts_to_run:
        script(pipeline_context)


def make_master_bias_console():
    run_end_of_night_from_console([make_master_bias])


def make_master_dark_console():
    run_end_of_night_from_console([make_master_dark])


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
        images = stage_to_run.run(images)


    output_files = image_utils.save_images(pipeline_context, images,
                                           master_calibration=calibration_maker)

    logger.info(str(output_files[0]), extra={'tags': {'raw_path': pipeline_context.raw_path}})
