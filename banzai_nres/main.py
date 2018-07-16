"""
main.py: Main driver script for the banzai-NRES pipeline.
    The main() function is a console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)
July 2018
    G. Mirek Brandt (gmbrandt@ucsb.edu)
July 2018
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import multiprocessing
import os
import traceback
import sys

from kombu import Connection, Queue, Exchange
from kombu.mixins import ConsumerMixin

import banzai.images
from banzai import bias, trim
from banzai import logs
from banzai.utils import image_utils

import banzai.tests.utils

logger = logs.get_logger(__name__)

ordered_stages = [bias.OverscanSubtractor,
                  trim.Trimmer,
                  bias.BiasSubtractor,
                  ]

class TestContext(object):
    """
    Picks out a frame to test.
    """
    def __init__(self,filename):
        self.processed_path = '/archive/engineering/lsc/nres01/20180328/tmp'
        self.raw_path = '/archive/engineering/lsc/nres01/20180328/raw'
        self.filename = filename

def test_one_image():
    test_image_context = TestContext('lscnrs01-fl09-20180328-0001-w00.fits.fz
')
    print(make_master_bias(test_image_context))
    return True


def get_stages_todo(last_stage=None, extra_stages=None):
    """
    Parameters
    ----------
    last_stage: banzai.stages.Stage
                Last stage to do
    extra_stages: Stages to do after the last stage
    Returns
    -------
    stages_todo: list of banzai.stages.Stage
                 The stages that need to be done
    Notes
    -----
    Extra stages can be other stages that are not in the ordered_stages list.
    """
    if extra_stages is None:
        extra_stages = []

    if last_stage is None:
        last_index = None
    else:
        last_index = ordered_stages.index(last_stage) + 1

    stages_todo = ordered_stages[:last_index] + extra_stages
    return stages_todo


class PipelineContext(object):
    """
    contains some basic information about the images to be processed.
    """
    def __init__(self, args):
        self.processed_path = args.processed_path
        self.raw_path = args.raw_path
        self.filename = args.filename



def make_master_bias(pipeline_context):
    """
    usually input is a PipelineContext object, however we define it herein for testing.
    :return:
    master bias and saves the images.
    """
    stages_to_do = get_stages_todo(trim.Trimmer, extra_stages=[bias.BiasMaker])
    run(stages_to_do, pipeline_context, image_types=['BIAS'], calibration_maker=True,
        log_message='Making Master BIAS')


def reduce_science_frames(pipeline_context):
    stages_to_do = get_stages_todo()
    reduce_frames_one_by_one(stages_to_do, pipeline_context)


def reduce_experimental_frames(pipeline_context):
    stages_to_do = get_stages_todo()
    reduce_frames_one_by_one(stages_to_do, pipeline_context, image_types=['EXPERIMENTAL'])


def reduce_frames_one_by_one(stages_to_do, pipeline_context, image_types=None):
    if image_types is None:
        image_types = ['EXPOSE', 'STANDARD']
    image_list = image_utils.make_image_list(pipeline_context)
    original_filename = pipeline_context.filename
    for image in image_list:
        pipeline_context.filename = os.path.basename(image)
        try:
            run(stages_to_do, pipeline_context, image_types=image_types)
        except Exception as e:
            logger.error('{0}'.format(e), extra={'tags': {'filename': pipeline_context.filename,
                                                          'filepath': pipeline_context.raw_path}})
    pipeline_context.filename = original_filename


def run(stages_to_do, pipeline_context, image_types=[], calibration_maker=False, log_message=''):
    """
    Main driver script for banzai-NRES
    """
    if len(log_message) > 0:
        logger.info(log_message, extra={'tags': {'raw_path': pipeline_context.raw_path}})

    image_list = image_utils.make_image_list(pipeline_context)
    image_list = image_utils.select_images(image_list, image_types)
    images = banzai.images.read_images(image_list, pipeline_context)

    for stage in stages_to_do:
        stage_to_run = stage(pipeline_context)
        images = stage_to_run.run(images)

    output_files = image_utils.save_images(pipeline_context, images,
                                           master_calibration=calibration_maker)
    return output_files
