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
import numpy as np
import mock
import traceback
import sys

from kombu import Connection, Queue, Exchange
from kombu.mixins import ConsumerMixin

from astropy.io import fits

import banzai.images
from banzai import bias, trim
from banzai import logs
from banzai.utils import image_utils
from banzai.utils.image_utils import save_pipeline_metadata
from banzai.images import Image

from banzai.munge import munge

from banzai.utils import file_utils
import banzai.tests.utils

from banzai.dbs import create_db, add_or_update_record

logger = logs.get_logger(__name__)

ordered_stages = [bias.OverscanSubtractor,
                  trim.Trimmer,
                  bias.BiasSubtractor,
                  ]

# begin functions which act as entry points

# end of entry point functions

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
        self.post_to_archive = args.post_to_archive
        self.post_to_elasticsearch = args.post_to_elasticsearch
        self.elasticsearch_url = args.elasticsearch_url
        self.fpack = args.fpack
        self.rlevel = args.rlevel
        self.db_address = args.db_address
        self.log_level = args.log_level
        self.preview_mode = args.preview_mode
        self.filename = args.filename
        self.max_preview_tries = args.max_preview_tries
        self.elasticsearch_doc_type = args.elasticsearch_doc_type
        self.elasticsearch_qc_index = args.elasticsearch_qc_index




def make_master_bias(pipeline_context):
    """
    Returns:
    master bias and saves the images.
    Note:
    image_types = ['BIAS'] selects only images which are bias type, naturally.
    """
    stages_to_do = get_stages_todo(trim.Trimmer, extra_stages=[bias.BiasMaker])
    output_files = run(stages_to_do, pipeline_context, image_types=['BIAS'], calibration_maker=True,
        log_message='Making Master BIAS')
    return output_files


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


def read_images_fixed(image_list, pipeline_context):
    """
    This is a copy of banzai.images.read_images
    which will properly handle images which already have a Bad Pixel Mask (BPM)
    as an extension in the fits file. Prior, if image.bpm existed, the main.run
    program would output an empty list. All that has been added is

            else:
                images.append(image)

    Parameters:
        pipeline_context: Object which contains attributes which describe the database etc.
        image_list: A list of path/filename to fits files.
    Returns:
        images: List of banzai.images.Image objects with attached bad pixel masks.
    """

    images = []
    for filename in image_list:
        try:
            image = Image(pipeline_context, filename=filename)
            munge(image, pipeline_context)
            if image.bpm is None:
                bpm = image_utils.get_bpm(image, pipeline_context)
                if bpm is None:
                    logger.error('No BPM file exists for this image.',
                                 extra={'tags': {'filename': image.filename}})
                else:
                    image.bpm = bpm
                    images.append(image)
            else:
                images.append(image)
        except Exception as e:
            logger.error('Error loading {0}'.format(filename))
            logger.error(e)
            continue
    return images



def run(stages_to_do, pipeline_context, image_types=[], calibration_maker=False, log_message=''):
    """
    Main driver script for banzai-NRES

    Note to self:     image_list does the following: given the pipeline_context object (file path info etc) we construct the list of images we will analyze.
                    based off of image_types.
    if pipeline_context.filename == None, then we iterate through all the files in the directory.
    """
    if len(log_message) > 0:
        logger.info(log_message, extra={'tags': {'raw_path': pipeline_context.raw_path}})

    image_list = image_utils.make_image_list(pipeline_context)

    image_list = image_utils.select_images(image_list, image_types)

    images = read_images_fixed(image_list, pipeline_context) #  in banzai.main this is banzai.images.read_images - but that function does nothing if image.bpm is not None

    for stage in stages_to_do:
        stage_to_run = stage(pipeline_context)  # isolate the stage that will be run
        images = stage_to_run.run(images)   # update the list of images after running the stage on them.


    output_files = image_utils.save_images(pipeline_context, images,
                                           master_calibration=calibration_maker)
    return output_files
