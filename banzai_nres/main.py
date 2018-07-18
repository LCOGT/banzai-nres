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

"""
TestContext and the functions in this following section are for
testing the partial pipeline without getting pipeline_context objects
from the actual array. e.g. from argparse.ArgumentParser etc.

elem_to_add = {lsc

"""
class TestContext(object):
    """
    Picks out a frame or a set of frames to test.
    Parameters
    ----------
    filename: None if you want just the path to be included in the context (only frames with OBSTYPE = 'BIAS' are used
    Returns
    -------
    stages_todo: list of banzai.stages.Stage
                 The stages that need to be done
    """
    def __init__(self,filename):
        _DEFAULT_DB = 'sqlite:////archive/engineering/test.db' #  from docker-compose file
        create_db('/archive/engineering/lsc/nres01/20180328/raw', db_address=_DEFAULT_DB,
                  configdb_address='http://configdb.lco.gtn/sites/')
        self.processed_path = '/tmp'
        self.raw_path = '/archive/engineering/lsc/nres01/20180328/raw'
        self.filename = filename
        self.post_to_archive = False
        self.db_address = _DEFAULT_DB

def test_making_master_biases():
    test_image_context = TestContext(None)
    image_list = image_utils.make_image_list(test_image_context)
    image_list = image_utils.select_images(image_list, image_types=['BIAS'])

    # spoofing the instrument name to assign each image a non-nres instrument name and
    # adding a bpm of zeros. - Monkey Patch
    for filename in image_list:
        fits.setval(filename, 'INSTRUME', value='ef06', ext=1)
        imagedata = fits.getdata(filename)
        bpm_array = np.zeros_like(imagedata)
        hdu_list = fits.open(filename)
        hdu_bpm = fits.ImageHDU(data=bpm_array, name='BPM')
        hdu_list.append(hdu_bpm)
        hdu_list.writeto(filename, overwrite=True)
        hdu_list.close()

    # End of monkey patch.

    print(make_master_bias(test_image_context))
    return True

"""
Start of the main.py that will remain after we no longer use the above.
"""

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

"""
Functions which are my horrible work around to not use mock at the moment.
An example of Mock is here: test_dbs.py

"""

def image_utils_no_db(pipeline_context, images, master_calibration = False):
    """
    An exact copy of image_utils.save_images without the part where it
    saves to a database.
    """
    output_files = []
    for image in images:
        output_directory = file_utils.make_output_directory(pipeline_context, image)

        if not master_calibration:
            image.filename = image.filename.replace('00.fits',
                                                    '{:02d}.fits'.format(int(pipeline_context.rlevel)))

        image_filename = os.path.basename(image.filename)
        filepath = os.path.join(output_directory, image_filename)
        output_files.append(filepath)

        save_pipeline_metadata(image, pipeline_context)
        image.writeto(filepath, pipeline_context.fpack)
        if pipeline_context.fpack:
            image_filename += '.fz'
            filepath += '.fz'

        if pipeline_context.post_to_archive:

            logger.info('Posting {filename} to the archive'.format(filename=image_filename))
            try:
                file_utils.post_to_archive_queue(filepath)
            except Exception as e:
                logger.error("Could not post {0} to ingester.".format(filepath))
                logger.error(e)
                continue
    return output_files


def read_images_zeros_BPM(image_list, pipeline_context):
    """
    This is an exact copy of banzai.images.read_images(image_list, pipeline_context)
    with a BPM of zeros.
    """
    images = []
    for filename in image_list:
        try:

            image = Image(pipeline_context, filename=filename)

            munge(image, pipeline_context)
            images.append(image)
        except Exception as e:
            logger.error('Error loading {0}'.format(filename))
            logger.error(e)
            continue
    return images

"""
end of horrible ad-hoc functions
"""

def run(stages_to_do, pipeline_context, image_types=[], calibration_maker=False, log_message=''):
    """
    Main driver script for banzai-NRES
    """
    if len(log_message) > 0:
        logger.info(log_message, extra={'tags': {'raw_path': pipeline_context.raw_path}})

    image_list = image_utils.make_image_list(pipeline_context)
    """
    image_list does the following: given the pipeline_context object (file path info etc) we construct the list of images we will analyze.
    if pipeline_context.filename == None, then we iterate through all the files in the directory.
    """
    image_list = image_utils.select_images(image_list, image_types)

    images = banzai.images.read_images(image_list, pipeline_context) # this makes a call to db_address only if site or instrument are both not None

    for stage in stages_to_do:
        stage_to_run = stage(pipeline_context)  # isolate the stage that will be run
        images = stage_to_run.run(images)   # update the list of images after running the stage on them.

    output_files = image_utils.save_images(pipeline_context, images,
                                           master_calibration=calibration_maker)
    # End of Monkey Patch
    return output_files
