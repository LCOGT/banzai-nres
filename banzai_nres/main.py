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
testing the partial pipeline on NRES frames.

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
    def __init__(self, filename):
        _DEFAULT_DB = 'sqlite:////archive/engineering/test.db' #  from docker-compose file
        create_db('/archive/engineering/lsc/nres01/20180328/raw', db_address=_DEFAULT_DB,
                  configdb_address='http://configdb.lco.gtn/sites/')
        self.processed_path = '/tmp'
        self.raw_path = '/archive/engineering/lsc/nres01/20180313/raw'
        self.filename = filename
        self.post_to_archive = False
        self.db_address = _DEFAULT_DB
        self.preview_mode = False
        self.rlevel = 0
        self.fpack = True


def test_making_master_biases():
    test_image_context = TestContext(None)
    image_list = image_utils.make_image_list(test_image_context)
    image_list = image_utils.select_images(image_list, image_types=['BIAS'])

    # Patching missing info so that files will pass through banzai's functions.
    for filename in image_list:
        # spoofing instrument name for one which banzai accepts has a database
        # this is used in building the image as a banzai.images.Image object.
        fits.setval(filename, 'INSTRUME', value='ef06', ext=1)
        # loading the image and building the null bad-pixel-mask if it needs it.
        need_null_bpm = False
        with fits.open(filename) as hdu_list:
            try:
                hdu_list['BPM']

            except:
                print('BPM key does not exist for %s \n creating and appending a zeros bad pixel mask' %filename)
                need_null_bpm = True

        hdu_list = fits.open(filename)
        if need_null_bpm:
            imagedata = fits.getdata(filename)
            bpm_array = np.zeros_like(imagedata)
            hdu_bpm = fits.ImageHDU(data=bpm_array, name='BPM')
            # Appending a bad pixel mask to the image.
            hdu_list.append(hdu_bpm)
        # loading the primary HDU header
        p_hdu_header = hdu_list[0].header
        # headers used in _trim_image
        p_hdu_header.set('CRPIX1', 0)
        p_hdu_header.set('CRPIX2', 0)
        p_hdu_header.set('L1STATTR', 0)
        # headers used in save_pipeline_metadata
        p_hdu_header.set('RLEVEL', 'TBD')  # reduction level
        p_hdu_header.set('PIPEVER', 'TBD')  # banzai pipeline version - Not banzai-nres.
        p_hdu_header.set('L1PUBDAT', 'TBD')  # when data will be made public.
        # Saving changes to the test files.
        hdu_list.writeto(filename, overwrite=True)
        hdu_list.close()
    print('finished patching keys to test fits files')
    # End of patching extravaganza.
    master_bias_path_and_filename = make_master_bias(test_image_context)
    test_master_bias = fits.getdata('~/' + master_bias_path_and_filename)
    print(test_master_bias)
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
    print(output_files)
    # End of Monkey Patch
    return output_files
