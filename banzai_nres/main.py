"""
main.py: Main driver script for the banzai-NRES pipeline.
    The main() function is a console entry point.
Authors
    Curtis McCully (cmccully@lcogt.net)
July 2018
    G. Mirek Brandt (gmbrandt@ucsb.edu)
July 2018
"""
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

from banzai_nres.utils.image_utils import read_images

from banzai import bias, trim, dark, gain
from banzai.qc import header_checker
from banzai import logs
from banzai.utils import image_utils
from banzai.main import get_stages_todo

from banzai.dbs import create_db, add_or_update_record, get_session, Site

logger = logs.get_logger(__name__)

ordered_stages = [header_checker.HeaderSanity,
                  bias.OverscanSubtractor,
                  gain.GainNormalizer,
                  trim.Trimmer,
                  bias.BiasSubtractor,
                  ]

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
    def __init__(self, filename, raw_path='/archive/engineering/lsc/nres01/20180313/raw'):
        _DEFAULT_DB = 'sqlite:////archive/engineering/test.db' #  from docker-compose file
        create_db('/archive/engineering/lsc/nres01/20180328/raw', db_address=_DEFAULT_DB,
                  configdb_address='http://configdb.lco.gtn/sites/')
        self.processed_path = '/tmp'
        self.raw_path = raw_path
        self.filename = filename
        self.post_to_archive = False
        self.db_address = _DEFAULT_DB
        self.preview_mode = False
        self.rlevel = 0
        self.fpack = True


def amend_nres_frames(pipeline_context, image_types = []):
    """
    Parameters:
        pipeline_context: pipeline context with a database already initialized.
        image_types: ['BIAS','DARK' etc...]
    This amends NRES frames to be able to pass through Banzai reduction steps.
    """
    image_list = image_utils.make_image_list(pipeline_context)
    image_list = image_utils.select_images(image_list, image_types)

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
                print('BPM key does not exist for %s \n creating and appending a zeros bad pixel mask' % filename)
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
        #p_hdu_header.set('CRPIX1', 0)
        #p_hdu_header.set('CRPIX2', 0)
        # Saving changes to the test files.
        hdu_list.writeto(filename, overwrite=True)
        hdu_list.close()
    print('finished patching keys to test fits files')


def parse_end_of_night_command_line_arguments():
    """
    :return: Directory where test NRES frames live. Eventually this would be hooked up to the
    pipeline, instead of giving a fixed directory.
    """
    return TestContext(filename=None,raw_path='/archive/engineering/lsc/nres01/20180313/raw')


def run_end_of_night_from_console(scripts_to_run):
    pipeline_context = parse_end_of_night_command_line_arguments()
    for script in scripts_to_run:
        script(pipeline_context)


def make_master_bias_console():
    run_end_of_night_from_console([make_master_bias])


def make_master_dark_console():
    run_end_of_night_from_console([make_master_dark])


def make_master_bias(pipeline_context):
    """
    Returns: None
    makes the master bias and saves the images.
    Note:
    image_types = ['BIAS'] selects only images which are bias type, naturally.
    """
    amend_nres_frames(pipeline_context, image_types=['BIAS'])

    stages_to_do = get_stages_todo(trim.Trimmer, extra_stages=[bias.BiasMaker])
    output_files = run(stages_to_do, pipeline_context, image_types=['BIAS'], calibration_maker=True,
        log_message='Making Master BIAS')
    return output_files


def make_master_dark(pipeline_context):
    amend_nres_frames(pipeline_context, image_types=['DARK'])

    stages_to_do = get_stages_todo(bias.BiasSubtractor, extra_stages=[dark.DarkNormalizer, dark.DarkMaker])
    run(stages_to_do, pipeline_context, image_types=['DARK'], calibration_maker=True,
        log_message='Making Master DARK')


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
    return output_files
