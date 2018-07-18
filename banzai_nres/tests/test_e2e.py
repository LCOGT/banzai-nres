from __future__ import absolute_import, division, print_function, unicode_literals

import pytest

import numpy as np

from astropy.io import fits

from banzai import logs
from banzai.utils import image_utils

from banzai.dbs import create_db, add_or_update_record

from banzai_nres.main import make_master_bias

logger = logs.get_logger(__name__)

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

@pytest.mark.e2e
def test_making_master_biases():
    """
    Returns:
        Shape of master bias frame. For now this test is simple.

    """
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
    if test_image_context.fpack:
        master_bias_path_and_filename = str(make_master_bias(test_image_context)[0] + '.fz')
    else:
        master_bias_path_and_filename = str(make_master_bias(test_image_context)[0])
    test_master_bias = fits.getdata(master_bias_path_and_filename)
    print(test_master_bias.shape)
    assert test_master_bias.shape is not None


@pytest.mark.e2e
def test_e2e():
    assert True
