import pytest
from banzai.dbs import create_db, populate_bpm_table
import os
import numpy as np
import shutil
from astropy.io import fits
from banzai_nres.main import make_master_bias_console, make_master_dark_console, TestContext


def setup_module(module):
    create_db('./', db_address=os.environ['DB_URL'],
              configdb_address='http://configdbdev.lco.gtn/sites/')

    # clearing the bad pixel mask folder if it exists
    if os.path.exists('/archive/engineering/lsc/nres01/bpm'):
        shutil.rmtree('/archive/engineering/lsc/nres01/bpm')
    if os.path.exists('/archive/engineering/elp/nres02/bpm'):
        shutil.rmtree('/archive/engineering/elp/nres02/bpm')

    # building bpm folder
    os.makedirs('/archive/engineering/lsc/nres01/bpm')
    os.makedirs('/archive/engineering/elp/nres02/bpm')

    # adding the bpm folder to database
    populate_bpm_table('/archive/engineering/lsc/nres01/bpm', db_address=os.environ['DB_URL'])
    populate_bpm_table('/archive/engineering/elp/nres02/bpm', db_address=os.environ['DB_URL'])

    # using a bias fits as a template for the bpm fits.
    fits_file_to_copy = '/archive/engineering/lsc/nres01/20180328/raw/lscnrs01-fl09-20180328-0066-b00.fits.fz'
    date_marker = '20180725'

    # generating the zeros bpm.
    with fits.open(fits_file_to_copy) as hdu_list:
        hdu_list[1].data = np.zeros(hdu_list[1].data.shape, dtype=np.uint8)
        hdu_list[1].header['OBSTYPE'] = 'BPM'
        hdu_list.writeto('/archive/engineering/lsc/nres01/bpm/lsc_fl09_BPM_' +
                         date_marker + '.fits.fz', overwrite=True)
        hdu_list[1].header['INSTRUME'] = 'fl17'
        hdu_list.writeto('/archive/engineering/elp/nres02/bpm/elp_fl17_BPM_'
                         + date_marker + '.fits.fz', overwrite=True)

@pytest.mark.e2e
def test_e2e():
    make_master_bias_console()
    make_master_dark_console()

    test_context = TestContext()
    instrument = 'nres01'
    epoch = test_context.raw_path[-12:-4]
    site = 'lsc'
    expected_dark_filename = 'dark_' + instrument + '_' + epoch + '_bin1x1.fits.fz'
    expected_bias_filename = 'bias_' + instrument + '_' + epoch + '_bin1x1.fits.fz'
    expected_processed_path = os.path.join(test_context.processed_path, site,
                                           instrument, epoch, 'processed')

    with fits.open(os.path.join(expected_processed_path, expected_bias_filename)) as hdu_list:
        assert hdu_list[1].data.shape is not None
        assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape

    with fits.open(os.path.join(expected_processed_path, expected_dark_filename)) as hdu_list:
        assert hdu_list[1].data.shape is not None
        assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape
