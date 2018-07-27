import pytest
from banzai.dbs import create_db, populate_bpm_table
import os
import numpy as np
import shutil
from astropy.io import fits
from banzai_nres.main import make_master_bias_console, make_master_dark_console, TestContext



def setup_module(module):
    create_db('./', db_address=os.environ['DB_URL'],
              configdb_address=os.environ['CONFIG_DB_URL'])

    # clearing the bad pixel mask folder if it exists
    if os.path.exists('/archive/engineering/lsc/nres01/bpm'):
        shutil.rmtree('/archive/engineering/lsc/nres01/bpm')
    if os.path.exists('/archive/engineering/elp/nres02/bpm'):
        shutil.rmtree('/archive/engineering/elp/nres02/bpm')

    # building bpm folder
    os.makedirs('/archive/engineering/lsc/nres01/bpm')
    os.makedirs('/archive/engineering/elp/nres02/bpm')

    # using an arbitrary fits as a template for the bpm fits.
    fits_file_to_copy = '/archive/engineering/lsc/nres01/20180228/raw/lscnrs01-fl09-20180228-0010-e00.fits'
    date_marker = '20180727'

    os.system('funpack {file}'.format(file=fits_file_to_copy + '.fz'))
    # generating the zeros bpm. Files need to start with bpm.
    with fits.open(fits_file_to_copy) as hdu_list:
        hdu_list[1].data = np.zeros(hdu_list[1].data.shape, dtype=np.uint8)
        hdu_list[1].header['OBSTYPE'] = 'BPM'
        hdu_list[1].header['EXTNAME'] = 'BPM'
        hdu_list[1].header['INSTRUME'] = 'nres01'
        output_filename = '/archive/engineering/lsc/nres01/bpm/bpm_lsc_fl09_' \
                          + date_marker + '.fits'
        hdu_list.writeto(output_filename, overwrite=True)
        # fpack the file (which creates an fpacked copy)
        os.system('fpack {file}'.format(file=output_filename))
        # delete the un-fpacked file.
        os.system('rm {file}'.format(file=output_filename))

        hdu_list[1].header['INSTRUME'] = 'nres02'
        hdu_list[1].header['SITEID'] = 'elp'
        hdu_list[1].header['TELESCOP'] = 'nres02'
        output_filename = '/archive/engineering/elp/nres02/bpm/bpm_elp_fl17_' \
                          + date_marker + '.fits'
        hdu_list.writeto(output_filename, overwrite=True)
        os.system('fpack {file}'.format(file=output_filename))
        os.system('rm {file}'.format(file=output_filename))

    # delete the unpacked file which was initially copied.
    os.system('rm {file}'.format(file=fits_file_to_copy))
    # adding the bpm folder to database and populating the sqlite tables.
    populate_bpm_table('/archive/engineering/lsc/nres01/bpm', db_address=os.environ['DB_URL'])
    populate_bpm_table('/archive/engineering/elp/nres02/bpm', db_address=os.environ['DB_URL'])

@pytest.mark.e2e
def test_e2e():
    test_context = TestContext()
    instrument = 'nres01'
    epoch = test_context.raw_path[-12:-4]
    site = 'lsc'
    expected_dark_filename = 'dark_' + instrument + '_' + epoch + '_bin1x1.fits.fz'
    expected_bias_filename = 'bias_' + instrument + '_' + epoch + '_bin1x1.fits.fz'
    expected_processed_path = os.path.join(test_context.processed_path, site,
                                           instrument, epoch, 'processed')

    make_master_bias_console()
    with fits.open(os.path.join(expected_processed_path, expected_bias_filename)) as hdu_list:
        assert hdu_list[1].data.shape is not None
        assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape

    make_master_dark_console()
    with fits.open(os.path.join(expected_processed_path, expected_dark_filename)) as hdu_list:
        assert hdu_list[1].data.shape is not None
        assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape
