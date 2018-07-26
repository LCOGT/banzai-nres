import pytest
from banzai.dbs import create_db, populate_bpm_table
import os
import datetime
import numpy as np
from astropy.io import fits
from  banzai_nres.main import make_master_bias_console, make_master_dark_console

def setup_module(module):
    create_db('./', db_address=os.environ['DB_URL'],
              configdb_address='http://configdb.lco.gtn/sites/')

    os.makedirs('/archive/engineering/lsc/nres01/bpm')
    os.makedirs('/archive/engineering/elp/nres02/bpm')

    populate_bpm_table('/archive/engineering/lsc/nres01/bpm', db_address=os.environ['DB_URL'])
    populate_bpm_table('/archive/engineering/elp/nres02/bpm', db_address=os.environ['DB_URL'])

    fits_file_to_copy = '/archive/engineering/lsc/nres01/20180328/raw/lscnrs01-fl09-20180328-0066-b00.fits.fz'
    now = datetime.datetime.now()
    date_marker = str(now.month) + '-' + str(now.day) + '-' + str(now.year)

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
    assert True

# check shape of loaded fits. check the fits key BPM, [1].data.