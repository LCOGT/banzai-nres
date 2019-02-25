import pytest
from banzai.dbs import create_db, populate_calibration_table_with_bpms
from banzai_nres.settings import NRESSettings
import os
import numpy as np
import shutil
from astropy.io import fits


def make_dummy_bpm(bpm_path, output_bpm_name_addition, fits_file_to_copy, date_marker, telescope_name, site_name):
    """
    Creates and saves a dummy bpm in the format of a real fits file.
    """
    # clearing the bad pixel mask folder if it exists
    if os.path.exists(bpm_path):
        shutil.rmtree(bpm_path)
    # building bpm folder
    os.makedirs(bpm_path)
    # unpacking a fits file via funpack. Astropy's unpack messes with the files.
    os.system('funpack {}'.format(fits_file_to_copy + '.fz'))
    # where to save the file
    output_filename = bpm_path + output_bpm_name_addition + date_marker + '.fits'
    # creating the bpm
    with fits.open(fits_file_to_copy) as hdu_list:
        hdu_list[0].data = np.zeros(hdu_list[0].data.shape, dtype=np.uint8)
        hdu_list[0].header['OBSTYPE'] = 'BPM'
        hdu_list[0].header['EXTNAME'] = 'BPM'
        hdu_list[0].header['INSTRUME'] = telescope_name
        hdu_list[0].header['SITEID'] = site_name
        hdu_list[0].header['TELESCOP'] = telescope_name
        hdu_list.writeto(output_filename, overwrite=True)

    # fpack the file and delete the funpacked input.
    os.system('fpack {0}'.format(output_filename))
    os.system('rm {0}'.format(output_filename))
    # delete the unpacked file which was initially copied into raw/ via funpack
    os.system('rm {0}'.format(fits_file_to_copy))


def setup_module(module):
    """
    :param module: Pytest placeholder argument.

    This function creates the sqlite database and populates it with
    telescopes and BPM's for the test data sets elp/nres02 and lsc/nres01.
    """
    create_db('./', db_address=os.environ['DB_URL'],
              configdb_address=os.environ['CONFIG_DB_URL'])

    # using an arbitrary fits as a template for the bpm fits. Then making and saving the bpm's
    fits_file_to_copy = '/archive/engineering/lsc/nres01/20180228/raw/lscnrs01-fl09-20180228-0010-e00.fits'
    date_marker = '20180727'

    make_dummy_bpm('/archive/engineering/lsc/nres01/bpm', '/bpm_lsc_fl09_',
                   fits_file_to_copy=fits_file_to_copy, date_marker=date_marker,
                   telescope_name='nres01', site_name='lsc')
    make_dummy_bpm('/archive/engineering/elp/nres02/bpm', '/bpm_elp_fl17_',
                   fits_file_to_copy=fits_file_to_copy, date_marker=date_marker,
                   telescope_name='nres02', site_name='elp')

    # adding the bpm folder to database and populating the sqlite tables.
    populate_calibration_table_with_bpms('/archive/engineering/lsc/nres01/bpm', db_address=os.environ['DB_URL'])
    populate_calibration_table_with_bpms('/archive/engineering/elp/nres02/bpm', db_address=os.environ['DB_URL'])


@pytest.mark.e2e
def test_e2e():
    db_address = os.environ['DB_URL']
    raw_data_path = '/archive/engineering/lsc/nres01/20180311/raw'
    instrument = 'nres01'
    site = 'lsc'
    epoch = '20180311'

    expected_bias_filename = 'lscnrs01-fl09-20180311-bias-bin1x1.fits'
    expected_dark_filename = 'lscnrs01-fl09-20180311-dark-bin1x1.fits'
    expected_flat_filenames = ['lscnrs01-fl09-20180311-lampflat-bin1x1-110.fits',
                               'lscnrs01-fl09-20180311-lampflat-bin1x1-011.fits']
    expected_trace_filenames = ['lscnrs01-fl09-20180311-trace-bin1x1-110.fits',
                                'lscnrs01-fl09-20180311-trace-bin1x1-011.fits']
    expected_processed_path = os.path.join('/tmp', site, instrument, epoch, 'processed')

    # executing the master bias maker as one would from the command line.
    os.system('reduce_bias_frames --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))
    os.system('stack_calibrations --site lsc --camera nres01 --frame-type BIAS --min-date 2018-03-11T00:00:00'
              ' --max-date 2018-03-12T23:59:59 --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))

    with fits.open(os.path.join(expected_processed_path, expected_bias_filename)) as hdu_list:
        assert hdu_list[0].data.shape is not None
        assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape

    # executing the master dark maker as one would from the command line.
    os.system('reduce_dark_frames --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))
    os.system('stack_calibrations --site lsc --camera nres01 --frame-type DARK --min-date 2018-03-11T00:00:00'
              ' --max-date 2018-03-12T23:59:59 --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))

    with fits.open(os.path.join(expected_processed_path, expected_dark_filename)) as hdu_list:
        assert hdu_list[0].data.shape is not None
        assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape

    # executing the master flat maker as one would from the command line.
    os.system('reduce_flat_frames --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))
    os.system('stack_calibrations --site lsc --camera nres01 --frame-type LAMPFLAT --min-date 2018-03-11T00:00:00'
              ' --max-date 2018-03-12T23:59:59 --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))

    for expected_flat_filename in expected_flat_filenames:
        with fits.open(os.path.join(expected_processed_path, expected_flat_filename)) as hdu_list:
            assert hdu_list[0].data.shape is not None
            assert hdu_list['BPM'].data.shape == hdu_list[1].data.shape

    # executing the master trace maker as one would from the command line
    trace_table_name = NRESSettings.TRACE_TABLE_NAME
    os.system('stack_calibrations --site lsc --camera nres01 --frame-type TRACE --min-date 2018-03-11T00:00:00'
              ' --max-date 2018-03-12T23:59:59 --db-address {0} --raw-path {1} --ignore-schedulability '
              '--processed-path /tmp --log-level debug'.format(db_address, raw_data_path))

    for filename in expected_trace_filenames:
        with fits.open(os.path.join(expected_processed_path, filename)) as hdu_list:
            assert hdu_list[trace_table_name] is not None
            assert hdu_list[trace_table_name].data['centers'].shape[0] > 100
