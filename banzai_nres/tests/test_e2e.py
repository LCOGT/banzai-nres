import pytest
from banzai.tests.utils import FakeResponse, get_min_and_max_dates
from banzai_nres import settings
from banzai.logs import get_logger
import os
import mock
import numpy as np
from banzai.utils import file_utils
import time
from glob import glob
from banzai.celery import app, schedule_calibration_stacking
from banzai.dbs import get_session
from banzai import dbs
from types import ModuleType
from datetime import datetime
from dateutil.parser import parse
from astropy.io import fits, ascii
from astropy.table import Table
from banzai.utils import fits_utils
import banzai_nres.dbs
import json
import importlib.resources

logger = get_logger()

TEST_PACKAGE = 'banzai_nres.tests'


TEST_PACKAGE = 'banzai_nres.tests'
TEST_FRAMES = ascii.read(os.path.join(importlib.resources.files(TEST_PACKAGE), 'data/test_data.dat'))
BPM_FILES = ascii.read(os.path.join(importlib.resources.files(TEST_PACKAGE), 'data/test_bpm.dat'))

DATA_ROOT = os.path.join(os.sep, 'archive', 'engineering')

SITES = set([frame[:3] for frame in TEST_FRAMES['filename']])
INSTRUMENTS = set([os.path.join(frame[:3], frame.split('-')[1]) for frame in TEST_FRAMES['filename']])

DAYS_OBS = set([os.path.join(frame[:3], frame.split('-')[1], frame.split('-')[2]) for frame in TEST_FRAMES['filename']])

CONFIGDB_FILENAME = os.path.join(importlib.resources.files('banzai_nres.tests'), 'data', 'configdb_example.json')
PHOENIX_FILENAME = os.path.join(importlib.resources.files('banzai_nres.tests'), 'data', 'phoenix.json')


def observation_portal_side_effect(*args, **kwargs):
    site = kwargs['params']['site']
    start = datetime.strftime(parse(kwargs['params']['start_after']).replace(tzinfo=None).date(), '%Y%m%d')
    filename = 'test_obs_portal_response_{site}_{start}.json'.format(site=site, start=start)
    filename = os.path.join(importlib.resources.files('banzai_nres.tests'), 'data', f'{filename}')
    return FakeResponse(filename)


def get_instrument_ids(db_address, names):
    with get_session(db_address) as db_session:
        instruments = []
        for name in names:
            criteria = dbs.Instrument.name == name
            instruments.extend(db_session.query(dbs.Instrument).filter(criteria).all())
    return [instrument.id for instrument in instruments]


def celery_join():
    celery_inspector = app.control.inspect()
    log_counter = 0
    while True:
        time.sleep(1)
        log_counter += 1
        if log_counter % 5 == 0:
            logger.info('Processing: ' + '. ' * (log_counter // 5))
        queues = [celery_inspector.active(), celery_inspector.scheduled(), celery_inspector.reserved()]
        if any([queue is None or 'celery@banzai-celery-worker' not in queue for queue in queues]):
            logger.warning('No valid celery queues were detected, retrying...', extra_tags={'queues': queues})
            # Reset the celery connection
            celery_inspector = app.control.inspect()
            continue
        if all([len(queue['celery@banzai-celery-worker']) == 0 for queue in queues]):
            break


def run_reduce_individual_frames(filename_pattern):
    logger.info('Reducing individual frames for filenames: {filenames}'.format(filenames=filename_pattern))
    for frame in TEST_FRAMES:
        if filename_pattern in frame['filename']:
            file_utils.post_to_archive_queue(frame['filename'], frame['frameid'],
                                             os.getenv('FITS_BROKER'),
                                             exchange_name=os.getenv('FITS_EXCHANGE'),
                                             SITEID=frame['site'], INSTRUME=frame['instrument'])
    celery_join()
    logger.info('Finished reducing individual frames for filenames: {filenames}'.format(filenames=filename_pattern))


def stack_calibrations(frame_type):
    logger.info('Stacking calibrations for frame type: {frame_type}'.format(frame_type=frame_type))
    for day_obs in DAYS_OBS:
        site, camera, dayobs = day_obs.split('/')
        timezone = dbs.get_timezone(site, db_address=os.environ['DB_ADDRESS'])
        min_date, max_date = get_min_and_max_dates(timezone, dayobs=dayobs)
        runtime_context = dict(processed_path=DATA_ROOT, log_level='debug', post_to_archive=False,
                               post_to_opensearch=False, fpack=True, reduction_level=92,
                               db_address=os.environ['DB_ADDRESS'], opensearch_qc_index='banzai_qc',
                               opensearch_url='https://opensearch.lco.global',
                               no_bpm=False, ignore_schedulability=True, use_only_older_calibrations=False,
                               preview_mode=False, max_tries=5, broker_url=os.getenv('FITS_BROKER'),
                               no_file_cache=False)
        for setting in dir(settings):
            if '__' != setting[:2] and not isinstance(getattr(settings, setting), ModuleType):
                runtime_context[setting] = getattr(settings, setting)
        schedule_calibration_stacking(site, runtime_context, min_date=min_date, max_date=max_date,
                                      frame_types=[frame_type])
    celery_join()
    logger.info('Finished stacking calibrations for frame type: {frame_type}'.format(frame_type=frame_type))


def mark_frames_as_good(raw_filenames):
    logger.info('Marking frames as good for filenames: {filenames}'.format(filenames=raw_filenames))
    for day_obs in DAYS_OBS:
        for filename in glob(os.path.join(DATA_ROOT, day_obs, 'processed', raw_filenames)):
            dbs.mark_frame(os.path.basename(filename), "good", db_address=os.environ['DB_ADDRESS'])
    logger.info('Finished marking frames as good for filenames: {filenames}'.format(filenames=raw_filenames))


def get_expected_number_of_calibrations(raw_filename_pattern, calibration_type):
    number_of_stacks_that_should_have_been_created = 0
    for day_obs in DAYS_OBS:
        site, instrument, dayobs = day_obs.split('/')
        raw_frames_for_this_dayobs = [
            frame for frame in TEST_FRAMES
            if site in frame['filename'] and instrument in frame['filename']
            and dayobs in frame['filename'] and raw_filename_pattern in frame['filename']
        ]
        if calibration_type.lower() == 'lampflat' or calibration_type.lower() == 'double':
            # Group by fibers lit if we are stacking lampflats or doubles (arc frames)
            observed_fibers = []
            for frame in raw_frames_for_this_dayobs:
                observed_fibers.append(frame['OBJECTS'])
            observed_fibers = set(observed_fibers)
            number_of_stacks_that_should_have_been_created += len(observed_fibers)
        else:
            # Just one calibration per night
            if len(raw_frames_for_this_dayobs) > 0:
                number_of_stacks_that_should_have_been_created += 1
    return number_of_stacks_that_should_have_been_created


def check_if_individual_frames_exist(filename_pattern):
    for frame in TEST_FRAMES:
        if filename_pattern in frame['filename']:
            processed_path = os.path.join(DATA_ROOT, frame['site'], frame['instrument'], frame['dayobs'], 'processed')
            assert os.path.exists(processed_path, frame['filename'].replace('00', '92'))


def run_check_if_stacked_calibrations_were_created(raw_filenames, calibration_type):
    created_stacked_calibrations = []
    number_of_stacks_that_should_have_been_created = get_expected_number_of_calibrations(raw_filenames,
                                                                                         calibration_type)
    for day_obs in DAYS_OBS:
        created_stacked_calibrations += glob(os.path.join(DATA_ROOT, day_obs, 'processed',
                                                          '*' + calibration_type.lower() + '*.fits*'))
    assert number_of_stacks_that_should_have_been_created > 0
    assert len(created_stacked_calibrations) == number_of_stacks_that_should_have_been_created


def run_check_if_stacked_calibrations_have_extensions(calibration_type, extensions_to_check):
    created_stacked_calibrations = []
    for day_obs in DAYS_OBS:
        created_stacked_calibrations += glob(os.path.join(DATA_ROOT, day_obs, 'processed',
                                                          '*' + calibration_type.lower() + '*.fits*'))
    for cal in created_stacked_calibrations:
        hdulist = fits.open(cal)
        extnames = [hdulist[i].header.get('extname', None) for i in range(len(hdulist))]
        for ext in extensions_to_check:
            logger.info(f'checking if {ext} is in the saved extensions of {cal}')
            assert ext in extnames


def check_extracted_spectra(raw_filename, spec_extname, columns):
    created_images = []
    for day_obs in DAYS_OBS:
        created_images += glob(os.path.join(DATA_ROOT, day_obs, 'processed', raw_filename))
    for filename in created_images:
        with fits.open(filename) as f:
            hdu = fits_utils.unpack(f)
        spectrum = Table(hdu[spec_extname].data)
        for colname in columns:
            assert colname in spectrum.colnames
            assert not np.allclose(spectrum[colname], 0)
        assert 'RV' in hdu[0].header


def run_check_if_stacked_calibrations_are_in_db(raw_filenames, calibration_type):
    number_of_stacks_that_should_have_been_created = get_expected_number_of_calibrations(raw_filenames,
                                                                                         calibration_type)
    with dbs.get_session(os.environ['DB_ADDRESS']) as db_session:
        calibrations_in_db = db_session.query(dbs.CalibrationImage).filter(
            dbs.CalibrationImage.type == calibration_type)
        calibrations_in_db = calibrations_in_db.filter(dbs.CalibrationImage.is_master).all()
    assert number_of_stacks_that_should_have_been_created > 0
    assert len(calibrations_in_db) == number_of_stacks_that_should_have_been_created


def mock_phoenix_models_in_db(db_address):
    with open(PHOENIX_FILENAME) as f:
        phoenix_data = json.load(f)
    with dbs.get_session(db_address) as db_session:
        db_session.bulk_insert_mappings(banzai_nres.dbs.PhoenixModel, phoenix_data)
        dbs.add_or_update_record(db_session, banzai_nres.dbs.ResourceFile, {'key': 'phoenix_wavelengths'},
                                 {'filename': 'phoenix_wavelength.fits',
                                  'location': 's3://banzai-nres-phoenix-models-lco-global',
                                  'key': 'phoenix_wavelengths'})


@pytest.mark.e2e
@pytest.fixture(scope='module')
@mock.patch('banzai.dbs.requests.get', return_value=FakeResponse(CONFIGDB_FILENAME))
def init(configdb):
    os.system(f'banzai_nres_create_db --db-address={os.environ["DB_ADDRESS"]}')
    dbs.populate_instrument_tables(db_address=os.environ["DB_ADDRESS"], configdb_address='http://fakeconfigdb')
    os.system((f'banzai_add_site --site elp --latitude 30.67986944 --longitude -104.015175'
              f' --elevation 2027 --timezone -6 --db-address={os.environ["DB_ADDRESS"]}'))
    os.system((f'banzai_add_site --site lsc --latitude -30.1673833333 --longitude -70.8047888889'
              f' --elevation 2198 --timezone -4 --db-address={os.environ["DB_ADDRESS"]}'))
    os.system((f'banzai_add_instrument --site lsc --camera fl09 --name nres01 --nx 4096 --ny 4096'
              f' --instrument-type 1m0-NRES-SciCam --db-address={os.environ["DB_ADDRESS"]}'))
    os.system((f'banzai_add_instrument --site elp --camera fl17 --name nres02 --nx 4096 --ny 4096'
              f' --instrument-type 1m0-NRES-SciCam --db-address={os.environ["DB_ADDRESS"]}'))

    mock_phoenix_models_in_db(os.environ["DB_ADDRESS"])
    for frame in BPM_FILES:
        logger.info(f'adding bpm {frame["filename"]} to the database')
        instrument = dbs.query_for_instrument(camera=frame['instrument'],
                                              site=frame['site'],
                                              db_address=os.environ['DB_ADDRESS'])
        calimage = dbs.CalibrationImage(
            type=frame['obstype'],
            filename=frame['filename'],
            frameid=f'{frame["frameid"]:d}',
            dateobs=datetime.strptime(frame['dateobs'], '%Y-%m-%d'),
            datecreated=datetime(2023, 11, 19),
            instrument_id=instrument.id,
            is_master=True, is_bad=False,
            attributes={'binning': frame['binning'], 'configuration_mode': frame['mode']}
        )
        with get_session(os.environ['DB_ADDRESS']) as db_session:
            db_session.add(calimage)
            db_session.commit()


@pytest.mark.e2e
@pytest.mark.master_bias
class TestMasterBiasCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.observation_utils.requests.get', side_effect=observation_portal_side_effect)
    def stack_bias_frames(self, mock_lake, init):
        run_reduce_individual_frames('b00.fits')
        mark_frames_as_good('*b92.fits*')
        stack_calibrations('bias')

    def test_if_stacked_bias_frame_was_created(self):
        check_if_individual_frames_exist('b00')
        run_check_if_stacked_calibrations_were_created('b00.fits', 'bias')
        run_check_if_stacked_calibrations_are_in_db('*b00.fits*', 'BIAS')


@pytest.mark.e2e
@pytest.mark.master_dark
class TestMasterDarkCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.observation_utils.requests.get', side_effect=observation_portal_side_effect)
    def stack_dark_frames(self, mock_lake):
        run_reduce_individual_frames('d00.fits')
        mark_frames_as_good('*d92.fits*')
        stack_calibrations('dark')

    def test_if_stacked_dark_frame_was_created(self):
        check_if_individual_frames_exist('d00')
        run_check_if_stacked_calibrations_were_created('d00.fits', 'dark')
        run_check_if_stacked_calibrations_are_in_db('*d00.fits*', 'DARK')


@pytest.mark.e2e
@pytest.mark.master_flat
class TestMasterFlatCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.observation_utils.requests.get', side_effect=observation_portal_side_effect)
    def stack_flat_frames(self, mock_lake):
        run_reduce_individual_frames('w00.fits')
        mark_frames_as_good('*w92.fits*')
        stack_calibrations('lampflat')

    def test_if_stacked_flat_frame_was_created(self):
        check_if_individual_frames_exist('w00')
        run_check_if_stacked_calibrations_were_created('w00.fits', 'lampflat')
        run_check_if_stacked_calibrations_have_extensions('lampflat', ['TRACES', 'PROFILE', 'BLAZE'])
        run_check_if_stacked_calibrations_are_in_db('*w00.fits*', 'LAMPFLAT')


@pytest.mark.e2e
@pytest.mark.master_arc
class TestMasterArcCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.observation_utils.requests.get', side_effect=observation_portal_side_effect)
    def stack_arc_frames(self, mock_lake):
        run_reduce_individual_frames('a00.fits')
        mark_frames_as_good('*a92.fits*')
        stack_calibrations('double')

    def test_if_stacked_arc_frame_was_created(self):
        check_if_individual_frames_exist('a00')
        run_check_if_stacked_calibrations_were_created('a00.fits', 'double')
        run_check_if_stacked_calibrations_have_extensions('double', ['WAVELENGTH', 'FEATURES'])
        run_check_if_stacked_calibrations_are_in_db('*a00.fits*', 'DOUBLE')

    def test_quality_of_wavelength_calibration(self, calibration_type='double', primaryextension=1):
        created_stacked_calibrations = []
        for day_obs in DAYS_OBS:
            created_stacked_calibrations += glob(os.path.join(DATA_ROOT, day_obs, 'processed',
                                                              '*' + calibration_type.lower() + '*.fits*'))
        for cal in created_stacked_calibrations:
            hdulist = fits.open(cal)
            quality_metrics = hdulist[primaryextension].header
            assert quality_metrics['RVPRECSN'] < 10
            assert quality_metrics['RVPRECSN'] > 1


@pytest.mark.e2e
@pytest.mark.science_frames
class TestScienceFrameProcessing:
    @pytest.fixture(autouse=True)
    # Note this requires the GAIA and SIMBAD services to be up. It's a little scary to depend on outside data source
    # for our tests. To mock this, we would have to write a control command and use broadcast() to get it to the workers
    # See https://stackoverflow.com/questions/30450468/mocking-out-a-call-within-a-celery-task
    def process_frames(self):
        run_reduce_individual_frames('*e00.fits*')

    def test_if_science_frames_were_created(self):
        for frame in TEST_FRAMES:
            if 'e00.fits' in frame['filename']:
                processed_path = os.path.join(DATA_ROOT, frame['site'], frame['instrument'],
                                              frame['dayobs'], 'processed')
                assert os.path.exists(processed_path, frame['filename'].replace('00', '92'))
                assert os.path.exists(processed_path, frame['filename'].replace('00', '92-1d'))
                assert os.path.exists(processed_path, frame['filename'].replace('00', '92-2d'))
                summary_filename = processed_path, frame['filename'].replace('00.fits', '92.fits').replace('.fz', '')
                assert os.path.exists(summary_filename)

        check_extracted_spectra('*e92-1d.fits*', 'SPECTRUM', ['wavelength', 'flux', 'uncertainty'])
