import pytest
import banzai_nres.settings as nres_settings
from banzai.tests.utils import FakeResponse
import os
import mock
from banzai import dbs
from banzai.tests import test_end_to_end
from glob import glob


import logging

logger = logging.getLogger(__name__)


DATA_ROOT = os.path.join(os.sep, 'archive', 'engineering')

SITES = [os.path.basename(site_path) for site_path in glob(os.path.join(DATA_ROOT, '???'))]
INSTRUMENTS = [os.path.join(site, os.path.basename(instrument_path)) for site in SITES
               for instrument_path in glob(os.path.join(os.path.join(DATA_ROOT, site, '*')))]

DAYS_OBS = [os.path.join(instrument, os.path.basename(dayobs_path)) for instrument in INSTRUMENTS
            for dayobs_path in glob(os.path.join(DATA_ROOT, instrument, '201*'))]


@pytest.mark.e2e
@pytest.fixture(scope='module')
@mock.patch('banzai.dbs.requests.get', return_value=FakeResponse('data/configdb_example.json'))
def init(configdb):
    dbs.create_db('.', db_address=os.environ['DB_ADDRESS'], configdb_address='http://configdbdev.lco.gtn/sites/')
    with dbs.get_session(db_address=os.environ['DB_ADDRESS']) as db_session:
        elp_nres = dbs.Instrument(name='nres02', camera='fl17', enclosure='igla',
                                  telescope='', type='NRES', schedulable=True, site='elp')
        lsc_nres = dbs.Instrument(name='nres01', camera='fl09', enclosure='igla',
                                  telescope='', type='NRES', schedulable=True, site='lsc')
        db_session.add(elp_nres)
        db_session.add(lsc_nres)
        db_session.commit()
    for instrument in INSTRUMENTS:
        dbs.populate_calibration_table_with_bpms(os.path.join(DATA_ROOT, instrument, 'bpm'),
                                                 db_address=os.environ['DB_ADDRESS'])


@pytest.mark.e2e
@pytest.mark.master_bias
class TestMasterBiasCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.lake_utils.requests.get', side_effect=test_end_to_end.lake_side_effect)
    def stack_bias_frames(self, mock_lake, init):
        test_end_to_end.run_reduce_individual_frames('*b00.fits*')
        test_end_to_end.mark_frames_as_good('*b91.fits*')
        test_end_to_end.run_stack_calibrations('bias')

    def test_if_stacked_bias_frame_was_created(self):
        test_end_to_end.run_check_if_stacked_calibrations_were_created('*b00.fits*', 'bias')
        test_end_to_end.run_check_if_stacked_calibrations_are_in_db('*b00.fits*', 'BIAS')


@pytest.mark.e2e
@pytest.mark.master_dark
class TestMasterDarkCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.lake_utils.requests.get', side_effect=test_end_to_end.lake_side_effect)
    def stack_dark_frames(self, mock_lake):
        test_end_to_end.run_reduce_individual_frames('*d00.fits*')
        test_end_to_end.mark_frames_as_good('*d91.fits*')
        test_end_to_end.run_stack_calibrations('dark')

    def test_if_stacked_dark_frame_was_created(self):
        test_end_to_end.run_check_if_stacked_calibrations_were_created('*d00.fits*', 'dark')
        test_end_to_end.run_check_if_stacked_calibrations_are_in_db('*d00.fits*', 'DARK')


@pytest.mark.e2e
@pytest.mark.master_flat
class TestMasterFlatCreation:
    @pytest.fixture(autouse=True)
    @mock.patch('banzai.utils.lake_utils.requests.get', side_effect=test_end_to_end.lake_side_effect)
    def stack_flat_frames(self, mock_lake):
        test_end_to_end.run_reduce_individual_frames('*w00.fits*')
        test_end_to_end.mark_frames_as_good('*w91.fits*')
        test_end_to_end.run_stack_calibrations('lampflat')

    def test_if_stacked_flat_frame_was_created(self):
        test_end_to_end.run_check_if_stacked_calibrations_were_created('*w00.fits*', 'lampflat')
        test_end_to_end.run_check_if_stacked_calibrations_are_in_db('*w00.fits*', 'LAMPFLAT')


# TODO add master traces
