import datetime
import pytest
import mock

from banzai_nres.main import ReductionCriterion
from banzai_nres.settings import NRESSettings

from banzai.tests.utils import FakeContext


def test_reduction_criterion_auto_fill():
    expected_max_date = datetime.datetime.utcnow()
    expected_min_date = expected_max_date - datetime.timedelta(hours=24)
    reduction_criterion = ReductionCriterion()
    assert reduction_criterion.frame_types == NRESSettings().REDUCE_NIGHT_FRAME_TYPES
    assert reduction_criterion.min_date > expected_min_date - datetime.timedelta(seconds=10)
    assert reduction_criterion.min_date < expected_min_date + datetime.timedelta(seconds=10)
    assert reduction_criterion.max_date > expected_max_date - datetime.timedelta(seconds=10)
    assert reduction_criterion.max_date < expected_max_date + datetime.timedelta(seconds=10)


def test_reduction_criterion_sets_min_date():
    fake_context = FakeContext(settings=NRESSettings())
    fake_context.max_date = datetime.datetime.utcnow()
    reduction_criterion = ReductionCriterion(pipeline_context=fake_context)
    assert reduction_criterion.min_date == fake_context.max_date - datetime.timedelta(hours=24)


def test_reduction_criterion_sets_max_date():
    fake_context = FakeContext(settings=NRESSettings())
    expected_max_date = datetime.datetime.utcnow()
    fake_context.min_date = datetime.datetime.utcnow()
    reduction_criterion = ReductionCriterion(pipeline_context=fake_context)
    assert reduction_criterion.max_date > expected_max_date - datetime.timedelta(seconds=10)
    assert reduction_criterion.max_date < expected_max_date + datetime.timedelta(seconds=10)


def test_reduction_criterion_outputs_single_frame_types():
    fake_context = FakeContext(settings=NRESSettings())
    fake_context.frame_type = 'TYPE'
    reduction_criterion = ReductionCriterion(pipeline_context=fake_context)
    assert reduction_criterion.frame_types == [fake_context.frame_type]


def test_reduction_criterion_keeps_good_raw_path():
    fake_context = FakeContext(settings=NRESSettings())
    fake_context.frame_type = 'TYPE'
    reduction_criterion = ReductionCriterion(raw_path='test/raw', pipeline_context=fake_context)
    assert reduction_criterion.raw_path == 'test/raw'


@mock.patch('banzai_nres.main.get_raw_path', return_value='test/raw')
def test_reduction_criterion_keeps_completes_raw_path(mock_get_raw_path):
    fake_context = FakeContext(settings=NRESSettings())
    fake_context.frame_type = 'TYPE'
    reduction_criterion = ReductionCriterion(raw_path='test/', pipeline_context=fake_context)
    assert reduction_criterion.raw_path == 'test/raw'


def test_raise_exception_if_min_date_greater_than_max_date():
    fake_context = FakeContext(settings=NRESSettings())
    fake_context.min_date = datetime.datetime.utcnow()
    fake_context.max_date = fake_context.min_date - datetime.timedelta(hours=1)
    with pytest.raises(Exception):
        ReductionCriterion(pipeline_context=fake_context)
