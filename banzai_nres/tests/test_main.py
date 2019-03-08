from banzai_nres.main import ReductionCriterion
from banzai.tests.utils import FakeContext
from banzai_nres.settings import NRESSettings
import datetime


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
