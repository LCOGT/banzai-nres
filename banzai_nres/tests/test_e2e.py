import pytest

from banzai import logs


logger = logs.get_logger(__name__)


@pytest.mark.e2e
def test_e2e():
    assert True
