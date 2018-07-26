import pytest

from banzai import logs
from banzai_nres.main import make_master_bias_console


logger = logs.get_logger(__name__)


@pytest.mark.e2e
def test_e2e():
    assert True

def test_making_master_bias():
    make_master_bias_console()
    assert True