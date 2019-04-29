import mock
import numpy as np
import tempfile
import os
import pytest

from banzai_nres.utils import db_utils
from banzai.tests.utils import FakeResponse
from banzai import dbs

from banzai.tests.utils import FakeContext, FakeImage


@mock.patch('banzai_nres.utils.db_utils.dbs.get_timezone', return_value='UTC')
@mock.patch('banzai_nres.utils.db_utils.date_utils.get_dayobs', return_value='42')
def test_get_raw_path(mockday, mockzone):
    fake_context = FakeContext()
    fake_context.site = 'site'
    fake_context.instrument_name = 'nres'
    new_path = db_utils.get_raw_path(base_raw_path='base/', runtime_context=fake_context)
    assert new_path == 'base/site/nres/42/raw'
    assert db_utils.get_raw_path(base_raw_path=None, runtime_context=fake_context) is None
