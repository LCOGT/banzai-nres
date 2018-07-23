from __future__ import absolute_import, division, print_function, unicode_literals

import pytest

import numpy as np

from astropy.io import fits

from banzai import logs


logger = logs.get_logger(__name__)


@pytest.mark.e2e
def test_e2e():
    assert True
