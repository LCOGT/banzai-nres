import pytest
import numpy as np

from banzai_nres.utils.array_utils import find_nearest


def test_find_nearest():
    array = np.arange(10)
    assert np.isclose(find_nearest(array, value=9.8), 9)
    assert np.isclose(find_nearest(array, value=7.2), 7)
    assert np.isclose(find_nearest(array, value=7.8), 8)
