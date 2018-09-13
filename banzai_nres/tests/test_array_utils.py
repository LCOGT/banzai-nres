import pytest
from banzai_nres.utils import array_utils


def test_finding_unique_elements_while_keeping_order():
    list = [1, 6, 7, 7, 2, 2]
    unique_list = array_utils.unique_elements_unordered(list)
    assert [1, 6, 7, 2] == unique_list