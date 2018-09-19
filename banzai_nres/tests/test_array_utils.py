import pytest
from banzai_nres.utils import array_utils


def test_finding_unique_elements_while_keeping_order():
    list_to_check = [1, 6, 7, 7, 2, 2]
    unique_list = array_utils.unique_elements_unordered(list_to_check)
    assert [1, 6, 7, 2] == unique_list