import numpy as np

from banzai_nres.utils.extract_utils import get_trace_xy_positions


def test_index_traces():
    traces = np.ones((3, 4))
    traces[[0, 0, 2, 0], [0, 1, 2, 3]] = 0
    idx, trace_idx = get_trace_xy_positions(traces)
    assert np.allclose(idx, [1])
    assert np.allclose(trace_idx[0], [[1, 1, 0, 1], [2, 2, 1, 2]])
    assert np.allclose(trace_idx[1], [[0, 1, 2, 3], [0, 1, 2, 3]])
