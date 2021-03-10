from banzai_nres.utils import trace_utils
import numpy as np


def test_get_region():
    traces = np.zeros((4, 4))
    traces[[2, 3], :] = 1
    trace_yx_pos = trace_utils.get_trace_region(traces == 1)
    assert np.allclose(traces[trace_yx_pos], 1)


def test_get_trace_region():
    nx, ny = 1001, 1003
    test_data = np.zeros((ny, nx))
    x0 = 498
    trace = np.round(2e-4 * (np.arange(nx) - x0) ** 2.0 + 0.01 * (np.arange(nx) - x0) + 502.0).astype(int)
    counter = 1
    for i in np.arange(nx, dtype=int):
        for j in range(-7, 8):
            test_data[trace[i] + j, i] = counter
            counter += 1

    traces = np.zeros((ny, nx), dtype=bool)
    traces[test_data > 0] = True
    trace_region = trace_utils.get_trace_region(traces)

    counter = 1
    expected_trace = np.zeros((15, nx), dtype=int)
    for i in np.arange(nx, dtype=int):
        for j in range(15):
            expected_trace[j, i] = counter
            counter += 1

    np.testing.assert_equal(test_data[trace_region], expected_trace)
