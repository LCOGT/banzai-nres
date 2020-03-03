import numpy as np


def get_region(traces):
    # TODO REPLACE WITH CURTIS' get_region FUNCTION
    """
    :param traces: ndarray.
           A 2d image where every pixel within a single trace is identified by a 1
           and every pixel in the background (and regions belonging to other traces) is 0.
           For example:
           trace = np.array([[0, 0, 0],
                             [1, 1, 1],
                             [0, 0, 0],
                             [0, 0, 0]])
           represents an image one trace running left to right
    :return: trace_yxpos: list of indices
             trace_yxpos is a list of coordinates suitable for numpy's fancy indexing. Specifically,
             trace_yxpos is [ypos, xpos], such that data[ypos, xpos] will yield a spectrum of shape.
             E.g. xpos.shape = N, M where M is the width (x_extent) of the trace
             and N is the vertical width of the trace. e.g. N, M = 15, 4096 pixels.
             ypos.shape == xpos.shape
    """
    if np.allclose(traces, 0):
        return []
    ny, nx = traces.shape
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    y, x = np.where(np.isclose(traces, 1))
    minx = np.min(x)
    trace_width = int(np.max(x) - minx) + 1
    trace_height = int(x.size / trace_width)
    xpos = (np.arange(np.min(x), np.max(x) + 1) * np.ones((trace_height, trace_width))).astype(int)
    ypos = np.zeros_like(xpos)
    for j in np.arange(trace_width):
        ypos[:, j] = y2d[:, minx + j][np.isclose(traces[:, minx + j], 1)]

    return [ypos, xpos]
