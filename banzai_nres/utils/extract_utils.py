import numpy as np


def get_trace_xy_positions(traces):
    """
    :param traces: ndarray.
           A 2d image where every pixel within a single trace is identified by one index, e.g. 1,
           and every pixel in the background is 0.
           For example:
           trace = np.array([[0, 0, 0],
                             [1, 1, 1],
                             [0, 0, 0],
                             [2, 2, 2]])
           represents an image with two traces running left to right, labelled trace 1 and trace 2
           respectively.
    :return: trace_ids, trace_xypos: list, list
            trace_ids = [1,2,....,N] where N is the total number of traces in the image.
            trace_xypos is a list of ndarrays, call the ith one A_i. Each A_i is such that
            data[A_i] will yield a spectrum of shape N_i, M where M is traces.shape[1]
            and N_i is the width of the trace in question.
    """
    if np.allclose(traces, 0):
        # if there are no traces, return empty lists.
        return [], []
    trace_ids = np.arange(1, int(np.max(traces)) + 1)
    trace_xypos = []
    # TODO make this faster. This is a placeholder for now.
    ny, nx = traces.shape
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))
    for trace_id in trace_ids:
        y, x = np.where(np.isclose(traces, trace_id))
        minx = np.min(x)
        trace_width = int(np.max(x) - minx) + 1
        trace_height = int(x.size / trace_width)
        xpos = (np.arange(np.min(x), np.max(x) + 1) * np.ones((trace_height, trace_width))).astype(int)
        ypos = np.zeros_like(xpos)
        for j in np.arange(trace_width):
            ypos[:, j] = y2d[:, minx + j][np.isclose(traces[:, minx + j], trace_id)]

        trace_xypos.append([ypos, xpos])

    return trace_ids, trace_xypos
