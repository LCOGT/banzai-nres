import numpy as np


def index_traces(traces):
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
    # TODO make this faster.
    for trace_id in trace_ids:
        y, x = np.where(np.isclose(traces, trace_id))
        sorted_indices = np.lexsort((y, x), axis=0)
        trace_width = int(max(x) - min(x)) + 1
        xpos = x[sorted_indices].reshape(-1, trace_width)
        ypos = y[sorted_indices].reshape(-1, trace_width)

        #xpos, ypos = x.reshape(-1, trace_width), y.reshape(-1, trace_width)

        trace_xypos.append([xpos, ypos])

    return trace_ids, trace_xypos
