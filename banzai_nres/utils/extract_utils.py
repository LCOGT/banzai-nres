import numpy as np


def index_traces(traces):
    """
    :param traces: ndarray.
           A 2d image where every pixel within a single trace is identified by one index, e.g. 1
           and every pixel in the background is 0.
           For example:
           trace = np.array([[0, 0, 0],
                             [1, 1, 1],
                             [0, 0, 0],
                             [2, 2, 2]])
           represents an image with two traces running left to right, labelled trace 1 and trace 2
           respectively.
    :return: idx, trace_idxs: list, list
            idx is a list of integers, labelling each entry in trace_idxs
            trace_idxs is a list of ndarrays, call the ith one A_i. Each A_i is such that
            data[A_i] will yield a spectrum of shape N_i, M where M is traces.shape[1]
            and N_i is the width of the trace in question.
    """
    if np.allclose(traces, 0):
        # if there are no traces, return
        return [], []
    nx = traces.shape[1]
    # TODO the below is a lazy and very slow version of turning the traces array into a set
    #  of arrays A_i such that image.data[A_i] gives a properly flattened spectrum for the trace.
    #  I or someone else will make something much faster later.
    #  Also it does not work right now.
    idx = []
    trace_idxs = []
    for i in range(1, int(np.max(traces)) + 1):
        idx.append(i)
        row, col = np.where(np.isclose(traces, i))
        row, col = row.reshape(-1, nx), col.reshape(-1, nx)
        trace_idxs.append([row, col])

    return idx, trace_idxs
