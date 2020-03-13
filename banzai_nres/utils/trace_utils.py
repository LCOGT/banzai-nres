import numpy as np


def get_trace_region(trace_mask: np.array):
    """
    Get the region of the image for a given trace

    :param trace_mask: boolean mask with pixels inside the trace of interest set to True
    :return: set of indexes that can be passed to any 2d array that is the same size as trace_mask, e.g. image.data

    :note: Using this return value e.g. image.data[get_trace_region(traces ==35)] will return a 2d array that can be
           written out directly to a fits file
    """
    pixel_indicies = np.where(trace_mask)

    # This ravels the data in the wrong direction, so we need to resort the indexes
    sorted_indices = np.lexsort((pixel_indicies[0], pixel_indicies[1]))

    # Assume that the traces are even bands (always the same number of pixels tall)
    nx = np.max(pixel_indicies[1]) - np.min(pixel_indicies[1]) + 1
    ny = len(sorted_indices) // nx
    return pixel_indicies[0][sorted_indices].reshape(nx, ny).T, pixel_indicies[1].reshape(nx, ny).T
