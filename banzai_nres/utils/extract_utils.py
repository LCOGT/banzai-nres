import numpy as np
import abc
from scipy.ndimage import map_coordinates


class Extract(object):
    @abc.abstractmethod
    def extract(self):
        pass

    @staticmethod
    def extract_order(twod_spectrum, weights=None):
        if weights is None:
            return np.sum(twod_spectrum, axis=0)
        else:
            return np.sum(twod_spectrum * weights, axis=0)


def rectify_orders(image_data, trace, half_window=10, debug=False):
    rectified_orders = {}
    num_orders = trace.num_traces_found()
    image_copy = np.copy(image_data)
    x_coordinates, y_coordinates = np.meshgrid(np.arange(image_data.shape[1]), np.arange(image_data.shape[0]))
    image_coordinates = {'x': x_coordinates, 'y': y_coordinates}
    for index in np.arange(num_orders):
        rectified_orders[trace.get_id(index)], image_copy = rectify_order(image_copy, image_coordinates,
                                                                          trace.get_centers(index),
                                                                          half_window=half_window,
                                                                          nullify_mapped_values=True)
    if debug:
        return rectified_orders, image_copy
    else:
        return rectified_orders


def rectify_order(image_data, image_coordinates, single_order_centers, half_window=10, nullify_mapped_values=True):
    """
    :param image_data:
    :param image_coordinates: dictionary with x and y coordinates for each pixel.
    :param single_order_centers: the y centers of the center of the rectification window
    :param half_window: the half width of the window about which one is rectifying the trace/order
    :param nullify_mapped_values: Default true, this sets to zero the values which have been mapped to the
                                  rectified grid.
    :return: An array of 2*half_window + 1 rows by width(image_data) about the order. The center of the order
             is right at the index of half_window, e.g. rectified_order[half_window] gives the flux down the
             center of the order.
    """
    # should test this function in two ways, one should verify that an already flattened trace is not modified.
    rectified_order = np.zeros((2*half_window + 1, image_data.shape[1]))
    x_coords = np.arange(image_data.shape[1])
    for offset, row in zip(np.arange(-half_window, half_window + 1), np.arange(rectified_order.shape[0])):
        mapped_y_values = map_coordinates(image_coordinates['y'], [single_order_centers + offset, x_coords],
                                          order=0, mode='constant', cval=0, prefilter=False)
        rectified_order[row] = image_data[(mapped_y_values, x_coords)]
        if nullify_mapped_values:
            image_data[(mapped_y_values, x_coords)] = 0
    # TODO we will need to adopt the true x, y positions relative to the trace center from this as well.
    return rectified_order, image_data