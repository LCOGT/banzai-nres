import numpy as np
from scipy.ndimage import map_coordinates


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
    # should test this function in two ways, one should verify that an already flattened trace is not modified.
    rectified_order = np.zeros((2*half_window + 1, image_data.shape[1]))
    x_coords = np.arange(image_data.shape[1])
    for offset, row in zip(np.arange(-half_window, half_window + 1), np.arange(rectified_order.shape[0])):
        mapped_y_values = map_coordinates(image_coordinates['y'], [single_order_centers + offset, x_coords],
                                          order=0, mode='constant', cval=0, prefilter=False)
        rectified_order[row] = image_data[(mapped_y_values, x_coords)]
        if nullify_mapped_values:
            image_data[(mapped_y_values, x_coords)] = 0
    return rectified_order, image_data
