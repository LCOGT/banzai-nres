import numpy as np
from scipy import ndimage
from scipy import interpolate
from banzai_nres.utils.trace_utils import get_trace_centroids_from_coefficients


def region_in_y_surrounding_trace(trace_centroid_y_values, max_possible_row, min_possible_row):
    min_y, max_y = max(np.min(trace_centroid_y_values) - 40, min_possible_row), min(np.max(trace_centroid_y_values) + 40, max_possible_row)
    return int(min_y), int(max_y)


def generate_delta_y_and_closest_trace_coordinates(image):
    X, Y = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
    trace_y_values, num_traces, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
    delta_Y_from_trace_center = Y - trace_y_values[0] * np.ones_like(image.data)
    closest_trace = np.zeros(image.data.shape)
    max_row = np.max(Y) + 1
    min_row = np.min(Y)
    for trace_number in range(1, num_traces):
        min_trace_pos, max_trace_pos = region_in_y_surrounding_trace(trace_y_values[trace_number], max_row, min_row)
        Y_temp = Y[min_trace_pos : max_trace_pos]
        delta_Y_temp = Y_temp - trace_y_values[trace_number]
        indices_of_values_to_change = np.where(np.abs(delta_Y_from_trace_center[min_trace_pos: max_trace_pos]) > np.abs(delta_Y_temp))

        delta_Y_from_trace_center[min_trace_pos : max_trace_pos][indices_of_values_to_change] = delta_Y_temp[indices_of_values_to_change]
        closest_trace[min_trace_pos : max_trace_pos][indices_of_values_to_change] = trace_number

    image.coordinates.closest_trace_in_y = closest_trace
    image.coordinates.delta_y_from_trace = delta_Y_from_trace_center
