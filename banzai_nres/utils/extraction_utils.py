import numpy as np
from banzai_nres.utils.coordinate_utils import region_in_y_surrounding_trace
from banzai_nres.utils.array_utils import sum_like_index_coords
from banzai_nres.utils.trace_utils import get_trace_centroids_from_coefficients


class VerticalExtraction(object):
    """
    general class of methods for extracting spectrum assuming fixed x coordinate means
    fixed lambda.
    """
    def __init__(self):
        None

    def extract(self, flux_data, image, mask, weights=None):
        return vertical_x_y(flux_data, image, mask, weights=weights)

    def get_weights(self, *args):
        return None


class BoxExtraction(VerticalExtraction):
    def __init__(self):
        super(BoxExtraction, self).__init__()

    def get_weights(self, image, *args):
        return np.ones_like(image.data)


class OptimalFiberProfileExtraction(VerticalExtraction):
    def __init__(self):
        super(OptimalFiberProfileExtraction, self).__init__()

    def get_weights(self, image, *args):
        weights = image.fiber_profile.normalized_fiber_profile_image * image.ivar
        squared_profile = image.fiber_profile.normalized_fiber_profile_image **2 * image.ivar
        trace_centroid_y_values_per_order, num_traces, x = get_trace_centroids_from_coefficients(image.trace.coefficients, image)
        min_row, max_row = np.min(image.coordinates.y), np.max(image.coordinates.y)
        mask = np.ones_like(squared_profile, dtype=bool)
        for trace_number in range(num_traces):
            min_y, max_y = region_in_y_surrounding_trace(trace_centroid_y_values_per_order[trace_number], max_row, min_row)
            # the denominator of the weights involves box extracting the spectrum of the square profile
            numerator_of_weights, x_coords = extract_spectrum_single_order(squared_profile, image,
                                                                             image.coordinates.delta_y_from_trace, mask,
                                                                             image.coordinates.closest_trace_in_y,
                                                                             trace_number, window=12,
                                                                             weights=np.ones_like(image.data),
                                                                             ExtractMethod=BoxExtraction)
            min_col, max_col = np.min(x_coords), np.max(x_coords) + 1
            mask_of_vals_to_change = (image.coordinates.closest_trace_in_y[min_y:max_y, min_col: max_col] == trace_number)
            weights[min_y:max_y, min_col: max_col][mask_of_vals_to_change] /= \
                (numerator_of_weights * np.ones_like(weights[min_y:max_y, min_col: max_col]))[mask_of_vals_to_change]
        return weights


def vertical_x_y(flux_data, image, mask, weights=None):
    """
    :param image: banzai_nres image object with Coordinates() object appended as an attribute
    :param mask: boolean mask such that image.data[mask] returns only the flux values of interest.
    :param weights: The associated weights for each flux value to be used in extraction. e.g. if an array of 1's, then
                    this is box extraction.
    :return: one dimensional spectrum from vertical extraction. If weights are None, this is rectangular/box extraction
            if given optimal extraction weights, then this is optimal extraction
    """
    fluxes = flux_data[mask]
    weights = weights[mask]
    x_coords = image.coordinates.x[mask]
    one_d_spectrum, unique_x_coordinates = sum_like_index_coords(fluxes * weights, x_coords)
    return one_d_spectrum, unique_x_coordinates


def extract_spectrum_single_order(flux_data, image, delta_from_trace_coordinates, mask, ownership_coordinates, trace_to_extract, window=5, weights=None, ExtractMethod=VerticalExtraction):
    near_trace = (ownership_coordinates == trace_to_extract)
    mask = np.logical_and(mask, near_trace)
    mask[near_trace] = np.logical_and(mask[near_trace], np.abs(delta_from_trace_coordinates[near_trace]) < window)
    del near_trace
    one_d_spectrum, horizontal_coordinates = ExtractMethod().extract(flux_data, image, mask, weights=weights)
    return one_d_spectrum, horizontal_coordinates


def extract_spectrum_full_image(image, delta_from_trace_coordinates, ownership_coordinates, window=5, ExtractMethod=VerticalExtraction):
    """
    :param image:
    :param delta_from_trace_coordinates:
    :param ownership_coordinates:
    :param window:
    :param extract_method:
    :return: list of tuples of (flux_vs_horizontal_coordinate, horizontal_coordinates), where the elements of that tuple are
            an ndarray. E.g. for rectangular_x_y extraction, spectrum_for_each_order[4] = (spectrum_of_trace_4_vs_x, x_coords)
    """
    num_traces = int(np.max(ownership_coordinates) + 1)
    mask = np.ones_like(image.data, dtype=bool) # this would be a bpm.
    spectrum_for_each_order = []
    #start1 = time.clock()
    weights = ExtractMethod().get_weights(image)
    if weights is None:
        weights = np.ones_like(image.data)
    #print(time.clock() - start1, 'time to get weights')
    for i in range(num_traces):
        spectrum_for_each_order.append(extract_spectrum_single_order(image.data, image, delta_from_trace_coordinates,
                                       mask, ownership_coordinates, i, window=window, weights=weights,
                                       ExtractMethod=ExtractMethod))
    return spectrum_for_each_order