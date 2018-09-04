import numpy as np


def find_nearest(value, array):
    """
    :param value: some scalar element to search for inside of array.
    :param array: array to search for elements near in value to element
    :return : array element which is the nearest to value.
    """
    nearest_array_value = array[np.abs(array - value) == np.min(np.abs(array - value))][0]
    return nearest_array_value


def sum_like_index_coords(elements_to_sum, index_coords):
    """
    :param elements_to_sum: ndarray or list length N, e.g. flux per x and y coordinates [flux_at1,2, flux_at1,5,...]
    :param index_coords: ndarray or list length N, e.g. x coordinates [1,1,1,1,2,3,4,4,4,...]
    :return: the sum of elements which have like index coordinates (ndarray) and the unique index coordinates.
                e.g. extracted flux at each x-pixel.
    This is substantially faster than using the simpler for loop version (about a factor of 10). And faster than
    boolean indexing version by a factor of 3 (I don't know why...)
    """
    index_coords = np.array(index_coords)
    sorted_elements_to_sum = elements_to_sum[index_coords.argsort()]
    unique_index_coords, counts = np.unique(index_coords, return_counts=True)
    positions = np.zeros(len(counts) + 1, dtype=np.int)
    positions[1:] = np.cumsum(counts).astype(np.int)
    sums = [sum(sorted_elements_to_sum[positions[i - 1]:positions[i]]) for i in range(1, len(positions))]
    return np.array(sums), unique_index_coords


def evaluate_functions_only_on_equal_index_coords(elements_to_eval, index_coords, func, function_coefficients_list):
    """
    :param elements_to_eval:
    :param index_coords:
    :param func: function of the type func(x, *params)
    :param function_coefficients_list: array such that the ith params = function_coefficients_list[:, i]
    :return: similar to sum_like_index_coords above, except we act a different function on the elements instead of summing
    them together. Also returns a tuple which is suitable for input into any array which gives the indices of the values
    to change to outputs.
    """
    index_coords = np.array(index_coords)
    argsort_indices = index_coords.argsort()
    sorted_elements_to_eval = elements_to_eval[argsort_indices]
    unique_index_coords, counts = np.unique(index_coords, return_counts=True)
    positions = np.zeros(len(counts) + 1, dtype=np.int)
    positions[1:] = np.cumsum(counts).astype(np.int)
    outputs = []
    for i in range(1, len(positions)):
        outputs.extend(list(func(sorted_elements_to_eval[positions[i - 1]:positions[i]], *function_coefficients_list[:, i-1])))
    return np.array(outputs), unique_index_coords, argsort_indices


def fill_with_nearest_left_value_if_flagged_as_false(boolean_array, list_to_patch):
    """
    :param boolean_array: ndarray. e.g. [False, True, True]
    :param list_to_patch: list. e.g. [None, real_value1, real_value2]
    :return: [real_value1, real_value1, real_value2]
    """
    assert not (boolean_array == False).all()
    reference_array = boolean_array.astype(np.int) * np.arange(1, len(boolean_array) + 1)
    reference_array[reference_array == 0] = -1000.0
    reference_array = reference_array.astype(np.float64)
    reference_array -= 0.9 # so we select values to the left of the missing value
    index_array = np.arange(len(boolean_array))
    for flagged_index in index_array[~boolean_array]:
        list_to_patch[flagged_index] = list_to_patch[int(find_nearest(flagged_index, reference_array))]
