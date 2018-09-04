"""
Sets of functions for extracting the fiber intensity profiles in the 'x' direction (not lambda). 'x' could
be either the physical x coordinate, or the arc length coordinate, or any choice.
"""
from scipy import ndimage, interpolate
from scipy import optimize, integrate

import numpy as np
from .array_utils import evaluate_functions_only_on_equal_index_coords
from .coordinate_utils import  region_in_y_surrounding_trace
import matplotlib.pyplot as plt
import time


class Shapelets(object):
    def __init__(self):
        None

    def normalize_profile(self, fit_coefficients):
        """
        :param fit_coefficients:
        :return: coefficients which yield a profile which integrates to one between -8 and 8 full width half maxes.
        """
        fit_coefficients[0] = 0  # removing offset so that pdf has finite area.
        fwhm = fit_coefficients[1]
        A, err = integrate.quad(self.profile_function, a=-8*fwhm, b=8*fwhm, args=(*fit_coefficients,))
        fit_coefficients[3:] /= A # normalizing the function.
        return fit_coefficients

    def full_width_half_max(self, fit_coefficients):
        """
        :param fit_coefficients:
        :return: a measure of the full width half max of the fiber profile.
        """
        return fit_coefficients[1]

    def profile_function(self, x, *params):
        """
        :param x: array which you wish to evaluate the sum of the shapelets along
        :param params: (linear_offset, fwhm, center, a0, a1, a2,...)
        :return: y values of sum, ndarray type - e.g. the value of the shapelet profile function at each x value.
        Note that fhwm here is not truly the full with half max.
        Note we require a linear offset because the shapelet series does not contain a 0th order term (as each shapelet is
        multiplied by np.exp(-x ** 2 / 2)
        """
        (offset, fwhm, center), shapelet_amplitudes = params[:3], params[3:]
        x = x * 1.
        x -= center
        x /= fwhm
        num_of_shapelets = int(len(shapelet_amplitudes))
        y = np.ones_like(x, dtype=np.float64) * offset
        for i in range(num_of_shapelets):
            y += shapelet_amplitudes[i] * nth_shapelet(x, i)
        return y

    def initial_guess(self, size_of_basis):
        return (0, 2, 0) + tuple([1 / size_of_basis] * size_of_basis)


def hermite_n(x, n):
    coefficients = [0]*n + [1]
    return np.polynomial.hermite.hermval(x, c=np.array(coefficients))


def nth_shapelet(x, n):
    # omitted n factorial component inside of the squareroot.
    return np.sqrt(2 ** n * np.sqrt(np.pi)) * hermite_n(x, n) * np.exp(-x ** 2 / 2)


def fit_fiber_intensity_function_over_window(filtered_flux_data, image, window, horiz_coords, max_distance_from_trace_center, delta_from_trace_coords, trace_ownership_coords, trace_number, trace_centroid_values, coeffs_guess, model_intensity_function, size_of_basis = 10):
    """
    :param horiz_coords: grid of horizontal coordinates, e.g. image.coordinates.x
    :param vertical_coords: grid of vertical coordinates, e.g. image.coordinates.y or image.coordinates.perpendiculars
    :param delta_from_trace_coords: grid of delta_y to the nearest trace
    :param flux_data: image.data
    :param trace_ownership_coords: grid of integer values where trace_ownership_coords[x,y] = index of the trace which this x,y point is closest to.
    :param trace_number: the index of the trace which
    :param trace_centroid_values: A grid of values like horiz_coords which has the y_values of the trace running across it
                                evaluated at each of the horiz_coords.
    :param window: min_column, max_column where you want to evaluate the fit (e.g.
    :param coeffs_guess: the initial guess for the fitting procedure.
    :param model_intensity_function: the function which models the fiber intensity profile along equipotential lines of vertical_coords
    :param size_of_basis: the number of basis functions to use in the fit. For shapelets, 10 is good.
    :return: the optimum fit coefficients for model_intensity_function which describe the line profile.
    """
    min_x, max_x, min_y, max_y = window
    belong_to_trace = (trace_ownership_coords == trace_number)
    closer_to_trace = (np.abs(delta_from_trace_coords) < max_distance_from_trace_center)
    # this could be improved further by delta_from_trace_coords[belong_to_trace] < max_
    # which would yield a smaller mask that one could obtain by arr_to_mask[belong_to_trace] then arr_to_mask[closer_to_trace]
    mask_of_values_to_extract = np.logical_and(belong_to_trace, closer_to_trace)

    max_fluxes = ndimage.map_coordinates(filtered_flux_data,
                                         [trace_centroid_values[mask_of_values_to_extract],
                                          horiz_coords[mask_of_values_to_extract]], mode='constant', cval=0.0,
                                         prefilter=False)
    fluxes = image.data.astype(np.float64)[min_y : max_y, min_x: max_x][mask_of_values_to_extract]
    vertical_values_to_fit = delta_from_trace_coords[mask_of_values_to_extract]
    ivar = image.ivar[min_y : max_y, min_x: max_x][mask_of_values_to_extract]
    var = np.reciprocal(ivar)
    fluxes /= max_fluxes
    var /= max_fluxes ** 2

    sigmas = np.sqrt(var)
    """
    plt.figure()
    plt.errorbar(vertical_values_to_fit, fluxes, yerr=sigmas, xerr=1E-2, fmt="none")
    plt.show()
    """
    # fitting with flux errors included!
    try:
        # currently 75% of the computation time is spent fitting.
        fit_coeffs, pcov = optimize.curve_fit(model_intensity_function, vertical_values_to_fit, fluxes, p0=coeffs_guess,
                                              sigma=sigmas, method='lm')
    except ValueError or RuntimeError:
        print('failed on trace number %s over x,x,y,y window %s, '
              'likely the lampflat signal is too dim to resolve the line profile.'%(trace_number, window))
        fit_coeffs = coeffs_guess

    """
    plt.figure()
    continuous_vals = np.linspace(-10, 10, 200)
    plt.plot(continuous_vals, model_intensity_function(continuous_vals, *fit_coeffs))
    plt.errorbar(vertical_values_to_fit, fluxes, yerr=sigmas, xerr=1E-2, fmt="none")
    plt.show()
    """
    return fit_coeffs


def fit_vertical_fiber_intensity_functions_over_given_horizontal_ranges(filtered_image_data, image, trace_centroid_values_per_order, trace_number, horizontal_ranges, Model=Shapelets, size_of_basis=10):
    """
    horizontal_windows: list of the form [(low, high), (low2, high2),...]
    """
    trace_centroid_values = trace_centroid_values_per_order[trace_number] * np.ones_like(image.data)
    horiz_coords = image.coordinates.x
    vertical_coords = image.coordinates.y
    delta_from_trace_coords = image.coordinates.delta_y_from_trace
    trace_ownership_coords = image.coordinates.closest_trace_in_y
    model = Model()
    coeffs_guess = model.initial_guess(size_of_basis)
    model_fit_coeffs = []
    max_possible_row, min_possible_row = np.max(vertical_coords), np.min(vertical_coords)
    for window in horizontal_ranges:
        min_y, max_y = region_in_y_surrounding_trace(trace_centroid_values_per_order[trace_number], max_possible_row, min_possible_row)
        window += (min_y, max_y)
        min_x, max_x, min_y, max_y = window
        max_distance_from_trace_center = 10
        fiber_fit_coeffs = fit_fiber_intensity_function_over_window(filtered_image_data, image, window, horiz_coords[min_y:max_y, min_x: max_x]
                                                                    , max_distance_from_trace_center, delta_from_trace_coords[min_y:max_y, min_x: max_x],
                                                                    trace_ownership_coords[min_y:max_y, min_x: max_x], trace_number,
                                                                    trace_centroid_values[min_y:max_y, min_x: max_x], coeffs_guess,
                                                                    model.profile_function, size_of_basis=size_of_basis)
        coeffs_guess = 1. * fiber_fit_coeffs
        model_fit_coeffs.append(fiber_fit_coeffs)

    return model_fit_coeffs


def normalize_fiber_fits(model_fit_coeffs, Model=Shapelets):
    model = Model()
    normalized_fiber_profile_coefficients = []
    for individual_fit_coeffs in model_fit_coeffs:
        individual_fit_coeffs = model.normalize_profile(individual_fit_coeffs)
        normalized_fiber_profile_coefficients.append(individual_fit_coeffs)
    return normalized_fiber_profile_coefficients


def interpolate_fiber_fits(normalized_coeffs, sample_points):
    """
    This builds linear interpolation objects for each for the shapelet coefficients across the detector.
    This uses coefficient values at boundary for any calls outside of the interpolated range.
    :param sample_points = ndarray 1 dimensional, corresponding to the coordinate points
                           where the normalized_coeffs were fit at.
    """
    normalized_coeffs = np.array(normalized_coeffs)
    coefficient_generating_functions = []
    for i in range(normalized_coeffs.shape[1]):
        coefficient_generating_functions.append(interpolate.UnivariateSpline(sample_points,
                                   normalized_coeffs[:, i], k=1, s=0, ext='const'))

    return coefficient_generating_functions


def sample_coefficients_from_generating_functions(coefficient_generating_functions, sample_points, Model=Shapelets, renormalize=False):
    """
    :param coefficient_generating_functions:
    :param sample_points:
    :return: [0thcoefficient at sample_point1, 0thcoefficient at sample_point2,..., nthcoefficient at sample_point1,...]
            normalized if normalize=True
    renormalization is often not necessary if you sampled the coefficients over a large enough region (e.g. 10 points across the
    detector). Typical errors using 10 points gives only 1E-5 absolute error, e.g. your peak value will be 1E-5 away from the truly normalized value.
    """
    num_of_coefficients = len(coefficient_generating_functions)
    sampled_coefficients = np.array([coefficient_generating_functions[i](sample_points) for i in range(num_of_coefficients)])
    if renormalize:
        model = Model()
        sampled_coefficients = np.array([model.normalize_profile(sampled_coefficients[:, i]) for i in range(len(sample_points))]).T
    return sampled_coefficients


def evaluate_normalized_fiber_profile_across_detector_x_y_per_trace(image, coefficient_generating_functions, trace_number, window, Model=Shapelets, renormalize=False):
    """
    :param image:
    :param coefficient_generating_functions:
    :param trace_number:
    :param window: max vertical extent we evaluate the fiber profile through, in pixels.
    :param Model:
    :return:
    """
    belong_to_trace = image.coordinates.closest_trace_in_y == trace_number
    delta_y_trace = image.coordinates.delta_y_from_trace[image.coordinates.closest_trace_in_y == trace_number]
    near_trace = delta_y_trace < window
    delta_y_trace = delta_y_trace[near_trace]
    # the above could be sped up by limiting closest_trace_in_y[min_trace_pos, max_trace_pos], but one would evaluate
    # that outside of the loop via trace positions.
    x_coords = image.coordinates.x[belong_to_trace][near_trace]
    y_coords = image.coordinates.y[belong_to_trace][near_trace]
    interpolated_coefficients = sample_coefficients_from_generating_functions(coefficient_generating_functions,
                                                                              np.unique(x_coords), Model=Model,
                                                                              renormalize=renormalize)
    normalized_fiber_profile_values, unique_x_coords, argsort_indices = evaluate_functions_only_on_equal_index_coords(delta_y_trace, x_coords,
                                                                                             Model().profile_function,
                                                                                             interpolated_coefficients)
    indices = (y_coords[argsort_indices], x_coords[argsort_indices])

    return normalized_fiber_profile_values, indices