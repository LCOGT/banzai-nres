import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import interpolate
from utils.trace_utils import get_trace_centroids_from_coefficients
import time


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
        indices_of_values_to_change = np.where(np.abs(delta_Y_from_trace_center[min_trace_pos : max_trace_pos]) > np.abs(delta_Y_temp))

        delta_Y_from_trace_center[min_trace_pos : max_trace_pos][indices_of_values_to_change] = delta_Y_temp[indices_of_values_to_change]
        closest_trace[min_trace_pos : max_trace_pos][indices_of_values_to_change] = trace_number

    image.coordinates.closest_trace_in_y = closest_trace
    image.coordinates.x, image.coordinates.y = X, Y
    image.coordinates.delta_y_from_trace = delta_Y_from_trace_center


"""
deprecated functions for W-S transformations
"""

def generate_legendre_mth_derivative_coefficients(image, derivatives=1):
    """
    :param image: Banzai image object with trace_coefficients attached
    :param x: points across the CCD where you are evaluating the derivative.
    :return: The coefficients c=[0th, 1st, ] which when fed to the legendre package will return the
            legendre series of the mth derivative.
    """
    coefficients = image.trace.coefficients[:, 1:]
    num_of_traces = coefficients.shape[0]
    legendre_coefficients_list = []
    for trace in range(0, num_of_traces):
        legendre_coefficients_list.append(np.polynomial.legendre.legder(c=coefficients[trace], m=derivatives))
    return legendre_coefficients_list


def evaluate_legendre_from_coefficients(image, x, coefficients):
    """
    :param image: Banzai image object with trace_coefficients attached
    :param x: pixel points across the CCD where you are evaluating the derivative (not normalized).
    :param coefficients: A single list of coefficients which is passable to legendre.legval (e.g. c)
    :return: The legendre polynomial evaluated at the normalized points x
    """
    max_x = image.data.shape[1] - 1
    # x normalized to agree with -1 to 1 maximum range which the legendre fits are defined over.
    xnorm = convert_pixels_to_normed_coordinates(x, max_x)
    return np.polynomial.legendre.legval(xnorm, c=coefficients, tensor=True)


def find_next_guess_NR(guess, x, y, first_derivative_at_guess, second_derivative_at_guess, value_at_guess):
    """
    :param guess: guess for perpendicular intercept of the point with the trace
    :param x: x coordinate of the point
    :param y: y coordinates of the point
    :param first_derivative_at_guess:
    :param second_derivative_at_guess: second derivative of the trace at guess
    :param value_at_guess: value of the trace at guess
    :return: the next guess at the intercept of the traces perpendicular offshoot which intersects the point x,y.
    """
    f_prime = second_derivative_at_guess * (y - value_at_guess) - first_derivative_at_guess ** 2 - 1
    f = (x - guess) + first_derivative_at_guess * (y - value_at_guess)
    return guess - f / f_prime


def find_k_coordinates(image, slice_region, first_derivative_coefficients, second_derivative_coefficients, trace_coefficients_no_indices):
    """
    :param image:
    :param slice_region:
    :param derivative_coefficients: coefficients for the legendre poly series which describes the first derivative of the trace.
                                    As usual, defined over the range -1 to 1.
    :param trace_coefficients_no_indices: image.trace.coefficients[:, 1:] e.g. without the index column.
    :return: the array of values the same size as the image region, where each value is the k coordinate for the pixel.
            The k coordinate is the pixel value such that the arc length integral from midpoint to k on the trace, yields
            the w (trace length) coordinate of the pixel.
    """
    xnorm = np.linspace(-1, 1, image.data.shape[1])
    conversion_factor = xnorm[1] - xnorm[0] # conversion factor between dy/dx_norm and dy/dx (where x are pixels now)
    min_row, max_row, min_column, max_column = slice_region
    convergence_criterion = 1E-10  # in pixels

    x_values, y_values = np.meshgrid(np.arange(min_column, max_column + 1), np.arange(min_row, max_row + 1))
    k_coordinates = np.zeros_like(x_values, dtype=np.float64)
    trace_y_vals_for_slice = np.zeros_like(x_values, dtype=np.float64)

    converged = False
    counter = 0
    previous_guess = np.copy(x_values)
    while not converged:
        counter += 1
        first_derivatives = evaluate_legendre_from_coefficients(image, previous_guess, first_derivative_coefficients) * conversion_factor
        second_derivatives = evaluate_legendre_from_coefficients(image, previous_guess, second_derivative_coefficients) * conversion_factor ** 2
        trace_y_values = evaluate_legendre_from_coefficients(image, previous_guess, trace_coefficients_no_indices)

        new_guess = find_next_guess_NR(previous_guess, x_values, y_values, first_derivatives, second_derivatives, trace_y_values)
        if np.max(np.abs(new_guess - previous_guess)) < convergence_criterion:
            converged = True
            k_coordinates = new_guess
            trace_y_vals_for_slice = trace_y_values

        previous_guess = np.copy(new_guess)

    return k_coordinates, trace_y_vals_for_slice


def find_perpendicular_distances_for_slice(slice_region, k_coordinates, trace_y_vals_for_slice):
    min_row, max_row, min_column, max_column = slice_region
    x_values, y_values = np.meshgrid(np.arange(min_column, max_column + 1), np.arange(min_row, max_row + 1))
    perpendicular_distances = np.hypot(trace_y_vals_for_slice - y_values, k_coordinates - x_values)
    above_or_below_trace = np.ones_like(y_values) * (y_values > trace_y_vals_for_slice) - np.ones_like(y_values) * (y_values < trace_y_vals_for_slice)
    return perpendicular_distances * above_or_below_trace


def find_arc_length_distances_for_slice(image, k_coordinates, first_derivative_coefficients, a=-2, b=2):
    """
    Arc lengths along the trace at 10000 points are computed via simpsons rule. The cumulative integration is then
    spline interpolated and called at each value of k_coordinates.

    :param image:
    :param k_coordinates:
    :param first_derivative_coefficients:
    :param a: the lower bound of the bulk integration. typically -2 which corresponds to -1/2*image.data.shape[1] pixels.
                Should be plenty sufficient. CHANGE ONLY WITH CAUTION
    :param b: the lower bound of the bulk integration. typically -2 which corresponds to 3/2*image.data.shape[1] pixels.
                Should be plenty sufficient. CHANGE ONLY WITH CAUTION
    :return: an array of arc lengths where the i,j entry corresponds to the arc length integrated to the
                i,j point in k_coordinates

    """
    xnorm = np.linspace(-1, 1, image.data.shape[1])
    conversion_factor = xnorm[1] - xnorm[0] # conversion factor between dy/dx_norm and dy/dx (where x are pixels now)

    max_x = image.data.shape[1] - 1
    integration, normed_points_integrated_at = cumulative_simpson_integration(trace_arc_length_integrand,
                                                            a=a, b=b, num=10000,
                                                            extraargs=(first_derivative_coefficients, conversion_factor))
    # Converting the values of the cumulative arc length to pixels, away from a weird mixed set of units.
    real_pixels_integration = integration * max_x / 2

    # spline interpolating the integral for accurate recalling later.
    spline_of_arc_length = interpolate.UnivariateSpline(normed_points_integrated_at,
                                                        real_pixels_integration, s=0, ext=1)

    # converting to normed coordinates so that we can call the spline.
    normed_coordinates = convert_pixels_to_normed_coordinates(k_coordinates, max_x)
    arc_length_coordinates = spline_of_arc_length(normed_coordinates.flatten()).reshape(normed_coordinates.shape)

    return arc_length_coordinates


def trace_arc_length_integrand(x, derivative_coefficients, conversion_factor):
    return np.sqrt(1 + conversion_factor ** 2 * np.polynomial.legendre.legval(x, c=derivative_coefficients) ** 2)


def cumulative_simpson_integration(func, a, b, num, extraargs=()):
    """
    Cumulative simpson integration over 2*num intervals
    :param func: callable function of the integration variable, with possible extra arguments
    :param extraargs: extraargs for func.
    :param a: lower bound of integration
    :param b: upper bound
    :param num: the number of integration samples to return
    :return: the cumulative integration (type depends on function),
                the points the integration was evaluated at (type depends on function),

    the cumulative integration of the function from a to b over the intervals num. And the intervals
    NOTE: This is adapted from https://stackoverflow.com/questions/18215163/cumulative-simpson-integration-with-scipy
    """

    num *= 2
    a = float(a)
    b = float(b)
    h = (b-a)/num
    output = 4*func(a+h*np.arange(1, num, 2), *extraargs)
    tmp = func(a+h*np.arange(2, num-1, 2), *extraargs)
    output[1:] += tmp
    output[:-1] += tmp
    output[0] += func(a, *extraargs)
    output[-1] += func(b, *extraargs)
    points_integration_evaluated_at = a+h*np.arange(2, num+2, 2)
    return np.cumsum(output*h/3), points_integration_evaluated_at


def convert_pixels_to_normed_coordinates(pixels, max_x_pixel):
    """
    :param pixels: an array of pixels you wish to normalize such that 0,1,2.,,,,max_x_pixel falls between -1 and 1
    :param max_x_pixel: the maximum x pixel value of the pixels. This assumes np.min(pixels) = 0. I.e. indexing from zero
    :return: ndarray
    """
    return 2. * pixels / max_x_pixel - 1


def convert_normed_coordinates_to_pixels(normed_coordinates, max_x_pixel):
    """
    the inverse transform for convert_pixels_to_normed_coordinates
    """
    return (normed_coordinates + 1) * max_x_pixel / 2


def plot_normalized_flux_versus_perpendicular_coordinate(image, slice_region, perpendicular_distances, k_coordinates, trace_y_vals_for_slice):
    min_row, max_row, min_column, max_column = slice_region
    image_slice = image.data[min_row : max_row + 1, min_column : max_column + 1]
    var = np.copy(image_slice)
    var = var * (var > 100) + 0 * (var < 100) # simple variance estimate which includes the 10 count read noise.
    var += 100 # read noise estimate.

    fluxes_along_trace_centers = ndimage.map_coordinates(image.data.astype(np.float64),
                                                         [trace_y_vals_for_slice.flatten(), k_coordinates.flatten()], mode='constant', cval=0.0,
                                                         prefilter=True).reshape(image_slice.shape)

    normed_flux_values = (image_slice - np.min(image_slice))/(fluxes_along_trace_centers - np.min(image_slice))
    var = var/(fluxes_along_trace_centers ** 2)
    var[fluxes_along_trace_centers < 1E-3] = 0
    error = np.sqrt(var)
    # selecting only values which are less than 10 pixels away from center
    x = perpendicular_distances.flatten()[np.abs(perpendicular_distances.flatten()) < 10]
    y = normed_flux_values.flatten()[np.abs(perpendicular_distances.flatten()) < 10]
    errs = error.flatten()[np.abs(perpendicular_distances.flatten()) < 10]
    plt.errorbar(x, y, yerr=errs, xerr=1E-2, alpha=0.4, fmt="none")
    plt.ylim((-1/4,5/4))
    # xerr is the uncertainty in trace centroiding.
    #plt.plot(perpendicular_distances.flatten(), normed_flux_values.flatten(), '+')
    plt.title('Flux versus perpendicular distance from Trace')
    plt.show()