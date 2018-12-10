#! /usr/bin/python3
"""
trace_utils.py: Routines for finding echelle orders across a CCD.
Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)

    Tim Brandt (tbrandt@physics.ucsb.edu)

Every function here is covered by either a unit test or an integration test inside of tests/test_trace_maker.py
"""

import numpy as np
from scipy import ndimage, optimize
from astropy.table import Table
import copy

from banzai_nres.utils.array_utils import unique_elements_unordered
import logging

logger = logging.getLogger(__name__)


class Trace(object):
    # TODO: add generate trace y values as a method of this class.- including image width and all necessary attributes
    # to do so.
    """
    Object for storing all the trace related attributes. This gets appended to each Image instance.
    """
    def __init__(self):
        self.coefficients = None
        self.fiber_order = None
        self.trace_center_table_name = 'trace_center'
        self.coefficients_table_name = 'coefficients'

    def get_trace_centroids_from_coefficients(self, image_width, coefficients_and_indices=None):
        """
        :param coefficients_and_indices: polynomial fit coefficients which describe the traces. Legendre polynomials
                                        normalized between -1 and 1.
        :param image_width: image.data.shape[1]
        :return: trace centroids for each trace, versus x pixel. E.g. trace_values_versus_xpixel[2,5] is the 3rd orders value
                at x=5.
                num_traces = num_orders*Num_fibers
                x = [0,1,2,...,image.data.shape[1]-1]
        """
        if coefficients_and_indices is None:
            coefficients_and_indices = self.coefficients
        coeflen, coefwidth = coefficients_and_indices.shape
        num_traces, order_of_poly_fits = coeflen, coefwidth - 2
        legendre_polynomial_array, x, xnorm = generate_legendre_array(image_width, order_of_poly_fits)
        trace_values_versus_xpixel = np.dot(coefficients_and_indices[:, 1:], legendre_polynomial_array)
        return trace_values_versus_xpixel, num_traces, x

    def construct_undesignated_fiber_order(self, num_lit_fibers):
        alphabetical_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        fiber_order = tuple(alphabetical_list[:num_lit_fibers])
        return fiber_order

    def convert_numpy_array_coefficients_to_astropy_table(self, num_lit_fibers, fiber_order=None, coefficients=None):
        if coefficients is None:
            coefficients = self.coefficients
        order_numbers = list(coefficients[:, 0].astype(np.int))
        coefficients_table = self.generate_astropy_table_from_numpy_array_and_orders(self.coefficients_table_name,
                                                                                     num_lit_fibers,
                                                                                     order_numbers,
                                                                                     coefficients[:, 1:], fiber_order)
        coefficients_table[self.coefficients_table_name].description = 'Legendre polynomials ' \
                                                                                   'normalized from -1 to 1 over the ' \
                                                                                   'number of columns of the image, ' \
                                                                                   'usually 0 to 4095'
        return coefficients_table

    def convert_numpy_array_trace_centroids_to_astropy_table(self, num_lit_fibers, trace_centroids, coefficients, fiber_order=None):
        order_numbers = list(coefficients[:, 0].astype(np.int))
        trace_centroids_table = self.generate_astropy_table_from_numpy_array_and_orders(self.trace_center_table_name,
                                                                                        num_lit_fibers,
                                                                                        order_numbers,
                                                                                        trace_centroids, fiber_order)
        trace_centroids_table[self.trace_center_table_name].unit = 'pixel'
        trace_centroids_table[self.trace_center_table_name].description = 'y pixel position for trace' \
                                                                                         'center at each x pixel' \
                                                                                         'from 0 to number of columns' \
                                                                                         'in the image'
        return trace_centroids_table

    def convert_astropy_table_coefficients_to_numpy_array(self, astropy_table_of_coefficients):
        coefficients_and_indices, fiber_order = self.recombine_values_from_table_into_nd_array_with_order_indices(
                                                                                     astropy_table_of_coefficients,
                                                                                     self.coefficients_table_name)
        return coefficients_and_indices, fiber_order

    def convert_astropy_table_trace_y_values_to_numpy_array(self, astropy_table_of_trace_centroids):
        trace_values_with_indices, fiber_order = self.recombine_values_from_table_into_nd_array_with_order_indices(
                                                                                       astropy_table_of_trace_centroids,
                                                                                       self.trace_center_table_name)
        trace_values_versus_xpixel = trace_values_with_indices[:, 1:]
        return trace_values_versus_xpixel, fiber_order

    def generate_astropy_table_from_numpy_array_and_orders(self, array_name, num_lit_fibers, order_numbers_list, array_of_values, fiber_order):
        if fiber_order is None:
            fiber_order = self.fiber_order
        if fiber_order is None:
            fiber_order = self.construct_undesignated_fiber_order(num_lit_fibers)

        total_num_of_orders = np.max(np.array(order_numbers_list)) + 1
        fiber_designations = []
        list_of_array_values_per_trace = []
        for fiber_name in fiber_order:
            fiber_designations.extend([fiber_name] * total_num_of_orders)
        array_of_array_values_per_trace = np.split(array_of_values, array_of_values.shape[0], axis=0)
        for two_d_array in array_of_array_values_per_trace:
            list_of_array_values_per_trace.append(list(two_d_array[0]))

        names = ('fiber', 'order', array_name)
        coefficients_table = Table([fiber_designations, order_numbers_list,
                                    list_of_array_values_per_trace],
                                   names=names)
        return coefficients_table

    def recombine_values_from_table_into_nd_array_with_order_indices(self, astropy_table, name_of_values):
        column_names = astropy_table.colnames
        order_numbers = np.array([np.array(astropy_table[column_names[1]])]).T
        list_of_array_values_per_trace = astropy_table[name_of_values]
        values_and_indices = np.hstack((order_numbers, np.vstack(list_of_array_values_per_trace)))
        fiber_order = tuple(unique_elements_unordered(list(astropy_table[column_names[0]])))
        if fiber_order == self.construct_undesignated_fiber_order(num_lit_fibers=len(fiber_order)):
            fiber_order = None
        return values_and_indices, fiber_order


def maxima(A, s, k, ref):
    # TODO: replace this with a function from a package. This works fine, but for maintainability reasons we
    # want a function from a package, like scipy signal peak finder.
    """
    A procedure for finding maxima which are not just peaks in the noise floor. E.g. searches for statistically
    significant maxima.
    :param A: An array of values where you wish to find maxima
    :param s: The window were the central-point's value must be larger than all the other values in the window to qualify
              as a maximum.
    :param k: multiplicative constant which qualifies a point as a real maxima (not just noise)
    :param ref: the reference value which k*r defines the value any real maximal point must exceed
    :return: The index of the a point near a maximum  and a boolean for whether or not
            a valid maximum exists.
    """

    i = s - 1
    l, r = 0, 0
    A = A - np.ones_like(A) * np.min(A)
    threshold = k * ref
    first_max_index = 0
    maximum_exists = False
    while i < len(A) - s and int(l + r) != 2:
        l, r = 1, 1
        for j in range(1, s + 1):
            if A[i] < A[i + j] or A[i] < threshold:
                r, j = 0, s + 1
            if A[i] < A[i - j] or A[i] < threshold:
                l, j = 0, s + 1
        if int(l + r) == 2:
            first_max_index = i
            maximum_exists = True
        i += 1
    return first_max_index, maximum_exists


def negative_flux_across_trace(legendre_polynomial_coefficients, imfilt, x, evaluated_legendre_polynomials):
    """
    :param legendre_polynomial_coefficients: list of legendre polynomial coefficients
    :param imfilt: ndarray, image_data passed through ndimage.spline_filter if you are using this with scipy optimize.
            If you just want the total flux, just pass image_data.
    :param evaluated_legendre_polynomials: Legendre polynomial evaluated over -1 to 1 , ndarray.
    :return: The negative of the total flux summed across the trace.
    NOTE: This requires that the input filtered_image_data has been spline filtered by ndimage.spline_filter
    """
    y = x * 0.

    for i in range(len(legendre_polynomial_coefficients)):
        y += evaluated_legendre_polynomials[i] * legendre_polynomial_coefficients[i]
    val = np.sum(ndimage.map_coordinates(imfilt, [y, x], prefilter=False))
    return -1 * val


def flux_across_trace_up_detector(testpoints, legendre_polynomial_coefficients, imfilt, x, evaluated_legendre_polynomials):
    """
    :param legendre_polynomial_coefficients: list of legendre polynomial coefficients excluding 0th order coefficients
    :param imfilt: ndarray. image_data passed through ndimage.spline_filter if you are using this with scipy optimize.
            If you just want the total flux, just pass image_data.
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param evaluated_legendre_polynomials: Legendre polynomial evaluated over -1 to 1 , ndarray.
    :return: The positive total flux at each point in testpoints
    """
    values = np.zeros_like(testpoints)
    for i in range(0, len(values)):
        coeffguesses = list([testpoints[i]] + legendre_polynomial_coefficients)
        values[i] = (-1) * negative_flux_across_trace(coeffguesses, imfilt, x, evaluated_legendre_polynomials)
    return values


def generate_initial_guess_for_trace_polynomial(image_data, x, evaluated_legendre_polynomials, order=2, second_order_coefficient_guess=90, lastcoef=None, direction='up'):
    """
    :param image_data: Ndarray, image. E.g. image.data for a banzai image object.
    :param x:
    :param evaluated_legendre_polynomials:
    :param order:
    :param second_order_coefficient_guess:
    :param lastcoef:
    :param direction:
    :return: guess for the coefficient of the polynomial fit for the following trace.

    If the last set of coefficients are None, then we use a manual guess for the first order to fit set by the extent
    of a trace across the detector.

    If a previous fit exists, we find a new zero-th order term for the polynomial fit by
    stepping along either up or down the grid until we find
    a point which is a maximum and nearest to the last order.
    bnds = [[None,None]]*(order+1)
    """

    if lastcoef is None:
        p = [0. for i in range(order + 1)]
        p[0] = int(image_data.shape[0]/3)
        if order >= 2:
            p[2] = second_order_coefficient_guess
        coeffsguess = copy.deepcopy(p)
        maximum_exists = True
        refflux = 0
    else:

        p = list(lastcoef[1:])
        p0 = int(lastcoef[0])

        if direction == 'up':
            testpoints = list(range(p0 + 6, p0 + 100))
        elif direction == 'down':
            testpoints = list(range(p0 - 100, p0 - 6))
            testpoints = testpoints[::-1]
        elif direction == 'inplace':
            testpoints = list(range(p0 - 10, p0 + 10))

        fluxvals = flux_across_trace_up_detector(testpoints, p, image_data, x, evaluated_legendre_polynomials)
        refflux = max((-1) * negative_flux_across_trace(lastcoef, image_data, x, evaluated_legendre_polynomials), max(fluxvals))

        deltap0guess, maximum_exists = maxima(fluxvals, 5, 1 / 20, refflux)

        if direction == 'up':
            p0 = deltap0guess + min(testpoints)
        elif direction == 'down':
            p0 = max(testpoints) - deltap0guess
        elif direction == 'inplace':
            p0 = deltap0guess + min(testpoints)

        coeffsguess = [p0] + list(p)
    return coeffsguess, maximum_exists, refflux


def find_order(image_data, imfilt, x, evaluated_legendre_polynomials, order=2, second_order_coefficient_guess=90, lastcoef=None, direction='up'):
    """
    :param image_data: Ndarray. Image_data e.g. image.data if image is a banzai image object.
    :param imfilt: ndarray, image.data passed through ndimage.spline_filter
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param evaluated_legendre_polynomials: Legendre polynomial evaluated over -1 to 1 , ndarray.
    :param order: order of the legendre polynomial fit
    :param second_order_coefficient_guess: coefficient guess for the second order legendre polynomial which
    describe the traces across the ccd. The 70-100 works well for all LCOGT NRES instruments as of 2018/09/07
    :param lastcoef: [0th_order_coeff, 1st_order_coeff,...] for the last polynomial fit.
    :param direction: Which direction we are proceeding along the detector.
    :return: optimized coefficients, value of the sum of the fluxes across the trace, 0 or 1.
    0 if another maximum exists (e.g. another potential trace), 1 if no reasonable maximum exists.

    WARNING: It is hardcoded that the traces are spaced by no more than 100 pixels. If for some future detector they
    become further spaced, we will have to change the two instances of 100 here to that larger value. Although 100
    is like 10x the current trace spacing, so I do not forsee any issues.

    NOTE: there is no unit test for this, rather this is tested under an integration test for
    the do_stage of order-by-order fitting
    """
    coeffsguess, maximum_exists, refflux = generate_initial_guess_for_trace_polynomial(image_data, x,
                                                evaluated_legendre_polynomials, order=order,
                                                second_order_coefficient_guess=second_order_coefficient_guess,
                                                lastcoef=lastcoef, direction=direction)

    if maximum_exists:
        p1 = optimize.minimize(negative_flux_across_trace, coeffsguess, (imfilt, x, evaluated_legendre_polynomials), method='Powell').x

        summed_flux_value = -1 * negative_flux_across_trace(p1, imfilt, x, evaluated_legendre_polynomials)
    else:
        p1, summed_flux_value = lastcoef, refflux

    return p1, summed_flux_value, maximum_exists


def validate_fit_and_trim_erroneous_fits(last_fit_coefficients, all_fit_coefficients, loop_counter, length, maximum_exists, done, direction='up'):
    """
    :param last_fit_coefficients: coefficients for the trace fit (without index) of the last fit.
    :param all_fit_coefficients: list of coefficients for all the previous fits (with index).
    :param loop_counter: the counter i inside of find_all_traces_marching_up_or_down
    :param length: the height of the image.data (i.e. vertical extent of the image)
    :param maximum_exists: True or False whether generate_initial_guess_for_trace_polynomial decided there exists
                        another trace to fit at all
    :param done: whether we are done searching and fitting for traces.
    :param direction: the direction we are heading along in the detector. increasing y corresponds to 'up'.
    :return: the number of orders found (length of allcoefs after trimming), the trimmed coefficients which have:
            1. repeated fits removed
            2. fits which fall off the detector removed
            3. fits which fell backwards removed
    and done = True if we are indeed done.
    """
    num_of_orders_found = None
    if direction == 'up':
        if (last_fit_coefficients[0] < all_fit_coefficients[-2][1] or last_fit_coefficients[0] > length) and maximum_exists:
            all_fit_coefficients = all_fit_coefficients[:-1]  # delete bad fit
            num_of_orders_found = loop_counter - 1
            done = True

    if direction == 'down':
        if (last_fit_coefficients[0] < 0 or last_fit_coefficients[0] > all_fit_coefficients[-2][1]) and maximum_exists:
            all_fit_coefficients = all_fit_coefficients[:-1]  # delete bad fit
            num_of_orders_found = loop_counter - 1
            done = True

    if loop_counter >= 2 and maximum_exists and not done:
        if abs(last_fit_coefficients[0] - all_fit_coefficients[-2][1]) < 1 and abs(last_fit_coefficients[0] - all_fit_coefficients[-3][1]) < 1:
            all_fit_coefficients = all_fit_coefficients[:-2]  # delete repeated fits
            num_of_orders_found = loop_counter - 2
            done = True
    if not maximum_exists:
        done = True
        all_fit_coefficients = all_fit_coefficients[:-1]  # delete duplicate fit
        num_of_orders_found = loop_counter - 1
    return num_of_orders_found, all_fit_coefficients, done


def find_all_traces_marching_up_or_down(image_data, imfilt, x, summed_flux_values, evaluated_legendre_polynomials, length, order_of_poly_fit, last_fit_coefficients, all_fit_coefficients, direction='up'):
    """
    :param image_data: Ndarray. image_data e.g. image.data if image is a banzai image object.
    :param imfilt:
    :param x:
    :param summed_flux_values:
    :param evaluated_legendre_polynomials:
    :param length:
    :param order_of_poly_fit:
    :param last_fit_coefficients:
    :param all_fit_coefficients:
    :param direction:
    :return:

    NOTE: there is no unit test for this, rather this is tested under an integration test for
    the do_stage of order-by-order fitting
    """
    num_of_orders_found = 0
    done = False
    i = 1
    while not done:
        last_fit_coefficients, summed_flux_value, maximum_exists = find_order(image_data, imfilt, x, evaluated_legendre_polynomials, order=order_of_poly_fit,
                                                                lastcoef=last_fit_coefficients, direction=direction)
        summed_flux_values += [summed_flux_value]
        if direction == 'up':
            all_fit_coefficients += [[i] + list(last_fit_coefficients)]
        if direction == 'down':
            all_fit_coefficients += [[-i] + list(last_fit_coefficients)]
        num_of_orders_found, all_fit_coefficients, done = validate_fit_and_trim_erroneous_fits(last_fit_coefficients, all_fit_coefficients, i, length,
                                                                                               maximum_exists, done, direction=direction)
        i += 1
    return num_of_orders_found, all_fit_coefficients, last_fit_coefficients, summed_flux_values


def find_all_traces(image_data, imfilt, order_of_poly_fit, second_order_coefficient_guess):

    """
    :param image_data: Ndarray. Image_data e.g. image.data if image is a banzai image object.
    :param imfilt: ndarray, image.data passed through ndimage.spline_filter
    :param order_of_poly_fit: order of the polynomial fit.
    :param second_order_coefficient_guess: coefficient guess for the second order legendre polynomial which
    describe the traces across the ccd.
    :return:

    NOTE: there is no unit test for this, rather this is tested under an integration test for
    the do_stage of order-by-order fitting
    """
    length, width = image_data.shape

    evaluated_legendre_polynomials, x, xnorm = generate_legendre_array(width, order_of_poly_fit)

    # For the first coefficients, fit a quadratic
    # unconstrained.  Use this fitted quadratic to fit higher order
    # terms as desired, adding one order at a time to place as few
    # demands as possible on the optimization routine. Future orders are fit with the full
    # nth order polynomial (using the last fit as the initial guess) all at once.
    coef = None
    for i in range(2, order_of_poly_fit + 1):
        if coef is not None:
            coef = list(coef) + [0]
        coef, summed_flux_value, maximum_exists = find_order(image_data, imfilt, x, evaluated_legendre_polynomials, order=i,
                                               second_order_coefficient_guess=second_order_coefficient_guess, lastcoef=coef,
                                               direction='inplace')

    summed_flux_values = [summed_flux_value]
    initcoef = copy.deepcopy(coef)
    all_trace_coefficients = [[0] + list(coef)]

    ordersabove, all_trace_coefficients, coef, summed_flux_values = find_all_traces_marching_up_or_down(image_data, imfilt, x, summed_flux_values,
                                                                           evaluated_legendre_polynomials, length,
                                                                           order_of_poly_fit, coef, all_trace_coefficients,
                                                                           direction='up')

    ordersbelow, all_trace_coefficients, coef, summed_flux_values = find_all_traces_marching_up_or_down(image_data, imfilt, x, summed_flux_values,
                                                                           evaluated_legendre_polynomials, length,
                                                                           order_of_poly_fit, initcoef, all_trace_coefficients,
                                                                           direction='down')

    return all_trace_coefficients, summed_flux_values, ordersabove + ordersbelow + 1


def generate_legendre_array(image_width, order_of_poly_fits):
    """
    :param image_width: image.data.shape[1], x_extent of pixels.
    :param order_of_poly_fits: order of the polynomials used to fit the traces across the CCD.
    :return:
    """
    x = np.arange(image_width)
    xnorm = x * 2. / x[-1] - 1  # x normalized to run from -1 to 1

    # Set up Legendre polynomials to avoid roundoff error from
    # explicitly computing polynomials from their coefficients

    legendre_polynomial_array = np.ones((order_of_poly_fits + 1, x.shape[0]))
    for i in range(1, order_of_poly_fits + 1):
        legendre_polynomial_array[i] = np.polynomial.legendre.legval(xnorm, [0 for j in range(i)] + [1])
    return legendre_polynomial_array, x, xnorm


def totalflux_all_traces(coefficients_and_indices, image):
    """
    :param coefficients_and_indices: polynomial fit to traces
    :param image: banzai image object
    :return: total flux summed across all traces.
    The exact same thing as crosscoef except it generates the polynomial array in here and prefilters the image.
    """
    order_of_poly_fits = coefficients_and_indices.shape[1]-2
    legendre_array, x, xnorm = generate_legendre_array(image.data.shape[1], order_of_poly_fits)
    X = list(x)*(coefficients_and_indices.shape[0])
    #  X = [0,1,...,4095,0,1,..,4095,..]
    TraceYvals = np.dot(coefficients_and_indices[:, 1:], legendre_array).flatten()
    totalflux = np.sum(ndimage.map_coordinates(image.data.astype(np.float64), [TraceYvals, X], prefilter=True))
    return totalflux


def get_coefficients_from_meta(allmetacoeffs, stpolyarr):
    """
    NOTE: This is used in the suite of meta fit procedures (which are not implemented into Banzai-NRES as of
     11/13/2018) AND for generating realistic test frames for unit tests.
    :param allmetacoeffs: meta coefficients which describe the polynomial coefficients for each trace as a function
    of order.
    :param stpolyarr: The poly array which is the basis for the meta fit. Should be a legendre polynomial array.
    :return:
    """
    return np.dot(allmetacoeffs, stpolyarr).T


def legpolynomial(normxaxis, *metacoeffs):
    """
    legendre polynomial suitable for use in scipy.curve_fit
    :param normxaxis: the normalized axis from -1 to 1 which forms the domain of the legendre polynomial
    :param metacoeffs: metacoefficients cast as a tuple
    :return:
    """
    polyorder = len(metacoeffs)
    y = normxaxis * 0.
    for i in range(polyorder):
        y += metacoeffs[i] * np.polynomial.legendre.legval(normxaxis, [0 for j in range(i)] + [1])
    return y


def extract_coeffs_entire_lampflat_frame(image, order_of_poly_fits, second_order_coefficient_guess):
    """
    This extracts the trace coefficients for each bright order of a frame. This is only stable for lampflat frames.
    It returns a list of the coefficients, ordered arbitrarily (fibers are not separated). It also returns the summed fluxed across each order
    called val, and the total number of orders found by the algorithm.
    Parameters:
        image : Banzai Image object.
        order_of_poly_fits : order of the polynomial fit per trace.
        second_order_coefficient_guess : coefficient guess for the second order legendre polynomial which
        describe the traces across the ccd.

    NOTE: there is no unit test for this, rather this is tested under an integration test for
    the do_stage of order-by-order fitting

    """
    imagefiltered = ndimage.spline_filter(image.data)

    # finding coefficients of traces which fit the echelle orders across the CCD.
    allcoef, vals, totalnumberoforders = find_all_traces(image.data, imagefiltered, order_of_poly_fits, second_order_coefficient_guess)
    sortedallcoefs = np.array(allcoef)[np.array(allcoef)[:, 0].argsort()]
    order_indices = np.arange(totalnumberoforders)
    # appending indices 0,1,2...,totalnumberoforders as the first column. prior it is indexed from negative numbers.

    coefficients_and_indices = np.insert(sortedallcoefs[:, 1:], obj=0, values=order_indices, axis=1)

    return coefficients_and_indices, vals, totalnumberoforders


def exclude_traces_which_jet_off_detector(coefficients_and_indices, image):
    """
    :param coefficients_and_indices: ndarray. list of trace polynomial coefficients
    :param image: Banzai image object.
    :return: coefficients_and_indices excluding any traces which have discontinuity because they fall off the detector.
    """
    order_of_poly_fits = coefficients_and_indices.shape[1] - 2
    legendre_polynomial_array, not_needed, not_needed_2 = generate_legendre_array(image.data.shape[1],
                                                                                  order_of_poly_fits)
    trace_values_versus_xpixel = np.dot(coefficients_and_indices[:, 1:], legendre_polynomial_array)
    # trim any traces which are not contiguously on the detector.
    coefficients_and_indices = coefficients_and_indices[np.all(trace_values_versus_xpixel > 0, axis=1)]
    return coefficients_and_indices


def split_and_sort_coefficients_for_each_fiber(coefficients_and_indices, num_lit_fibers):
    """
    :param coefficients_and_indices: ndarray. list of trace polynomial coefficients as-is from blind-fitting (e.g. 134 traces,
    ordered by occurrence, starting from the bottom (red half) of the detector.
    :param num_lit_fibers: Int. The number of fibers which are lit on the detector (e.g. 2)
    :return:
    """
    # cutting of order index column
    coefficients_no_indices = coefficients_and_indices[:, 1:]

    # ensuring an even number of traces in the coefficients so that they can be split easily
    while coefficients_no_indices.shape[0] % num_lit_fibers != 0:
        coefficients_no_indices = coefficients_no_indices[:-1]
    num_orders = int(coefficients_no_indices.shape[0] / num_lit_fibers)
    order_indices = [i for i in range(num_orders)]*num_lit_fibers

    fiber_coefficients = (coefficients_no_indices[::num_lit_fibers],)
    for i in range(1, num_lit_fibers):
        fiber_coefficients += (coefficients_no_indices[i::num_lit_fibers],)

    ordered_coefficients = np.vstack(fiber_coefficients)
    coefficients_and_indices = np.insert(ordered_coefficients, obj=0, values=order_indices, axis=1)

    return coefficients_and_indices


def fit_traces_order_by_order(image, second_order_coefficient_guess, order_of_poly_fits=4, num_lit_fibers=2):
    """
    :param image: Banzai image object
    :param second_order_coefficient_guess: guess for the coefficient of the second order legendre polynomial for
    the blind fit.
    :param order_of_poly_fits: Highest order of the polynomial fit to each trace. 4 is good. Do not change needlessly.
    :return array of trace fit coefficients arranged such that those for the first fiber are first.
    the first 67 rows of the array. fiber designation is arbitrary at this point.
    """
    coefficients_and_indices, vals, totalnumberoftraces = extract_coeffs_entire_lampflat_frame(image, order_of_poly_fits,
                                                                                          second_order_coefficient_guess)

    coefficients_and_indices = exclude_traces_which_jet_off_detector(coefficients_and_indices, image)
    coefficients_and_indices = split_and_sort_coefficients_for_each_fiber(coefficients_and_indices, num_lit_fibers)

    logger.debug('%s traces found' % coefficients_and_indices.shape[0])

    return coefficients_and_indices


def get_number_of_lit_fibers(image):
    """
    :param image: banzai image
    :return: the number of lit fibers (e.g. the number of entries in tung&tung&.... which are not 'none')
    """
    if image.header.get('OBJECTS') is None:
        logger.error('header keyword OBJECTS not found, cannot get the number of lit fibers.')
        return None
    fiber_info = image.header.get('OBJECTS').split('&')
    num_unlit_fibers = fiber_info.count('none')
    num_lit_fibers = int(len(fiber_info) - num_unlit_fibers)
    return num_lit_fibers
