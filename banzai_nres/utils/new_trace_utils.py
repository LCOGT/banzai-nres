"""
trace_utils.py: Routines for finding echelle orders across a CCD.
Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)

    Tim Brandt (tbrandt@physics.ucsb.edu)
"""

import numpy as np
from scipy import ndimage, optimize
from astropy.table import Table
import copy

import logging

logger = logging.getLogger(__name__)


class Trace(object):
    """
    :param trace_centers = {'id': ndarray, 'centers': ndarray}. 'centers' gives a 2d array, where
    the jth row are the y centers across the detector for the trace with identification trace_centers['id'][j]
    :param design_matrix = the 2d array such that coefficients dot design_matrix gives the trace centers for all the orders.
    """
    def __init__(self, trace_centers=None, trace_table_name=None, second_order_coefficient_guess=None, poly_fit_order=None):
        self.trace_table_name = trace_table_name # will be taken from settings. when instantiated in TraceMaker
        self.trace_centers = trace_centers
        self.second_order_coefficient_guess = second_order_coefficient_guess
        self.poly_fit_order = poly_fit_order
        self.design_matrix = None # consider putting this into a TraceFItter class.
        self.x_norm = None # consider putting this into a TraceFItter class.
        self.x = None # consider putting this into a TraceFItter class.

    def get_trace_centers(self, row):
        return self.trace_centers['centers'][row]

    @staticmethod
    def fit_traces(cls, image, poly_fit_order, second_order_coefficient_guess):
        trace = cls(trace_centers={'id': [], 'centers': []},
                    poly_fit_order=poly_fit_order,
                    second_order_coefficient_guess=second_order_coefficient_guess)
        start_point = image.data.shape[0]/3
        trace_fitter = SingleTraceFitter(image_data=image.data, start_point=start_point,
                                         second_order_coefficient_guess=second_order_coefficient_guess,
                                         poly_fit_order=poly_fit_order)
        return trace

    def write(self, filename):
        table = Table(self.trace_centers, name=self.trace_table_name)
        table['id'].description = 'Blah'
        table.write(filename)


class SingleTraceFitter(object):
    #TODO add the other non necessary arguments as keyword arguments for extraargs.
    def __init__(self, image_data=None, poly_fit_order=2, start_point=None, second_order_coefficient_guess=None,
                 extraargs={}):
        self.second_order_coefficient_guess = second_order_coefficient_guess
        self.start_point = start_point
        self.image_data = image_data
        self.poly_fit_order = poly_fit_order
        self.filtered_image_data = extraargs.get('filtered_image_data')
        self.coefficients = extraargs.get('coefficients')
        self.initial_guess_next_fit = extraargs.get('initial_guess_next_fit')
        self.x = extraargs.get('x')
        self.x_norm = extraargs.get('xnorm')
        self.design_matrix = extraargs.get('design_matrix')
        if extraargs.get('initialize_fit_objects', True) is True:
            self._initialize_fit_objects()

    def match_filter_to_refine_initial_guess(self):
        return None

    def generate_initial_guess_after_turn_around(self):
        self.initial_guess_next_fit = self.coefficients[0]

    def generate_initial_guess(self):
        if self.start_point is None or self.second_order_coefficient_guess is None:
            logger.error('Starting y position up the detector nor the second order guess have been specified '
                         ', failed to generate an initial guess for trace fitting.')
        else:
            self.initial_guess_next_fit = np.zeros(self.poly_fit_order+1).astype(np.float64)
            self.initial_guess_next_fit[0] = self.start_point
            self.initial_guess_next_fit[2] = self.second_order_coefficient_guess

    def _centers_from_coefficients(self, single_trace_coefficients):
        #trace_centers = self.x * 0.
        #for i in range(len(self.design_matrix)):
        #    trace_centers += self.design_matrix[i] * single_trace_coefficients[i]
        trace_centers = np.dot(single_trace_coefficients, self.design_matrix)
        return trace_centers

    def _flux_across_trace(self, trace_centers):
        total_flux = np.sum(ndimage.map_coordinates(self.filtered_image_data, [trace_centers, self.x], prefilter=False))
        return total_flux

    def _initialize_fit_objects(self):
        self._normalize_domain_coordinates()
        self.filtered_image_data = self._prefilter_image_data(self.image_data)
        self.design_matrix = self._generate_design_matrix(self.x_norm, self.poly_fit_order)

    def _normalize_domain_coordinates(self):
        self.x = np.arange(self.image_data.shape[1])
        self.x_norm = self.x * 2. / self.x[-1] - 1  # x normalized to run from -1 to 1

    @staticmethod
    def _prefilter_image_data(image_data):
        return ndimage.spline_filter(image_data)

    @staticmethod
    def _generate_design_matrix(normalized_domain, poly_fit_order):
        design_matrix = np.ones((poly_fit_order + 1, normalized_domain.shape[0]))
        for i in range(poly_fit_order + 1):
            design_matrix[i] = legendre(order=i, values=normalized_domain)
        return design_matrix


def legendre(order, values):
    #TODO replace with scipy legendre(order)(values) if it works the same.
    if order == 0:
        return np.ones_like(values)
    else:
        return np.polynomial.legendre.legval(values, [0 for j in range(order)] + [1])


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


def get_trace_centers_from_coefficients(image_width, coefficients_and_indices):
    #TODO this will likely be a test_util function since we will not use
    """
    :param coefficients_and_indices: polynomial fit coefficients which describe the traces. Legendre polynomials
                                    normalized between -1 and 1.
    :param image_width: image.data.shape[1]
    :return: trace centroids for each trace, versus x pixel. E.g. trace_values_versus_xpixel[2,5] is the 3rd orders value
            at x=5.
            num_traces = num_orders*Num_fibers
            x = [0,1,2,...,image.data.shape[1]-1]
    """
    coeflen, coefwidth = coefficients_and_indices.shape
    num_traces, order_of_poly_fits = coeflen, coefwidth - 2
    legendre_polynomial_array, x, xnorm = generate_legendre_array(image_width, order_of_poly_fits)
    trace_values_versus_xpixel = np.dot(coefficients_and_indices[:, 1:], legendre_polynomial_array)
    return trace_values_versus_xpixel, num_traces, x


def extract_coeffs_entire_lampflat_frame(image_data, order_of_poly_fits, second_order_coefficient_guess):
    """
    This extracts the trace coefficients for each bright order of a frame. This is only stable for lampflat frames.
    It returns a list of the coefficients, ordered arbitrarily (fibers are not separated). It also returns the summed fluxed across each order
    called val, and the total number of orders found by the algorithm.
    Parameters:
        image_data: Image.data of banzai image object (e.g. the data portion of a fits frame)
        order_of_poly_fits : order of the polynomial fit per trace.
        second_order_coefficient_guess : coefficient guess for the second order legendre polynomial which
        describe the traces across the ccd.

    NOTE: there is no unit test for this, rather this is tested under an integration test for
    the do_stage of order-by-order fitting

    """
    imagefiltered = ndimage.spline_filter(image_data)

    # finding coefficients of traces which fit the echelle orders across the CCD.
    allcoef, vals, totalnumberoforders = find_all_traces(image_data, imagefiltered, order_of_poly_fits, second_order_coefficient_guess)
    sortedallcoefs = np.array(allcoef)[np.array(allcoef)[:, 0].argsort()]
    order_indices = np.arange(totalnumberoforders)
    # appending indices 0,1,2...,totalnumberoforders as the first column. prior it is indexed from negative numbers.

    coefficients_and_indices = np.insert(sortedallcoefs[:, 1:], obj=0, values=order_indices, axis=1)

    return coefficients_and_indices, vals, totalnumberoforders


def exclude_traces_which_jet_off_detector(coefficients_and_indices, image_data):
    """
    :param coefficients_and_indices: ndarray. list of trace polynomial coefficients
    :param image_data: Image.data of banzai image object (e.g. the data portion of a fits frame)
    :return: coefficients_and_indices excluding any traces which have discontinuity because they fall off the detector.
    """
    order_of_poly_fits = coefficients_and_indices.shape[1] - 2
    legendre_polynomial_array, not_needed, not_needed_2 = generate_legendre_array(image_data.shape[1],
                                                                                  order_of_poly_fits)
    trace_values_versus_xpixel = np.dot(coefficients_and_indices[:, 1:], legendre_polynomial_array)
    # trim any traces which are not contiguously on the detector.
    coefficients_and_indices = coefficients_and_indices[np.all(trace_values_versus_xpixel > 0, axis=1)]
    return coefficients_and_indices


def fit_traces_order_by_order(image_data, second_order_coefficient_guess, order_of_poly_fits=4):
    """
    :param image_data: Image.data of banzai image object (e.g. the data portion of a fits frame)
    :param second_order_coefficient_guess: guess for the coefficient of the second order legendre polynomial for
    the blind fit.
    :param order_of_poly_fits: Highest order of the polynomial fit to each trace. 4 is good. Do not change needlessly.
    :return array of trace fit coefficients arranged such that those for the first fiber are first.
    the first 67 rows of the array. fiber designation is arbitrary at this point.
    """
    coefficients_and_indices, vals, totalnumberoftraces = extract_coeffs_entire_lampflat_frame(image_data, order_of_poly_fits,
                                                                                          second_order_coefficient_guess)

    coefficients_and_indices = exclude_traces_which_jet_off_detector(coefficients_and_indices, image_data)

    logger.debug('%s traces found' % coefficients_and_indices.shape[0])

    return coefficients_and_indices


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
    #TODO this should be a test_util function.
    """
    NOTE: This is used in the suite of meta fit procedures (which are not implemented into Banzai-NRES as of
     11/13/2018) AND for generating realistic test frames for unit tests.
    :param allmetacoeffs: meta coefficients which describe the polynomial coefficients for each trace as a function
    of order.
    :param stpolyarr: The poly array which is the basis for the meta fit. Should be a legendre polynomial array.
    :return:
    """
    return np.dot(allmetacoeffs, stpolyarr).T



