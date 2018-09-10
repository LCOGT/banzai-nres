#! /usr/bin/python3
"""
trace_utils.py: Routines for finding echelle orders across a CCD.
Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)

    Tim Brandt (tbrandt@physics.ucsb.edu)
"""

import numpy as np
from scipy import ndimage, optimize, interpolate
from scipy.optimize import curve_fit
import copy
import itertools
from banzai import logs

logger = logs.get_logger(__name__)


def maxima(A, s, k, ref):
    """
    A procedure for finding maxima which are not just peaks in the noise floor. E.g. searches for statistically
    significant maxima.
    :param A: An array of values where you wish to find maxima
    :param s: The window were the central-point's value must be larger than all the other values in the window to qualify
              as a maximum.
    :param k: multiplicative constant which qualifies a point as a real maxima (not just noise)
    :param ref: the reference value which k*r defines the value any real maximal point must exceed
    :return: The index of the a point near a maximum and the value at the maximum, and a boolean for whether or not
            a valid maximum exists.
    """

    i = s - 1
    l, r = 0, 0
    A = A - np.ones_like(A) * np.min(A)
    threshold = k * ref
    firstmax = [0, 0]
    maximum_exists = False
    while i < len(A) - s and int(l + r) != 2:
        l, r = 1, 1
        for j in range(1, s + 1):
            if A[i] < A[i + j] or A[i] < threshold:
                r, j = 0, s + 1
            if A[i] < A[i - j] or A[i] < threshold:
                l, j = 0, s + 1
        if int(l + r) == 2:
            firstmax = [i, A[i]]
            maximum_exists = True
        i += 1
    return firstmax, maximum_exists


def crosscoef(legendre_polynomial_coefficients, imfilt, x, evaluated_legendre_polynomials):
    """
    :param legendre_polynomial_coefficients: list of legendre polynomial coefficients
    :param imfilt: ndarray, image.data passed through ndimage.spline_filter
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param evaluated_legendre_polynomials: Legendre polynomial evaluated over -1 to 1 , ndarray.
    :return: The negative of the total flux summed across the trace.
    """
    y = x * 0.

    for i in range(len(legendre_polynomial_coefficients)):
        y += evaluated_legendre_polynomials[i] * legendre_polynomial_coefficients[i]
    val = np.sum(ndimage.map_coordinates(imfilt, [y, x], prefilter=False))
    return -1 * val


def fluxvalues(testpoints, legendre_polynomial_coefficients, imfilt, x, evaluated_legendre_polynomials):
    """
    :param legendre_polynomial_coefficients: list of legendre polynomial coefficients excluding 0th order coefficients
    :param imfilt: ndarray, image.data passed through ndimage.spline_filter
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param evaluated_legendre_polynomials: Legendre polynomial evaluated over -1 to 1 , ndarray.
    :return: The positive total flux at each point in testpoints
    """
    values = np.zeros_like(testpoints)
    for i in range(0, len(values)):
        coeffguesses = list([testpoints[i]] + legendre_polynomial_coefficients)
        values[i] = (-1) * crosscoef(coeffguesses, imfilt, x, evaluated_legendre_polynomials)
    return values


def generate_initial_guess_for_trace_polynomial(image, imfilt, x, evaluated_legendre_polynomials, order=2, second_order_coefficient_guess=90, lastcoef=None, direction='up'):
    """
    :param image:
    :param imfilt:
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
        p[0] = int(imfilt.shape[0]/3)
        if order >= 2:
            p[2] = second_order_coefficient_guess
        coeffsguess = copy.deepcopy(p)
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

        fluxvals = fluxvalues(testpoints, p, image.data, x, evaluated_legendre_polynomials)
        refflux = max((-1) * crosscoef(lastcoef, image.data, x, evaluated_legendre_polynomials), max(fluxvals))

        deltap0guess, maximum_exists = maxima(fluxvals, 5, 1 / 20, refflux)[0]

        if direction == 'up':
            p0 = deltap0guess + min(testpoints)
        elif direction == 'down':
            p0 = max(testpoints) - deltap0guess
        elif direction == 'inplace':
            p0 = deltap0guess + min(testpoints)

        coeffsguess = [p0] + list(p)
    return coeffsguess, maximum_exists, refflux


def findorder(image, imfilt, x, evaluated_legendre_polynomials, order=2, second_order_coefficient_guess=90, lastcoef=None, direction='up'):
    """
    :param image: banzai image object.
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
    """
    coeffsguess, maximum_exists, refflux = generate_initial_guess_for_trace_polynomial(image, imfilt, x,
                                                evaluated_legendre_polynomials, order=order,
                                                second_order_coefficient_guess=second_order_coefficient_guess,
                                                lastcoef=lastcoef, direction=direction)

    if maximum_exists:
        p1 = optimize.minimize(crosscoef, coeffsguess, (imfilt, x, evaluated_legendre_polynomials), method='Powell').x

        val = -1 * crosscoef(p1, imfilt, x, evaluated_legendre_polynomials)
    else:
        p1, val = lastcoef, refflux

    return p1, val, maximum_exists


def tracesacrossccd(image, imfilt, order_of_poly_fit, second_order_coefficient_guess):
    """
    :param image: banzai image object
    :param imfilt: ndarray, image.data passed through ndimage.spline_filter
    :param order_of_poly_fit: order of the polynomial fit.
    :param second_order_coefficient_guess: coefficient guess for the second order legendre polynomial which
    describe the traces across the ccd.
    :return:
    """
    length, width = image.data.shape

    evaluated_legendre_polynomials, x, xnorm = generate_legendre_array(image, order_of_poly_fit)

    # For the first coefficients, fit a quadratic
    # unconstrained.  Use this fitted quadratic to fit higher order
    # terms as desired, adding one order at a time to place as few
    # demands as possible on the optimization routine. Future orders are fit with the full
    # nth order polynomial (using the last fit as the initial guess) all at once.
    coef = None
    for i in range(2, order_of_poly_fit + 1):
        if coef is not None:
            coef = list(coef) + [0]
        coef, val, maximum_exists = findorder(image, imfilt, x, evaluated_legendre_polynomials, order=i,
                                     second_order_coefficient_guess=second_order_coefficient_guess, lastcoef=coef,
                                     direction='inplace')

    vals = [val]
    initcoef = copy.deepcopy(coef)
    allcoef = [[0] + list(coef)]

    # March up in the orders.

    # the search stops if the 0th order coefficient computed is less than 0 or is less than the previous coefficient computed (i.e. down step).
    # Or is within a pixel of the previous two computed coefficients.
    # So 3 zero-order coefficients are equal. In the latter case, we delete two of the last fits, and then stop the search.
    # if the sum across the trace for the found order k is less than 1/10 of the sum of the k-2 order (belonging
    # to the same fiber, then we stop the search.

    done = 0
    i = 1
    while done == 0:
        coef, val, maximum_exists = findorder(image, imfilt, x, evaluated_legendre_polynomials, order=order_of_poly_fit,
                                     lastcoef=coef, direction='up')
        vals += [val]
        allcoef += [[i] + list(coef)]
        if coef[0] < allcoef[-2][1] or coef[0] > length and maximum_exists:
            allcoef = allcoef[:-1]  # delete bad fit
            ordersabove = i - 1
            done = 1
        if i >= 2 and done == 0:
            if abs(coef[0] - allcoef[-2][1]) < 1 and abs(coef[0] - allcoef[-3][1]) < 1 and maximum_exists:
                allcoef = allcoef[:-2]  # delete repeated fits
                print(allcoef[-1][1], coef[0])
                ordersabove = i - 2
                done = 1
        if not maximum_exists:
            done = 1
            allcoef = allcoef[:-1]  # delete duplicate fit
            ordersabove = i - 1
        i += 1

    # Now march down to the bottom of the array.  Resume from the
    # initial coefficients from earlier.

    # the search stops if the 0th order coefficient computed is less than 0 or is greater than the previous coefficient computed (i.e. up step).
    # Or is within a pixel of the previous two computed coefficients,
    # So 3 zero-order coefficients are equal. In the latter case, we delete two of the last fits, and then stop the search.

    done = 0
    i = 1
    coef = initcoef
    while done == 0:
        coef, val, maximum_exists = findorder(image, imfilt, x, evaluated_legendre_polynomials,
                                     order=order_of_poly_fit, lastcoef=coef,
                                     direction='down')
        vals += [val]
        allcoef += [[-i] + list(coef)]

        if coef[0] < 0 or coef[0] > allcoef[-2][1] and maximum_exists:
            allcoef = allcoef[:-1]  # delete bad fit
            ordersbelow = i - 1
            done = 1
        if i >= 2 and done == 0 and maximum_exists:
            if abs(coef[0] - allcoef[-2][1]) < 1 and abs(coef[0] - allcoef[-3][1]) < 1:
                allcoef = allcoef[:-2]  # delete repeated fits
                print(allcoef[-1][1], coef[0])
                ordersbelow = i - 2
                done = 1
        if not maximum_exists:
            done = 1
            allcoef = allcoef[:-1]  # delete duplicate fit
            ordersbelow = i - 1
        i += 1

    return allcoef, vals, ordersabove + ordersbelow + 1


def generate_legendre_array(image, order_of_poly_fits):
    """
    :param image: Banzai image object.
    :param order_of_poly_fits: order of the polynomials used to fit the traces across the CCD.
    :return:
    """
    x = np.arange(image.data.shape[1])
    xnorm = x * 2. / x[-1] - 1  # x normalized to run from -1 to 1

    # Set up Legendre polynomials to avoid roundoff error from
    # explicitly computing polynomials from their coefficients

    legendre_polynomial_array = np.ones((order_of_poly_fits + 1, x.shape[0]))
    for i in range(1, order_of_poly_fits + 1):
        legendre_polynomial_array[i] = np.polynomial.legendre.legval(xnorm, [0 for j in range(i)] + [1])
    return legendre_polynomial_array, x, xnorm


def check_for_close_fit(coefficients_and_indices_list, images, num_lit_fibers, max_pixel_error=1E-1):
    """
    :param coefficients_and_indices_list: list of trace_coefficients across the detector for multiple different fits
     . I.e. a list of ndarray with the first column 0,1,2,..66,0,1.. the fiber indexes, and the second column
            the 0th order coefficients for that order trace. The fibers are arranged fiber_order[0] then fiber_order[1].
            as listed in the attribute Image().fiber_order
    :param images: List of Banzai Image objects
    :param max_pixel_error: Max allowed y pixel error between the two sets of traces. Computed at every x value But only for
                        the central orders. E.g. orders 10-50 for 67 orders per fiber.
    :return: True if close, False if not.
    """

    trace_values_versus_xpixel_list = []
    num_traces_list = []
    for coefficients, image in zip(coefficients_and_indices_list, images):
        trace_values_versus_xpixel, num_traces, x = get_trace_centroids_from_coefficients(coefficients, image)
        num_traces_list += [num_traces]
        trace_values_versus_xpixel_list += [trace_values_versus_xpixel]
        image_shape = image.data.shape
    assert num_traces_list[1:] == num_traces_list[:-1]
    assert num_traces % num_lit_fibers == 0
    num_orders = int(num_traces/num_lit_fibers)
    order_buffer = int(num_orders/5)
    x_buffer = int(image_shape[1]/4)

    # computing absolute differences between the trace centroid locations at every x value.
    trace_values_versus_xpixel_arr = np.array(trace_values_versus_xpixel_list)
    error_between_fits = np.abs((trace_values_versus_xpixel_arr - np.mean(trace_values_versus_xpixel_arr, axis=0)))
    # restricting region where we care about differences (center of detector)
    select_errors_first_fiber = error_between_fits[:, order_buffer:(num_orders-order_buffer), x_buffer:(image_shape[1] - x_buffer)]
    select_errors_second_fiber = error_between_fits[:, num_orders + order_buffer:(num_traces-order_buffer), x_buffer:(image_shape[1] - x_buffer)]

    max_error_between_fits = max(np.max(select_errors_first_fiber), np.max(select_errors_second_fiber))
    if max_error_between_fits < max_pixel_error:
        return True
    else:
        print('warning! central trace centroids between test lampflats disagreed \n '
              'beyond max allowed error of %s pixels'%max_pixel_error)
        return False


def check_flux_change(coefficients_and_indices_new, coefficients_and_indices_initial, image):
    """
    :param coefficients_and_indices_new: polynomial fit to traces that is new
    :param coefficients_and_indices_initial: polynomial fit to traces that is old (e.g. from a master file)
    :param image: Banzai image object
    :return: True if the flux change is less than 1 percent.
    """
    initial_flux = totalflux_all_traces(coefficients_and_indices_initial, image)
    delta_fraction_flux = (totalflux_all_traces(coefficients_and_indices_new, image) - initial_flux)/initial_flux
    logger.debug('(new_flux - master_cal_flux)/master_cal_flux) = %s' % delta_fraction_flux)
    if np.abs(delta_fraction_flux) < 1E-2:
        return True
    else:
        return False


def totalflux_all_traces(coefficients_and_indices, image):
    """
    :param coefficients_and_indices: polynomial fit to traces
    :param image: banzai image object
    :return: total flux summed across all traces.
    """
    order_of_poly_fits = coefficients_and_indices.shape[1]-2
    legendre_array, x, xnorm = generate_legendre_array(image, order_of_poly_fits)
    X = list(x)*(coefficients_and_indices.shape[0])
    #  X = [0,1,...,4095,0,1,..,4095,..]
    TraceYvals = np.dot(coefficients_and_indices[:, 1:], legendre_array).flatten()
    totalflux = np.sum(ndimage.map_coordinates(image.data.astype(np.float64), [TraceYvals, X], prefilter=True))
    return totalflux


def get_trace_centroids_from_coefficients(coefficients_and_indices, image):
    """
    :param coefficients_and_indices: polynomial fit to traces
    :param image: banzai image object
    :return: trace centroids for each trace, versus x pixel. E.g. trace_values_versus_xpixel[2,5] is the 3rd orders value
            at x=5.
            num_traces = num_orders*Num_fibers
            x = [0,1,2,...,image.data.shape[1]-1]
    """
    coeflen, coefwidth = coefficients_and_indices.shape
    num_traces, order_of_poly_fits = coeflen, coefwidth - 2
    legendre_polynomial_array, x, xnorm = generate_legendre_array(image, order_of_poly_fits)
    trace_values_versus_xpixel = np.dot(coefficients_and_indices[:, 1:], legendre_polynomial_array)
    return trace_values_versus_xpixel, num_traces, x


def find_nearest(array, value):
    """
    :param array: 1d numpy array or list.
    :param value: value with which you wish to find the i such that array[i] is closest to value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class ImageSplines(object):
    """
    Stores the arrays of scipy-spline objects which give the derivatives and values at points
    in the image.
    """
    def __init__(self, spline, first_derivative, second_derivative):
        self.spline = spline
        self.first_derivative = first_derivative
        self.second_derivative = second_derivative


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


def get_coefficients_from_meta(allmetacoeffs, stpolyarr):
    """
    :param allmetacoeffs: meta coefficients
    :param stpolyarr: The poly array which is the basis for the meta fit. Should be a legendre polynomial array.
    :return:
    """
    return np.dot(allmetacoeffs, stpolyarr).T


def metacoefficients(orderarray, trace_coefficients, metapolyorder):
    """
    :param orderarray: The order index list, nd array
    :param trace_coefficients: Coefficients (without  the order index list as first column)
    :param metapolyorder: order of the polynomial fit to how the coefficients vary as a function of order number
    :return:
    """
    allmetacoeffs = []
    metacoeffsinit = tuple([0] * (metapolyorder + 1))
    for j in range(0, trace_coefficients.shape[1]):
        popt, pcov = curve_fit(legpolynomial, orderarray, trace_coefficients[:, j], p0=metacoeffsinit)
        allmetacoeffs += [list(popt)]
    return np.array(allmetacoeffs)


def neg_totalflux_for_scipy(coeffs_vector, *extraargs):
    """
    :param coeffs_vector: Vector of meta coefficients.
    :param extraargs: (image_splines, stpolyarr, legpolyarr, pixelxarray, imfilt)
    :returns
        the negative of the total flux. Calculated using the scipy spline object stored in image_splines. Not ndimage
        like totalflux_all_traces.
    """
    image_splines, stpolyarr, legpolyarr, x = extraargs

    # get order of poly coefficients.
    tracepolyorder = legpolyarr.shape[0] - 1
    metapolyorder = stpolyarr.shape[0] - 1

    # reshape meta coeffs.
    metacoeffs_temp = coeffs_vector.reshape((tracepolyorder + 1), (metapolyorder + 1))
    # construct y values of the polynomial traces from the metacoefficients
    tracecoeffs = get_coefficients_from_meta(metacoeffs_temp, stpolyarr)
    # evaluate each trace at every x pixel
    traces = np.dot(tracecoeffs, legpolyarr)
    pixelyarray = np.copy(x)

    flux_values_along_traces = np.array([image_splines.spline[i](traces[:, i]) for i in pixelyarray]).T
    return (-1) * np.sum(flux_values_along_traces)


def NegativeHessian(coeffs_vector, *args):
    """
    :param coeffs_vector: nd-array list of the n-meta coefficients
    :param args: tuple type, (image_splines, stpolyarr, legpolyarr, pixelxarray) .
    See Newtonssmethod_meta_fit below for parameter designations
    :return: the negative of the Hessian, an nd array shape (n,n)
    """
    image_splines, stpolyarr, legpolyarr, pixelxarray = args

    # get order of poly coefficients.
    tracepolyorder = legpolyarr.shape[0] - 1
    metapolyorder = stpolyarr.shape[0] - 1

    # reshape meta coeffs.
    metacoeffs_temp = coeffs_vector.reshape((tracepolyorder + 1), (metapolyorder + 1))  # first column is order number

    # construct y values of the polynomial traces from the metacoefficients for this iteration
    tracecoeffs = get_coefficients_from_meta(metacoeffs_temp, stpolyarr)  # first column is order number
    # evaluate each trace at every x pixel
    traces = np.dot(tracecoeffs, legpolyarr)
    pixelyarray = np.copy(pixelxarray)

    secd = np.array([image_splines.second_derivative[i](traces[:, i]) for i in pixelyarray]).T

    # construct the filled Hessian for all traces
    Hessians = np.array([np.dot(legpolyarr * secd[i], legpolyarr.T) for i in
                         range(traces.shape[0])])

    Hessian = []
    for j in range(tracepolyorder + 1):
        for p in range(metapolyorder + 1):
            for k in range(tracepolyorder + 1):
                for q in range(metapolyorder + 1):
                    Hessian += [np.sum(stpolyarr[p] * stpolyarr[q] * Hessians[:, j, k])]
                    # note: from np.dot, there is an absolute error of 1e-9
                    # between some Hessians[:,j,k] - Hessians[:,k,j]

    Hessian = np.array(Hessian)
    Hessian = Hessian.reshape((metapolyorder + 1) * (tracepolyorder + 1), (metapolyorder + 1) * (tracepolyorder + 1))

    # Setting so that the Hessian is perfectly symmetric.
    Hessian = np.tril(Hessian, k=0)  # zero the elements in the upper triangle - they differ due to np.dot errors.
    Hessian = Hessian + Hessian.T - np.diag(np.diag(Hessian))  # calculate the fully symmetric Hessian.
    return (-1) * Hessian


def NegativeGradient(coeffs_vector, *args):
    """
    :param coeffs_vector: nd-array list of the n-meta coefficients
    :param args: tuple type, (image_splines, stpolyarr, legpolyarr, pixelxarray) .
    See Newtonssmethod_meta_fit below for parameter designations
    :return: the negative of the gradient.
    """
    image_splines, stpolyarr, legpolyarr, pixelxarray = args

    # get order of poly coefficients.
    tracepolyorder = legpolyarr.shape[0] - 1
    metapolyorder = stpolyarr.shape[0] - 1

    # reshape meta coeffs.
    metacoeffs_temp = coeffs_vector.reshape((tracepolyorder + 1), (metapolyorder + 1))
    # construct y values of the polynomial traces from the metacoefficients
    tracecoeffs = get_coefficients_from_meta(metacoeffs_temp, stpolyarr)
    # evaluate each trace at every x pixel
    traces = np.dot(tracecoeffs, legpolyarr)
    pixelyarray = np.copy(pixelxarray)

    firstd = np.array([image_splines.first_derivative[i](traces[:, i]) for i in pixelyarray]).T

    # construct the filled gradient
    gradients = np.dot(firstd, legpolyarr.T)  # evaluating the gradient for all traces
    grad = []
    for k in range(tracepolyorder + 1):
        for j in range(metapolyorder + 1):
            grad += [np.sum(stpolyarr[j] * gradients[:, k])]
    grad = np.array(grad)
    return (-1) * grad


def meta_fit(metacoeffs, stpolyarr, legpolyarr, image_splines, pixelxarray):
    """
    Evaluates meta-coefficients which optimize the flux across the entire detector.
    Best method is Newton-CG because it is well suited theoretically to this problem.
    Performance notes : On a 500x500 test frame with 20 traces, this function takes 10 ms to evaluate the Hessian,
    and 5 ms to evaluate the gradient. It calls the grad and hessian ~5 and ~20 times respectively, and evalutes the
    function ~20 times. function evaluations take about 5 ms.
    The total meta fit time is roughly 150 ms on said test frame.
    """
    metacoeffsinitial = np.copy(metacoeffs).reshape(metacoeffs.size)
    extraargs = (image_splines, stpolyarr, legpolyarr, pixelxarray)
    tracepolyorder = legpolyarr.shape[0] - 1
    metapolyorder = stpolyarr.shape[0] - 1
    metacoeffsnew = optimize.minimize(neg_totalflux_for_scipy, metacoeffsinitial, args=extraargs, method='NEWTON-CG'
                                      , jac=NegativeGradient, hess=NegativeHessian, options={'disp': False}).x

    return metacoeffsnew.reshape((tracepolyorder + 1), (metapolyorder + 1))


def cross_correlate_image_indices(images, cross_correlate_num):
    """
    :param images: Banzai image objects
    :param cross_correlate_num: Number of images to pair together in the list of combinations
    :return: All unique combinations of the indices of images from image
    """
    assert len(images) >= 1

    image_indices_to_try = list(range(len(images)))
    try_combinations_of_images = False

    if len(images) >= cross_correlate_num >= 2:
        image_indices_to_try = list(itertools.combinations(range(len(images)), cross_correlate_num))
        try_combinations_of_images = True
    return image_indices_to_try, try_combinations_of_images


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

    """
    imagefiltered = ndimage.spline_filter(image.data)

    # finding coefficients of traces which fit the echelle orders across the CCD.
    allcoef, vals, totalnumberoforders = tracesacrossccd(image, imagefiltered, order_of_poly_fits, second_order_coefficient_guess)
    sortedallcoefs = np.array(allcoef)[np.array(allcoef)[:, 0].argsort()]
    order_indices = np.arange(totalnumberoforders)
    # appending indices 0,1,2...,totalnumberoforders as the first column. prior it is indexed from negative numbers.

    coefficients_and_indices = np.insert(sortedallcoefs[:, 1:], obj=0, values=order_indices, axis=1)

    return coefficients_and_indices, vals, totalnumberoforders


def exclude_traces_which_jet_off_detector(coefficients_and_indices, image):
    order_of_poly_fits = coefficients_and_indices.shape[1] - 2
    legendre_polynomial_array, not_needed, not_needed_2 = generate_legendre_array(image, order_of_poly_fits)
    trace_values_versus_xpixel = np.dot(coefficients_and_indices[:, 1:], legendre_polynomial_array)
    # trim any traces which are not contiguously on the detector.
    coefficients_and_indices = coefficients_and_indices[np.all(trace_values_versus_xpixel > 0, axis=1)]
    return coefficients_and_indices


def split_and_sort_coefficients_for_each_fiber(coefficients_and_indices, num_lit_fibers):
    # cutting of order index column
    coefficients_and_indices = coefficients_and_indices[:, 1:]

    # ensuring an even number of traces in the coefficients so that they can be split easily
    while coefficients_and_indices.shape[0] % num_lit_fibers != 0:
        coefficients_and_indices = coefficients_and_indices[:-1]
    num_orders = int(coefficients_and_indices.shape[0] / num_lit_fibers)
    order_indices = [i for i in range(num_orders)]*num_lit_fibers

    fiber_coefficients = (coefficients_and_indices[::num_lit_fibers],)
    for i in range(1, num_lit_fibers):
        fiber_coefficients += (coefficients_and_indices[i::num_lit_fibers],)
    return fiber_coefficients, order_indices


def split_already_sorted_coefficients_into_each_fiber(coefficients_and_indices, num_lit_fibers):
    split_indices = np.where(coefficients_and_indices[:, 0] == np.max(coefficients_and_indices[:, 0]))[0]
    assert len(split_indices) == num_lit_fibers
    fiber_coefficients = (coefficients_and_indices[:split_indices[0] + 1],)
    for i in range(1, len(split_indices)):
        fiber_coefficients += (coefficients_and_indices[split_indices[i-1]+1: split_indices[i]+1],)
    return fiber_coefficients


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
    coeffs_per_fiber, order_indices = split_and_sort_coefficients_for_each_fiber(coefficients_and_indices, num_lit_fibers)
    ordered_coefficients = np.vstack(coeffs_per_fiber)
    coefficients_and_indices = np.insert(ordered_coefficients, obj=0, values=order_indices, axis=1)
    logger.info('%s traces found' % coefficients_and_indices.shape[0])

    return coefficients_and_indices


def optimize_coeffs_entire_lampflat_frame(coefficients_and_indices, image, num_of_lit_fibers=2, order_of_meta_fit=6, bpm=None):
    """
    coefficients which are loaded must be in the format where the first 0,1,2,3,N rows are for the first fiber, and the
    next set for the second fiber. The first column must be the order number designation.
    :param bpm: The binary mask whereby bad pixels are set to 0.
    """

    coeflen, coefwidth = coefficients_and_indices.shape
    number_of_traces, trace_poly_order = coeflen, coefwidth - 2

    if bpm == None:
        bpm = np.ones_like(image.data)

    pixel_x_array = np.arange(image.data.shape[0])
    pixel_y_array = np.arange(image.data.shape[1])

    # generating spline interpolations which incorporate only the good pixels
    f = [interpolate.UnivariateSpline(pixel_y_array[bpm[:, xx] != 0], image.data[:, xx][bpm[:, xx] != 0],
                                      k=3, s=0, ext=1) for xx in pixel_x_array]
    fp = [f[xx].derivative(n=1) for xx in pixel_x_array]
    fpp = [f[xx].derivative(n=2) for xx in pixel_x_array]
    image_splines = ImageSplines(f, fp, fpp)

    legendre_polynomial_array, x, xnorm = generate_legendre_array(image, trace_poly_order)

    order_indices = coefficients_and_indices[:, 0]
    fiber_coefficients = split_already_sorted_coefficients_into_each_fiber(coefficients_and_indices, num_of_lit_fibers)
    optimized_coeffs_per_fiber = []
    for coefficients in fiber_coefficients:
        #  simple order number arrays which are normalized from -1 to 1 to control numerical errors.
        single_fiber_coeffs = coefficients[:, 1:]
        order_norm = np.arange(single_fiber_coeffs.shape[0])
        order_norm = order_norm * 2. / order_norm[-1] - 1

        # initial meta-fit
        meta_coefficients = metacoefficients(order_norm, single_fiber_coeffs, order_of_meta_fit)

        #  building legendre polynomial arrays to use for the meta-fit basis.
        legendre_polynomial_array_meta = np.ones((order_of_meta_fit + 1, order_norm.shape[0]))

        for i in range(1, order_of_meta_fit + 1):
            legendre_polynomial_array_meta[i] = np.polynomial.legendre.legval(order_norm, [0 for j in range(i)] + [1])

        # optimizing the meta-fit
        optimum_meta_coefficients = meta_fit(meta_coefficients, legendre_polynomial_array_meta,
                                                   legendre_polynomial_array, image_splines, x)

        # outputting optimized coeffs with order indices
        optimized_coefficients_no_indices = get_coefficients_from_meta(optimum_meta_coefficients,
                                                                       legendre_polynomial_array_meta)

        optimized_coeffs_per_fiber.append(optimized_coefficients_no_indices)

    all_optimized_coefficients_no_indices = np.vstack(optimized_coeffs_per_fiber)
    all_optimized_coefficients_and_indices = np.insert(all_optimized_coefficients_no_indices, obj=0, values=order_indices,
                                                   axis=1)
    return all_optimized_coefficients_and_indices


def get_number_of_lit_fibers(coefficients_and_indices):
    """
    :param coefficients_and_indices: trace coefficients with indices (0,...,67,0,,..67) as the first column
    :return: the number of lit fibers. This works only if the trace coefficients are all in one array as described.
    This is a temporary function and should be replaced by len(fiber_order) at some point (when we have fiber_order)
    figured out.
    """
    unique, counts = np.unique(coefficients_and_indices[:, 0], return_counts=True)
    number_of_traces = coefficients_and_indices[:, 0].shape[0]
    assert (counts == counts[0]).all()
    num_lit_fibers = int(counts[0])
    assert number_of_traces % num_lit_fibers == 0
    return num_lit_fibers
