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
    :param A: An array of values where you wish to find maxima
    :param s: The window were a points value must be larger than all the other values in the window.
    :param k: multiplicative constant which qualifies a point as a real maxima (not just noise)
    :param ref: the reference value which k*r defines the value any real maximal point must exceed
    :return: The index of the a point near a maximum and the value at the maximum
    """

    i = s - 1
    l, r = 0, 0
    A = A - np.ones_like(A) * np.min(A)
    threshold = k * ref
    firstmax = [-1337, 0]
    while i < len(A) - s and int(l + r) != 2:
        l, r = 1, 1
        for j in range(1, s + 1):
            if A[i] < A[i + j] or A[i] < threshold:
                r, j = 0, s + 1
            if A[i] < A[i - j] or A[i] < threshold:
                l, j = 0, s + 1
        if int(l + r) == 2:
            firstmax = [i, A[i]]

        i += 1
    return firstmax


def crosscoef(p, imfilt, x, arr):
    """
    :param p: list of legendre polynomial coefficients
    :param imfilt: ndattay, image.data passed through ndimage.spline_filter
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param arr: Legendre polynomial array.
    :return: The negative of the total flux summed across the trace.
    """
    y = x * 0.

    for i in range(len(p)):
        y += arr[i] * p[i]

    val = np.sum(ndimage.map_coordinates(imfilt, [y, x], prefilter=False))
    return -1 * val


def fluxvalues(testpoints, p, imfilt, x, arr):
    """
    :param p: list of legendre polynomial coefficients
    :param imfilt: ndattay, image.data passed through ndimage.spline_filter
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param arr: Legendre polynomial array.
    :return: The positive total flux at each point in testpoints
    """
    values = np.zeros_like(testpoints)
    for i in range(0, len(values)):
        coeffguesses = list([testpoints[i]] + p)
        values[i] = (-1) * crosscoef(coeffguesses, imfilt, x, arr)
    return values


def findorder(im, imfilt, x, arr, order=2, lastcoef=None, direction='up'):
    """
    :param im: ndarray, data of the image
    :param imfilt: ndattay, image.data passed through ndimage.spline_filter
    :param x: array of the x pixels from [0,1,2,...,im.shape[1]]
    :param arr: Legendre polynomial array.
    :param order: order of the legendre polynomial fit
    :param lastcoef: [0th_order_coeff, 1st_order_coeff,...] for the last polynomial fit.
    :param direction: Which direction we are proceeding along the detector.
    :return:

    NOTE: This requires that the traces curve upwards on the detector. If for some new detector,
    things reverse, you must change p[2] to -90.
    """

    # Manual guess for the first order to fit--basically just choose a
    # position in the detector and give a good guess for the quadratic
    # coefficient. In this cases 80-100 is good.
    deltap0guess = 0
    if lastcoef is None:
        # bnds = [[None,None]]*(order+1)
        p = [0. for i in range(order + 1)]
        p[0] = int(imfilt.shape[0]/3)
        if order >= 2:
            p[2] = 90
        coeffsguess = copy.deepcopy(p)
    else:
        # Find a new zero-th order term for the polynomial fit by
        # stepping along either up or down the grid until we find
        # a point which is a maximum and nearest to the last order.
        # bnds = [[None,None]]*(order+1)
        p = list(lastcoef[1:])
        p0 = int(lastcoef[0])

        if direction == 'up':
            testpoints = list(range(p0 + 6, p0 + 100))
        elif direction == 'down':
            testpoints = list(range(p0 - 100, p0 - 6))
            testpoints = testpoints[
                         ::-1]  # important that we reverse list to make it greatest to smallest for this down sweep (needed in maxima func)
        elif direction == 'inplace':
            testpoints = list(range(p0 - 10, p0 + 10))

        fluxvals = fluxvalues(testpoints, p, im, x, arr)
        refflux = max((-1)*crosscoef(lastcoef, im, x, arr), max(fluxvals))

        deltap0guess = maxima(fluxvals, 5, 1 / 20, refflux)[0]

        if direction == 'up':
            p0 = deltap0guess + min(testpoints)
        elif direction == 'down':
            p0 = max(testpoints) - deltap0guess
        elif direction == 'inplace':
            p0 = deltap0guess + min(testpoints)

        coeffsguess = [p0] + list(p)

    if deltap0guess != -1337:  # if deltap0guess = -1337, no next order exists.
        p1 = optimize.minimize(crosscoef, coeffsguess, (imfilt, x, arr), method='Powell').x

        val = -1 * crosscoef(p1, imfilt, x, arr)

        return p1, val, 0
    else:
        return lastcoef, refflux, 1


def tracesacrossccd(im, imfilt, order):
    """
    :param im: ndarray, data of the image
    :param imfilt: ndattay, image.data passed through ndimage.spline_filter
    :param order: order of the polynomial fit.
    :return:
    """
    length, width = im.shape

    x = np.arange(imfilt.shape[1])
    xnorm = x * 2. / x[-1] - 1  # x normalized to run from -1 to 1

    # Set up Legendre polynomials to avoid roundoff error from
    # explicitly computing polynomials from their coefficients

    arr = np.ones((order + 1, x.shape[0]))
    for i in range(1, order + 1):
        arr[i] = np.polynomial.legendre.legval(xnorm, [0 for j in range(i)] + [1])

    # For the first coefficients, fit a quadratic
    # unconstrained.  Use this fitted quadratic to fit higher order
    # terms as desired, adding one order at a time to place as few
    # demands as possible on the optimization routine. Future orders are fit with the full
    # nth order polynomial (using the last fit as the initial guess) all at once.
    coef = None
    for i in range(2, order + 1):
        if coef is not None:
            coef = list(coef) + [0]
        coef, val, nomax = findorder(im, imfilt, x, arr, order=i, lastcoef=coef,
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
        coef, val, nomax = findorder(im, imfilt, x, arr, order=order, lastcoef=coef, direction='up')
        vals += [val]
        allcoef += [[i] + list(coef)]
        if coef[0] < allcoef[-2][1] or coef[0] > length and nomax == 0:
            allcoef = allcoef[:-1]  # delete bad fit
            ordersabove = i - 1
            done = 1
        if i >= 2 and done == 0:
            if abs(coef[0] - allcoef[-2][1]) < 1 and abs(coef[0] - allcoef[-3][1]) < 1 and nomax == 0:
                allcoef = allcoef[:-2]  # delete repeated fits
                print(allcoef[-1][1], coef[0])
                ordersabove = i - 2
                done = 1
        if nomax == 1:
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
        coef, val, nomax = findorder(im, imfilt, x, arr, order=order, lastcoef=coef,
                                     direction='down')
        vals += [val]
        allcoef += [[-i] + list(coef)]

        if coef[0] < 0 or coef[0] > allcoef[-2][1] and nomax == 0:
            allcoef = allcoef[:-1]  # delete bad fit
            ordersbelow = i - 1
            done = 1
        if i >= 2 and done == 0 and nomax == 0:
            if abs(coef[0] - allcoef[-2][1]) < 1 and abs(coef[0] - allcoef[-3][1]) < 1:
                allcoef = allcoef[:-2]  # delete repeated fits
                print(allcoef[-1][1], coef[0])
                ordersbelow = i - 2
                done = 1
        if nomax == 1:
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


def recognize_fibers_and_split_coefficients(coefficients, position_zero_of_first_fiber_at_image_center, first_fiber_designation, second_fiber_designation, num_of_orders, image):
    coefficients = coefficients[:, 1:]
    # cutting off order index column.
    assert num_of_orders % 2 == 0
    assert coefficients.shape[0] >= num_of_orders

    order_indices = [i for i in range(int(num_of_orders/2))]*2
    order_of_poly_fits = coefficients.shape[1] - 1

    legendre_polynomial_array, not_needed, not_needed_2 = generate_legendre_array(image, order_of_poly_fits)

    trace_values_versus_xpixel = np.dot(coefficients, legendre_polynomial_array)

    zeroth_order_row = find_nearest(trace_values_versus_xpixel[:, int(image.data.shape[1]/2)], position_zero_of_first_fiber_at_image_center)
    coefficients = coefficients[zeroth_order_row : zeroth_order_row + num_of_orders]


    first_fiber_coefficients = coefficients[::2]
    second_fiber_coefficients = coefficients[1::2]
    fiber_order = tuple((first_fiber_designation, second_fiber_designation))
    ordered_coefficients = np.vstack((first_fiber_coefficients, second_fiber_coefficients))
    coefficients_and_indices = np.insert(ordered_coefficients, obj=0, values=order_indices, axis=1)
    return coefficients_and_indices, fiber_order


def check_for_close_fit(coefficients_and_indices_list, images, max_pixel_error=1E-1):
    """
    :param coefficients_and_indices_list: ndarray with the first column 0,1,2,..66,0,1.. the fiber indexes, and the second column
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
    assert num_traces %2 == 0
    num_orders = int(num_traces/2)
    order_buffer = int(num_orders/5)
    x_buffer = int(image_shape[1]/4)

    # computing absolute differences between the trace centroid locations at every x value.
    trace_values_versus_xpixel_arr = np.array(trace_values_versus_xpixel_list)
    error_between_fits = np.abs((trace_values_versus_xpixel_arr - np.mean(trace_values_versus_xpixel_arr, axis=0)))
    # restricting region where we care about differences (center of detector)
    select_errors_first_fiber = error_between_fits[:,order_buffer:(num_orders-order_buffer), x_buffer:(image_shape[1] - x_buffer)]
    select_errors_second_fiber = error_between_fits[:,num_orders + order_buffer:(num_traces-order_buffer), x_buffer:(image_shape[1] - x_buffer)]

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
    print('(new_flux - master_cal_flux)/master_cal_flux) = %s' % delta_fraction_flux)
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
    legendre_array, xnorm, x = generate_legendre_array(image, order_of_poly_fits)
    X = list(x)*(coefficients_and_indices.shape[0])
    #  X = [0,1,...,4095,0,1,..,4095,..]
    TraceYvals = np.dot(coefficients_and_indices[:,1:],legendre_array).flatten()
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
    assert num_traces % 2 == 0
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
    if len(images) >= cross_correlate_num >= 2:
        image_indices_to_try = list(itertools.combinations(range(len(images)), cross_correlate_num))
        try_combinations_of_images = True
    if len(images) < cross_correlate_num or cross_correlate_num == 1:
        image_indices_to_try = list(range(len(images)))
        try_combinations_of_images = False
    return image_indices_to_try, try_combinations_of_images


def extract_coeffs_entire_lampflat_frame(image, order_of_poly_fits):
    """
    This extracts the trace coefficients for each bright order of a frame. This is only stable for lampflat frames.
    It returns a list of the coefficients, ordered arbitrarily (fibers are not separated). It also returns the summed fluxed across each order
    called val, and the total number of orders found by the algorithm.
    Parameters:
        image = Banzai Image object.
        orderofpolyfits = order of the polynomial fit per trace.
        xvals = list of x coordinates where you want to sample the center of the traces
        debug = will plot sampled points if True

    """

    image_data = image.data
    imagefiltered = ndimage.spline_filter(image_data)

    # finding coefficients of traces which fit the echelle orders across the CCD.
    allcoef, vals, totalnumberoforders = tracesacrossccd(image_data, imagefiltered, order_of_poly_fits)
    sortedallcoefs = np.array(allcoef)[np.array(allcoef)[:, 0].argsort()]
    order_indices = np.arange(totalnumberoforders)
    # appending indices 0,1,2...,totalnumberoforders as the first column. prior it is indexed from negative numbers.

    coefficients_and_indices = np.insert(sortedallcoefs[:, 1:], obj=0, values=order_indices, axis=1)

    return coefficients_and_indices, vals, totalnumberoforders


def fit_traces_order_by_order(image, order_of_poly_fits=4):
    """
    :param image: Banzai image object
    :param order_of_poly_fits: Highest order of the polynomial fit to each trace. 4 is good. Do not change needlessly.
    :internal params: CHANGE WITH CAUTION
            num_of_orders = 67 per fiber, giving 134 traces. Only change if for some reason you want to exclude
                            orders near the bottom of the detector
            position_zero_of_uppermost_fiber_at_image_center = Pixel position of fiber_one near the bottom of the detector
                            called uppermost because python displays CCD images with pixels 0,0 at the bottom.
                            E.g. 47 would be a good one. The correct region is where the orders are most densely packed.
                            The order that you select for this, will become what you designate uppermost_fiber_designation.
                            The fiber right below it will be lowermost_fiber_designation.
                uppermost_fiber_designation = 0
                lowermost_fiber_designation = 1 explained above.
    :return array of trace fit coefficients arranged such that those for uppermost_fiber_designation are first, e.g.
    the first 67 rows of the array.
    """

    coefficients_and_indices, vals, totalnumberoftraces = extract_coeffs_entire_lampflat_frame(image, order_of_poly_fits)
    # do not overscan trim image.
    num_of_orders = 67
    position_zero_of_uppermost_fiber_at_image_center = 47
    uppermost_fiber_designation = 0
    lowermost_fiber_designation = 1
    image_center = int(image.data.shape[1]/2)

    logger.info('%s traces found'%totalnumberoftraces)
    logger.info('selecting only %s of them, with the fiber %s centroid at y=%s x=%s'%(int(num_of_orders*2), uppermost_fiber_designation,
                                                                             position_zero_of_uppermost_fiber_at_image_center, image_center))

    coefficients_and_indices, fiber_order = recognize_fibers_and_split_coefficients(coefficients_and_indices,
                                            position_zero_of_uppermost_fiber_at_image_center, uppermost_fiber_designation,
                                            lowermost_fiber_designation, int(num_of_orders*2), image)
    logger.info(str(fiber_order))
    logger.info(coefficients_and_indices.shape)
    logger.info('inside fit_traces_order_by_order')
    return coefficients_and_indices, fiber_order


def optimize_coeffs_entire_lampflat_frame(trace_coefficients, image, order_of_meta_fit=6, bpm=None):
    """
    coefficients which are loaded must be in the format where the first 0,1,2,3,N rows are for the first fiber, and the
    next set for the second fiber. The first column must be the order number designation.
    :param bad_pixel_mask: The binary mask whereby bad pixels are set to 0.
    """

    coeflen, coefwidth = trace_coefficients.shape
    number_of_traces, trace_poly_order = coeflen, coefwidth - 2

    assert number_of_traces % 2 == 0

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

    # requires input coefficients to be ordered correctly.
    first_coefficients = trace_coefficients[:int(number_of_traces / 2), 1:]
    second_coefficients = trace_coefficients[int(number_of_traces / 2):, 1:]
    #  simple order number arrays which are normalized from -1 to 1 to control numerical errors.
    first_norm = np.arange(first_coefficients.shape[0])
    second_norm = np.arange(second_coefficients.shape[0])
    first_norm = first_norm * 2. / first_norm[-1] - 1
    second_norm = second_norm * 2. / second_norm[-1] - 1

    # initial meta-fit
    first_coefficients_meta = metacoefficients(first_norm, first_coefficients, order_of_meta_fit)
    second_coefficients_meta = metacoefficients(second_norm, second_coefficients, order_of_meta_fit)

    #  building legendre polynomial arrays to use for the meta-fit basis.
    legendre_polynomial_array_meta_first = np.ones((order_of_meta_fit + 1, first_norm.shape[0]))
    legendre_polynomial_array_meta_second = np.ones((order_of_meta_fit + 1, second_norm.shape[0]))

    for i in range(1, order_of_meta_fit + 1):
        legendre_polynomial_array_meta_first[i] = np.polynomial.legendre.legval(first_norm, [0 for j in range(i)] + [1])
        legendre_polynomial_array_meta_second[i] = np.polynomial.legendre.legval(second_norm,
                                                                                 [0 for j in range(i)] + [1])

    # optimizing the meta-fit
    first_coefficients_meta_optimum = meta_fit(first_coefficients_meta, legendre_polynomial_array_meta_first,
                                               legendre_polynomial_array, image_splines, x)
    second_coefficients_meta_optimum = meta_fit(second_coefficients_meta, legendre_polynomial_array_meta_second,
                                                legendre_polynomial_array, image_splines, x)

    # outputting optimized coeffs with order indices
    first_coefficients_optimum = get_coefficients_from_meta(first_coefficients_meta_optimum,
                                                            legendre_polynomial_array_meta_first)

    second_coefficients_optimum = get_coefficients_from_meta(second_coefficients_meta_optimum,
                                                             legendre_polynomial_array_meta_second)

    optimized_coefficients = np.vstack((first_coefficients_optimum, second_coefficients_optimum))

    order_indices = trace_coefficients[:, 0]
    optimized_coefficients_and_indices = np.insert(optimized_coefficients, obj=0, values=order_indices, axis=1)

    return optimized_coefficients_and_indices
