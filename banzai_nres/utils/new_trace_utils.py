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
    :param data = {'id': ndarray, 'centers': ndarray}. 'centers' gives a 2d array, where
    the jth row are the y centers across the detector for the trace with identification trace_centers['id'][j]
    :param design_matrix = the 2d array such that coefficients dot design_matrix gives the trace centers for all the orders.
    """
    def __init__(self, data=None, trace_table_name=None, second_order_coefficient_guess=None, poly_fit_order=None):
        self.trace_table_name = trace_table_name # will be taken from settings. when instantiated in TraceMaker
        self.data = data
        self.second_order_coefficient_guess = second_order_coefficient_guess
        self.poly_fit_order = poly_fit_order

    def get_trace_centers(self, row):
        return self.data['centers'][row]

    @staticmethod
    def fit_traces(image, poly_fit_order, second_order_coefficient_guess):
        trace = Trace(data={'id': [], 'centers': []},
                      poly_fit_order=poly_fit_order,
                      second_order_coefficient_guess=second_order_coefficient_guess)
        start_point = image.data.shape[0]/3
        trace_fitter = SingleTraceFitter(image_data=image.data, start_point=start_point,
                                         second_order_coefficient_guess=second_order_coefficient_guess,
                                         poly_fit_order=poly_fit_order)
        trace_fitter.generate_initial_guess()
        at_edge = False
        direction = 'up'
        trace_id = 0
        while not at_edge:
            single_trace_centers = trace_fitter.fit_trace()
            trace_fitter.use_previous_fit_as_initial_guess()
            no_more_traces = trace_fitter.match_filter_to_refine_initial_guess(single_trace_centers,
                                                                               direction=direction)
            trace.data['centers'].append(single_trace_centers)
            trace.data['id'].append(trace_id)
            trace_id += 1

            beyond_edge = Trace._beyond_edge(single_trace_centers, image_data=image.data)
            a_repeated_fit = trace._repeated_fit(image.data)
            a_bad_fit = trace._bad_fit(image.data, direction=direction)
            if no_more_traces or beyond_edge or a_repeated_fit or a_bad_fit:
                if beyond_edge or a_repeated_fit or a_bad_fit:
                    del trace.data['centers'][-1]
                    del trace.data['id'][-1]
                if direction == 'up':
                    direction = 'down'
                    trace_fitter.use_very_first_fit_as_initial_guess()
                    at_edge = trace_fitter.match_filter_to_refine_initial_guess(trace.data['centers'][0],
                                                                                direction=direction)
                else:
                    at_edge = True
        trace._sort_traces_hierarchically(image.data)
        return trace

    def write(self, filename):
        table = Table(self.data, name=self.trace_table_name)
        table['id'].description = 'Blah'
        table.write(filename)

    def _repeated_fit(self, image_data):
        center = int(image_data.shape[1] / 2)
        a_repeated_fit = False
        if len(self.data['id']) < 2:
            a_repeated_fit = False
        elif np.isclose(self.data['centers'][-1][center], self.data['centers'][-2][center], atol=2, rtol=0):
            a_repeated_fit = True
        return a_repeated_fit

    def _bad_fit(self, image_data, direction='up'):
        center = int(image_data.shape[1] / 2)
        a_bad_fit = False
        if len(self.data['id']) < 2:
            a_bad_fit = False
        elif direction == 'up' and self.data['centers'][-1][center] < self.data['centers'][-2][center]:
            a_bad_fit = True
        elif direction == 'down' and self.data['centers'][-1][center] > self.data['centers'][-2][center]:
            a_bad_fit = True
        return a_bad_fit

    @staticmethod
    def _beyond_edge(data, image_data):
        trace_protrudes_off_detector = False
        if np.max(data) > image_data.shape[0]:
            trace_protrudes_off_detector = True
        if np.min(data) < 0:
            trace_protrudes_off_detector = True
        return trace_protrudes_off_detector

    def _sort_traces_hierarchically(self, image_data):
        center = int(image_data.shape[1] / 2)
        self.data['centers'] = np.array(self.data['centers'])
        self.data['centers'] = self.data['centers'][self.data['centers'][:, center].argsort()]
        self.data['id'] = np.arange(self.data['centers'].shape[0])




class SingleTraceFitter(object):
    #TODO add the other non necessary arguments as keyword arguments for extraargs.
    def __init__(self, image_data=None, poly_fit_order=2, start_point=None, second_order_coefficient_guess=None,
                 march_parameters=None, extraargs={}):
        if march_parameters is None:
            march_parameters = {'window': 100, 'step_size': 6}
        if extraargs.get('coefficients') is None:
            extraargs['coefficients'] = []
        self.second_order_coefficient_guess = second_order_coefficient_guess
        self.start_point = start_point
        self.image_data = image_data
        self.poly_fit_order = poly_fit_order
        self.filtered_image_data = extraargs.get('filtered_image_data')
        self.initial_guess_next_fit = extraargs.get('initial_guess_next_fit')
        self.coefficients = extraargs.get('coefficients')
        self.x = extraargs.get('x')
        self.x_norm = extraargs.get('xnorm')
        self.design_matrix = extraargs.get('design_matrix')
        self.march_parameters = march_parameters
        if extraargs.get('initialize_fit_objects', True) is True:
            self._initialize_fit_objects()

    def fit_trace(self):
        initial_guess = self.initial_guess_next_fit
        refined_coefficients = optimize.minimize(SingleTraceFitter._trace_merit_function, initial_guess,
                                                 args=(self,), method='Powell').x
        self.coefficients.append(refined_coefficients)
        return self._centers_from_coefficients(refined_coefficients)

    def match_filter_to_refine_initial_guess(self, current_trace_centers, direction='up'):
        shifted_trace_centers, offsets = self._centers_shifting_traces_up_or_down(current_trace_centers,
                                                                                  direction=direction)
        flux_vs_shift = self._flux_as_trace_shifts_up_or_down(shifted_trace_centers)
        reference_flux = max(self._flux_across_trace(current_trace_centers), np.max(flux_vs_shift))
        index_of_max_flux, maximum_exists = maxima(flux_vs_shift, 5, 1 / 20, reference_flux)
        self.initial_guess_next_fit[0] += offsets[index_of_max_flux]
        no_more_traces = not maximum_exists
        return no_more_traces

    def _flux_as_trace_shifts_up_or_down(self, shifted_traces):
        flux_vs_shift = []
        for trace in shifted_traces:
            flux_vs_shift.append(self._flux_across_trace(trace))
        return flux_vs_shift

    def _centers_shifting_traces_up_or_down(self, current_trace_centers, direction='up'):
        window = self.march_parameters['window']
        step = self.march_parameters['step_size']
        if direction == 'up':
            offsets = np.arange(step, window, 1)
        if direction == 'down':
            offsets = np.arange((-1)*window, (-1)*step, 1)[::-1]
        shifted_trace_centers = np.ones((offsets.shape[0], current_trace_centers.shape[0]))
        shifted_trace_centers *= current_trace_centers
        shifted_trace_centers += np.array([offsets]).T
        return shifted_trace_centers, offsets

    def use_very_first_fit_as_initial_guess(self):
        self.initial_guess_next_fit = self.coefficients[0]

    def use_previous_fit_as_initial_guess(self):
        self.initial_guess_next_fit = self.coefficients[-1]

    def generate_initial_guess(self):
        if self.start_point is None or self.second_order_coefficient_guess is None:
            logger.error('Starting y position up the detector nor the second order guess have been specified '
                         ', failed to generate an initial guess for trace fitting.')
        else:
            self.initial_guess_next_fit = np.zeros(self.poly_fit_order+1).astype(np.float64)
            self.initial_guess_next_fit[0] = self.start_point
            self.initial_guess_next_fit[2] = self.second_order_coefficient_guess

    def _centers_from_coefficients(self, coefficients):
        """
        :param coefficients: a 1d array for the coefficients of a single trace in order of: [0th, 1st,.. mth degree coefficient]
        or a 2d array of shape (N, m+1) for N traces.
        :return: If coefficients is 1d of size m+1, this will return an array of shape image_data.shape[1] with the trace
        centers at each of those points. If coefficients is for many traces (shape (N, m+1)) then this will return a 2d
        array of shape (N, image_data.shape[1]) e.g. the trace centers for each of the traces.
        """
        trace_centers = np.dot(coefficients, self.design_matrix)
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
    def _trace_merit_function(single_trace_coefficients, cls):
        trace_centers = cls._centers_from_coefficients(single_trace_coefficients)
        flux = cls._flux_across_trace(trace_centers)
        return (-1)*flux

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



