"""
trace_utils.py: Routines for finding echelle orders across a CCD.
Authors
    G. Mirek Brandt (gmbrandt@ucsb.edu)

    Tim Brandt (tbrandt@physics.ucsb.edu)
"""

import numpy as np
from scipy import ndimage, optimize, signal
from astropy.table import Table
from astropy.io import fits

import logging

logger = logging.getLogger(__name__)


class Trace(object):
    """
    :param data = {'id': ndarray, 'centers': ndarray}. 'centers' gives a 2d array, where
    the jth row are the y centers across the detector for the trace with identification trace_centers['id'][j]
    :param design_matrix = the 2d array such that coefficients dot design_matrix gives the trace centers for all the orders.
    """
    def __init__(self, data=None, trace_table_name=None):
        if data is None:
            data = {'id': [], 'centers': []}
        self.trace_table_name = trace_table_name
        self.data = Table(data)
        self.data['id'].description = 'Identification tag for trace'
        self.data['centers'].description = 'Vertical position of the center of the trace as a function of horizontal pixel'
        self.data['centers'].unit = 'pixel'

    def get_centers(self, row):
        return self.data['centers'][row]

    def add_centers(self, trace_centers, trace_id):
        # TODO fix the fact that astropy cannot add a row to an empty table
        if len(self.data['id']) == 0:
            self.data = Trace(data={'id': [trace_id], 'centers': [trace_centers]}).data
        else:
            self.data.add_row([trace_id, trace_centers])

    def _del_last_fit(self):
        self.data.remove_row(-1)

    def num_traces_found(self):
        return len(self.data['id'])

    @staticmethod
    def fit_traces(image, poly_fit_order, second_order_coefficient_guess,
                   fit_march_parameters=None, match_filter_parameters=None):
        start_point = image.data.shape[0]/3
        trace = Trace(data=None)
        trace_fitter = SingleTraceFitter(image_data=image.data, start_point=start_point,
                                         second_order_coefficient_guess=second_order_coefficient_guess,
                                         poly_fit_order=poly_fit_order,
                                         march_parameters=fit_march_parameters,
                                         match_filter_parameters=match_filter_parameters)
        at_edge = False
        trace_id = 0
        trace_fitter.generate_initial_guess()
        for direction in ['up', 'down']:
            if direction == 'down':
                trace_fitter.use_very_first_fit_as_initial_guess()
                at_edge = trace_fitter.match_filter_to_refine_initial_guess(trace.get_centers(0), direction=direction)
            while not at_edge:
                trace_centers = trace_fitter.fit_trace()
                trace.add_centers(trace_centers, trace_id)
                trace_id += 1

                trace_fitter.use_previous_fit_as_initial_guess()
                at_edge = trace_fitter.match_filter_to_refine_initial_guess(trace_centers, direction=direction)
                if trace._bad_fit(image.data, direction):
                    trace._del_last_fit()
                    at_edge = True
            trace._sort_traces()
        return trace

    def write(self, filename):
        hdu = fits.BinTableHDU(self.data, name=self.trace_table_name)
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(filename)

    @staticmethod
    def load(hdu_list, trace_table_name):
        return Trace(data=hdu_list[trace_table_name].data)

    def _repeated_fit(self):
        center_pixel = int(len(self.get_centers(0)) / 2)
        a_repeated_fit = False
        if self.num_traces_found() < 2:
            a_repeated_fit = False
        elif np.isclose(self.get_centers(-1)[center_pixel], self.get_centers(-2)[center_pixel], atol=2, rtol=0):
            a_repeated_fit = True
        return a_repeated_fit

    def _bad_fit(self, image_data, direction='up'):
        return any((self._bad_shift(direction), self._beyond_edge(image_data), self._repeated_fit()))

    def _bad_shift(self, direction='up'):
        center_pixel = int(len(self.get_centers(0)) / 2)
        a_bad_shift = False
        if self.num_traces_found() < 2:
            a_bad_shift = False
        elif direction == 'up' and self.get_centers(-1)[center_pixel] < self.get_centers(-2)[center_pixel]:
            a_bad_shift = True
        elif direction == 'down' and self.get_centers(-1)[center_pixel] > self.get_centers(-2)[center_pixel]:
            a_bad_shift = True
        return a_bad_shift

    def _beyond_edge(self, image_data):
        """
        :param image_data:
        :return: True or False whether the y value at the center of the most recent trace fit is <0 or greater than
        the maximum y coordinate of the image (i.e. image_data.shape[0])
        """
        center_pixel = int(len(self.get_centers(0)) / 2)
        return any((self.get_centers(-1)[center_pixel] < 0, self.get_centers(-1)[center_pixel] > image_data.shape[0]))

    def _sort_traces(self):
        center = int(self.data['centers'].shape[1] / 2)
        self.data['centers'] = np.array(self.data['centers'])
        self.data['centers'] = self.data['centers'][self.data['centers'][:, center].argsort()]
        self.data['id'] = np.arange(self.data['centers'].shape[0])


class SingleTraceFitter(object):
    def __init__(self, image_data=None, poly_fit_order=2, start_point=None, second_order_coefficient_guess=None,
                 march_parameters=None, match_filter_parameters=None, extraargs=None):
        if extraargs is None:
            extraargs = {}
        if march_parameters is None:
            march_parameters = {'window': 100, 'step_size': 6}
        if match_filter_parameters is None:
            match_filter_parameters = {'min_peak_spacing': 5, 'neighboring_peak_flux_ratio': 5}
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
        self.match_filter_parameters = match_filter_parameters
        if extraargs.get('initialize_fit_objects', True):
            self._initialize_fit_objects()

    def fit_trace(self):
        initial_guess = self.initial_guess_next_fit
        refined_coefficients = optimize.minimize(SingleTraceFitter._trace_merit_function, initial_guess,
                                                 args=(self,), method='Powell').x
        self.coefficients.append(refined_coefficients)
        return self._centers_from_coefficients(refined_coefficients)

    def match_filter_to_refine_initial_guess(self, current_trace_centers, direction='up'):
        no_more_traces = False
        shifted_trace_centers, offsets = self._centers_shifting_traces_up_or_down(current_trace_centers,
                                                                                  direction=direction)
        flux_vs_shift = self._flux_as_trace_shifts_up_or_down(shifted_trace_centers)
        current_trace_flux = self._flux_across_trace(current_trace_centers)
        reference_flux = max(current_trace_flux, np.max(flux_vs_shift))
        min_peak_height = abs(reference_flux)/self.match_filter_parameters['neighboring_peak_flux_ratio']
        flux_vs_shift[flux_vs_shift < min_peak_height] = min_peak_height
        peak_indices = signal.find_peaks(flux_vs_shift,
                                         distance=self.match_filter_parameters['min_peak_spacing'])[0]
        if len(peak_indices) == 0:
            peak_indices = [0]
            no_more_traces = True
        index_of_first_peak = peak_indices[0]
        self.initial_guess_next_fit[0] += offsets[index_of_first_peak]
        return no_more_traces

    def _flux_as_trace_shifts_up_or_down(self, shifted_traces):
        flux_vs_shift = []
        for trace in shifted_traces:
            flux_vs_shift.append(self._flux_across_trace(trace))
        return np.array(flux_vs_shift)

    def _centers_shifting_traces_up_or_down(self, current_trace_centers, direction='up'):
        window = self.march_parameters['window']
        step = self.march_parameters['step_size']
        if direction == 'up':
            offsets = np.arange(step, window+step, 1)
        if direction == 'down':
            offsets = np.arange((-1)*step, (-1)*(window+step), -1)
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

    @staticmethod
    def _generate_design_matrix(normalized_domain, poly_fit_order):
        design_matrix = np.ones((poly_fit_order + 1, normalized_domain.shape[0]))
        for i in range(poly_fit_order + 1):
            design_matrix[i] = legendre(order=i, values=normalized_domain)
        return design_matrix

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


def legendre(order, values):
    #TODO replace with scipy legendre(order)(values) if it works the same.
    if order == 0:
        return np.ones_like(values)
    else:
        return np.polynomial.legendre.legval(values, [0 for j in range(order)] + [1])
