import pytest
import numpy as np
from astropy.table import Table
from unittest import mock
import logging

from banzai_nres.tests.utils import array_with_two_peaks
from banzai_nres.utils.new_trace_utils import Trace, SingleTraceFitter

logger = logging.getLogger(__name__)


class TestTrace:
    """
    Unit tests for the Trace() class.
    """
    def test_class_attributes(self):
        trace = Trace()
        assert trace.trace_table_name is None
        assert trace.data.colnames == ['id', 'centers']
        assert len(trace.data['id']) == 0
        assert len(trace.data['centers']) == 0

    def test_getting_trace_centers(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert np.allclose(trace.get_centers(0), [0, 1])

    def test_getting_num_found_traces(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert trace.num_traces_found() == 2

    def test_add_centers(self):
        trace = Trace(data=None)
        trace.add_centers(trace_centers=np.array([1, 2, 3]), trace_id=1)
        assert np.allclose(trace.data['centers'], [[1, 2, 3]])
        assert np.allclose(trace.data['id'], [1])
        trace = Trace(data={'id': [1], 'centers': [[1, 2, 3]]})
        trace.add_centers(trace_centers=np.array([1, 2, 3]), trace_id=2)
        assert np.allclose(trace.data['centers'], [[1, 2, 3], [1, 2, 3]])
        assert np.allclose(trace.data['id'], [1, 2])

    def test_write(self):
        #TODO
        assert True

    @mock.patch('banzai_nres.utils.new_trace_utils.Trace._bad_shift')
    @mock.patch('banzai_nres.utils.new_trace_utils.Trace._repeated_fit')
    @mock.patch('banzai_nres.utils.new_trace_utils.Trace._beyond_edge')
    def test_bad_fit(self, beyond_edge, repeated_fit, bad_shift):
        beyond_edge.return_value, repeated_fit.return_value, bad_shift.return_value = True, True, True
        assert Trace()._bad_fit(image_data=None, direction=None)
        beyond_edge.return_value, repeated_fit.return_value, bad_shift.return_value = False, False, False
        assert not Trace()._bad_fit(image_data=None, direction=None)
        beyond_edge.return_value, repeated_fit.return_value, bad_shift.return_value = True, False, True
        assert Trace()._bad_fit(image_data=None, direction=None)

    def test_detecting_repeated_fit(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2], 'centers': [centers, centers+0.1]}
        trace = Trace(data=data)
        assert trace._repeated_fit()

    def test_detecting_bad_shifts(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1], 'centers': np.array([centers])}
        trace = Trace(data=data)
        assert not trace._bad_shift(direction='up')

        trace.data = {'id': [1, 2], 'centers': np.array([centers, centers-0.1])}
        assert trace._bad_shift(direction='up')

        trace.data = {'id': [1, 2], 'centers': np.array([centers, centers+0.1])}
        assert trace._bad_shift(direction='down')

    def test_detecting_beyond_edge(self):
        fake_image_data = np.zeros((3, 3))
        data = {'id': [1], 'centers': [[3.1, 3.1, 3.1]]}
        assert Trace(data=data)._beyond_edge(fake_image_data)
        data = {'id': [1], 'centers': [[-0.1, -0.1, -0.1]]}
        assert Trace(data=data)._beyond_edge(fake_image_data)
        data = {'id': [1], 'centers': [[2, -0.1, 2]]}
        assert Trace(data=data)._beyond_edge(fake_image_data)

    def test_sorting_trace_centers(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers-10, centers+2]}
        trace = Trace(data=data)
        trace._sort_traces()
        assert np.allclose(trace.data['id'], np.arange(4))
        assert np.allclose(trace.data['centers'],
                           np.array([centers-10, centers, centers+2, centers+5]))

    def test_delete_centers(self):
        centers = np.array([1, 2, 3, 4])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers+10, centers+11]}
        trace = Trace(data=data)
        trace._del_centers([])
        assert np.allclose(trace.data['id'], data['id'])
        assert np.allclose(trace.data['centers'], data['centers'])
        trace._del_centers(-1)
        assert np.allclose(trace.data['id'], [1, 2, 3])
        assert np.allclose(trace.data['centers'], [centers, centers+5, centers+10])
        trace._del_centers([-1, -2])
        assert np.allclose(trace.data['id'], [1])
        assert np.allclose(trace.data['centers'], [centers])


class TestSingleTraceFitter:
    def test_class_attributes(self):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        assert fitter.second_order_coefficient_guess is None
        assert fitter.start_point is None
        assert fitter.image_data is None
        assert fitter.filtered_image_data is None
        assert fitter.initial_guess_next_fit is None
        assert fitter.x is None
        assert fitter.x_norm is None
        assert fitter.design_matrix is None

    def test_default_class_attributes(self):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        assert fitter.march_parameters['window'] == 100
        assert fitter.march_parameters['step_size'] == 6
        assert fitter.poly_fit_order == 2
        assert fitter.coefficients == []

    def test_fit_initilization(self):
        """
        tests that SingleTraceFitter calls _initialize_fit_objects correctly upon
        instantiation of the class.
        """
        poly_fit_order = 2
        fitter = SingleTraceFitter(image_data=np.zeros((2, 2)),
                                   poly_fit_order=poly_fit_order,
                                   start_point=1,
                                   second_order_coefficient_guess=90)
        assert np.allclose(fitter.x_norm, np.array([-1, 1]))
        assert np.allclose(fitter.x, np.arange(2))
        assert fitter.filtered_image_data is not None
        assert fitter.design_matrix.shape == (poly_fit_order+1, 2)

    def test_generating_initial_guess(self):
        fitter = SingleTraceFitter(image_data=np.zeros((2, 2)),
                                   poly_fit_order=2,
                                   start_point=1,
                                   second_order_coefficient_guess=90)
        fitter.generate_initial_guess()
        assert np.allclose(fitter.initial_guess_next_fit, np.array([1, 0, 90]))

    def test_changing_initial_guesses(self):
        coefficients = np.array([[0, 0], [1, 1]])
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'coefficients': coefficients})
        fitter.use_previous_fit_as_initial_guess()
        assert np.allclose(fitter.initial_guess_next_fit, coefficients[-1])
        fitter.use_very_first_fit_as_initial_guess()
        assert np.allclose(fitter.initial_guess_next_fit, coefficients[0])

    def test_generating_initial_guess_fail(self):
        fitter = SingleTraceFitter(image_data=np.zeros((2, 2)),
                                   poly_fit_order=2,
                                   start_point=None,
                                   second_order_coefficient_guess=90)
        fitter.generate_initial_guess()
        assert fitter.initial_guess_next_fit is None

    def test_centers_from_coefficients(self):
        design_matrix = np.ones((2, 5))
        design_matrix[1] = np.linspace(-1, 1, 5)
        offset, linear_coefficient = 1, 0.5
        x = np.arange(design_matrix.shape[1])
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'design_matrix': design_matrix,
                                              'x': x})
        single_trace_coefficients = np.array([offset, linear_coefficient])
        trace_centers = fitter._centers_from_coefficients(single_trace_coefficients)
        a_line = np.linspace(offset - linear_coefficient, offset + linear_coefficient, 5)
        assert np.allclose(trace_centers, a_line)

    def test_flux_across_trace(self):
        x = np.arange(5)
        fake_data = np.zeros((9, len(x)))
        fake_data[3] += 1
        trace_centers = np.ones(len(x))*3
        filtered_fake_data = SingleTraceFitter._prefilter_image_data(fake_data)
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'x': x,
                                              'filtered_image_data': filtered_fake_data})
        expected_flux = len(x)
        flux = fitter._flux_across_trace(trace_centers)
        assert np.isclose(flux, expected_flux)

    def test_trace_merit_function_returns_negative_flux(self):
        x = np.arange(5)
        fake_data = np.zeros((9, len(x)))
        fake_data[3] += 1
        coefficients_for_line = np.array([3, 0])
        design_matrix = np.array([np.ones(len(x)),
                                  np.linspace(-1, 1, len(x))])
        filtered_fake_data = SingleTraceFitter._prefilter_image_data(fake_data)
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'x': x,
                                              'design_matrix': design_matrix,
                                              'filtered_image_data': filtered_fake_data})
        negative_flux = (-1)*len(x)
        merit = fitter._trace_merit_function(single_trace_coefficients=coefficients_for_line,
                                             cls=fitter)
        assert np.isclose(merit, negative_flux)

    def test_normalizing_coordinates(self):
        x = np.arange(5)
        x_norm = np.linspace(-1, 1, len(x))
        fake_data = np.zeros((9, len(x)))
        fitter = SingleTraceFitter(image_data=fake_data, extraargs={'initialize_fit_objects': False})
        fitter._normalize_domain_coordinates()
        assert np.allclose(fitter.x, x)
        assert np.allclose(fitter.x_norm, x_norm)

    def test_generating_legendre_design_matrix(self):
        x = np.arange(5)
        xnorm = np.linspace(-1, 1, len(x))
        design_matrix = np.array([np.ones(len(x)),
                                  xnorm])
        assert np.allclose(design_matrix,
                           SingleTraceFitter._generate_design_matrix(xnorm, poly_fit_order=1))

    def test_fit_trace(self):
        #TODO
        assert True


class TestMatchFilter:
    """
    Tests for the functions in single trace fitter which are used to find
    the approximate locations of the next trace during a march up the detector.
    """
    def test_centers_shifting_traces_up_or_down(self):
        w, s = 2, 6
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False},
                                   march_parameters={'window': w, 'step_size': s})
        current_trace_centers = np.arange(10)
        shifted_trace_centers, offsets = fitter._centers_shifting_traces_up_or_down(current_trace_centers,
                                                                                    direction='up')
        assert np.allclose(shifted_trace_centers, np.array([current_trace_centers + s,
                                                            current_trace_centers + s + 1]))
        assert np.allclose(offsets, [6, 7])

        shifted_trace_centers, offsets = fitter._centers_shifting_traces_up_or_down(current_trace_centers,
                                                                                    direction='down')
        assert np.allclose(shifted_trace_centers, np.array([current_trace_centers - s,
                                                            current_trace_centers - s - 1]))
        assert np.allclose(offsets, [-6, -7])

    @mock.patch('banzai_nres.utils.new_trace_utils.SingleTraceFitter._flux_across_trace', side_effect=np.max)
    def test_getting_flux_as_trace_shifts_up_or_down(self, flux_across_trace):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        shifted_traces = np.array([[1, 2], [2, 3], [4, 5]])
        flux_per_trace = fitter._flux_as_trace_shifts_up_or_down(shifted_traces)
        assert np.allclose(flux_per_trace, np.max(shifted_traces, axis=1))

    @mock.patch('banzai_nres.utils.new_trace_utils.SingleTraceFitter._flux_across_trace')
    @mock.patch('banzai_nres.utils.new_trace_utils.SingleTraceFitter._flux_as_trace_shifts_up_or_down')
    @mock.patch('banzai_nres.utils.new_trace_utils.SingleTraceFitter._centers_shifting_traces_up_or_down')
    def test_match_filter_to_refine_initial_guess(self, shift_centers, flux_vs_shift, reference_flux):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        fitter.match_filter_parameters = {'min_peak_spacing': 5, 'neighboring_peak_flux_ratio': 20}
        positive_trace_signal, centroids, x_coords = array_with_two_peaks()
        no_trace_signal = np.random.normal(loc=0, scale=np.max(positive_trace_signal)/100,
                                           size=len(positive_trace_signal))
        shift_centers.return_value = None, x_coords
        reference_flux.return_value = np.max(positive_trace_signal)
        for trace_signal, outcome, prediction in zip([positive_trace_signal, no_trace_signal],
                                                     [centroids[0], 0],
                                                     [False, True]):
            fitter.initial_guess_next_fit = [0]
            flux_vs_shift.return_value = trace_signal
            no_more_traces = fitter.match_filter_to_refine_initial_guess(current_trace_centers=None, direction=None)
            assert np.isclose(fitter.initial_guess_next_fit[0], outcome, atol=3, rtol=0)
            assert no_more_traces is prediction


