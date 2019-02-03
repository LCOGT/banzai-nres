import pytest
import numpy as np
from unittest import mock
import logging

from banzai_nres.utils.new_trace_utils import Trace, SingleTraceFitter

logger = logging.getLogger(__name__)


class TestTrace:
    """
    Unit tests for the Trace() class.
    """
    def test_class_attributes(self):
        trace = Trace()
        assert trace.trace_table_name is None
        assert trace.data is None
        assert trace.second_order_coefficient_guess is None
        assert trace.poly_fit_order is None

    def test_getting_trace_centers(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert trace.get_trace_centers(0) == [0, 1]

    def test_write(self):
        #TODO
        assert True

    def test_detecting_repeated_fit(self):
        fake_image_data = np.zeros((3, 3))
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2], 'centers': [centers, centers+0.1]}
        trace = Trace(data=data)
        assert trace._repeated_fit(fake_image_data)

    def test_detecting_bad_fits(self):
        fake_image_data = np.zeros((3, 3))
        centers = np.array([1, 2, 3])
        data = {'id': [1], 'centers': [centers]}
        trace = Trace(data=data)
        assert not trace._bad_fit(fake_image_data, direction='up')

        trace.data = {'id': [1, 2], 'centers': [centers, centers-0.1]}
        assert trace._bad_fit(fake_image_data, direction='up')

        trace.data = {'id': [1, 2], 'centers': [centers, centers+0.1]}
        assert trace._bad_fit(fake_image_data, direction='down')

    def test_detecting_beyond_edge(self):
        fake_image_data = np.zeros((3, 3))
        single_trace_centers = [3.1, 3.1, 3.1]
        assert Trace._beyond_edge(single_trace_centers, fake_image_data)
        single_trace_centers = [-0.1, -0.1, -0.1]
        assert Trace._beyond_edge(single_trace_centers, fake_image_data)

    def test_sorting_trace_centers(self):
        fake_image_data = np.zeros((3, 3))
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers-10, centers+2]}
        trace = Trace(data=data)
        trace._sort_traces_hierarchically(fake_image_data)
        assert np.allclose(trace.data['id'], np.arange(4))
        assert np.allclose(trace.data['centers'],
                           np.array([centers-10, centers, centers+2, centers+5]))

    def test_trimming_bad_fits(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2, 3],
                'centers': [centers, centers+5, centers+10]}
        trace = Trace(data=data)
        trace._trim_last_fit_based_on_criterion(True)
        assert np.allclose(trace.data['id'], [1, 2])
        assert np.allclose(trace.data['centers'], [centers, centers+5])

        trace._trim_last_fit_based_on_criterion(False)
        assert np.allclose(trace.data['id'], [1, 2])
        assert np.allclose(trace.data['centers'], [centers, centers+5])


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
        coefficients = np.array([[0, 0],[1, 1]])
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

    def test_trace_centers_from_coefficients(self):
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

    def test_flux_across_image(self):
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

    def test_normalizing_coordinates(self):
        x = np.arange(5)
        x_norm = np.linspace(-1, 1, len(x))
        fake_data = np.zeros((9, len(x)))
        fitter = SingleTraceFitter(image_data=fake_data, extraargs={'initialize_fit_objects': False})
        fitter._normalize_domain_coordinates()
        assert np.allclose(fitter.x, x)
        assert np.allclose(fitter.x_norm, x_norm)
