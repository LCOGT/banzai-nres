import pytest
import tempfile
import numpy as np
import os

from astropy.table import Table
from astropy.io import fits
from scipy.optimize import OptimizeResult
from unittest import mock

from banzai.tests.utils import FakeContext

from banzai_nres.traces import TraceMaker, LoadTrace
from banzai_nres.tests.utils import array_with_two_peaks, FakeImage, noisify_image
from banzai_nres.utils.trace_utils import Trace, SingleTraceFitter, AllTraceFitter
from banzai_nres.tests.utils import fill_image_with_traces
import banzai_nres.settings

import logging


logger = logging.getLogger(__name__)


class FakeTraceImage(FakeImage):
    """
    Image of 500x500 is recommended. Drastic changes to that dimension may break trace testing.
    """
    def __init__(self, nx=500, ny=502, *args, **kwargs):
        super(FakeTraceImage, self).__init__(*args, **kwargs)
        self.caltype = 'TRACE'
        self.header = fits.Header()
        self.header['OBSTYPE'] = 'LAMPFLAT'
        self.header['OBJECTS'] = 'tung&tung&none'
        self.nx = nx
        self.ny = ny
        self.bpm = np.zeros((self.ny, self.nx), dtype=np.uint8)
        self.data = np.zeros((self.ny, self.nx))
        self.fiber0_lit, self.fiber1_lit, self.fiber2_lit = False, True, True


class TestTrace:
    """
    Unit tests for the Trace class.
    """
    def test_trace_instantiates_from_num_centers(self):
        trace = Trace(num_centers_per_trace=5)
        assert trace.trace_table_name is None
        assert trace.data.colnames == ['id', 'centers']
        assert len(trace.data['id']) == 0
        assert trace.data['centers'].shape == (0, 5)
        assert trace.data['centers'].description is not None
        assert trace.data['id'].description is not None

    def test_trace_instantiates_from_data(self):
        data = Table({'id': [1], 'centers': [[1, 2, 3]]})
        data['id'].description = 'test'
        data['centers'].description = 'test_2'
        trace = Trace(data=data)
        assert trace.trace_table_name is None
        for name in ['id', 'centers']:
            assert name in trace.data.colnames
        assert len(trace.data.colnames) == 2
        assert len(trace.data['id']) == 1
        assert trace.data['centers'].shape == (1, 3)
        assert trace.data['centers'].description == 'test_2'
        assert trace.data['id'].description == 'test'

    def test_trace_raises_exception(self):
        with pytest.raises(Exception):
            Trace(data=None, num_centers_per_trace=0)

    def test_getting_trace_centers(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert np.allclose(trace.get_centers(0), [0, 1])

    def test_getting_trace_id(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert trace.get_id(-1) == 1
        assert trace.get_id(0) == 0

    def test_getting_num_found_traces(self):
        trace = Trace(data={'id': [0, 1], 'centers': [[0, 1], [1, 2]]})
        assert trace.num_traces_found() == 2

    def test_add_centers(self):
        centers = np.arange(3)
        trace = Trace(data=None, num_centers_per_trace=len(centers))
        trace.add_centers(trace_centers=centers, trace_id=1)
        assert np.allclose(trace.data['centers'], [centers])
        assert np.allclose(trace.data['id'], [1])
        trace = Trace(data={'id': [1], 'centers': [centers]})
        trace.add_centers(trace_centers=centers, trace_id=2)
        assert np.allclose(trace.data['centers'], [centers, centers])
        assert np.allclose(trace.data['id'], [1, 2])

    def test_load_and_write(self):
        name = 'trace'
        trace = Trace(data={'id': [1], 'centers': [np.arange(3)]}, trace_table_name=name)
        with tempfile.TemporaryDirectory() as tmp_directory_name:
            filename = os.path.join(tmp_directory_name, 'test_trace_table.fits')
            trace.write(filename=filename)
            hdu_list = fits.open(filename)
            loaded_trace = Trace.load(hdu_list=hdu_list, trace_table_name=name)
        assert np.allclose(loaded_trace.get_centers(0), trace.get_centers(0))
        assert np.allclose(loaded_trace.get_id(0), trace.get_id(0))

    def test_sorting_trace_centers(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers-10, centers+2]}
        trace = Trace(data=data)
        trace.sort()
        assert np.allclose(trace.data['id'], np.arange(4))
        assert np.allclose(trace.data['centers'],
                           np.array([centers-10, centers, centers+2, centers+5]))

    def test_delete_centers(self):
        centers = np.array([1, 2, 3, 4])
        data = {'id': [1, 2, 3, 4],
                'centers': [centers, centers+5, centers+10, centers+11]}
        trace = Trace(data=data)
        trace.del_centers([])
        assert np.allclose(trace.data['id'], data['id'])
        assert np.allclose(trace.data['centers'], data['centers'])
        trace.del_centers(-1)
        assert np.allclose(trace.data['id'], [1, 2, 3])
        assert np.allclose(trace.data['centers'], [centers, centers+5, centers+10])
        trace.del_centers([-1, -2])
        assert np.allclose(trace.data['id'], [1])
        assert np.allclose(trace.data['centers'], [centers])


class TestAllTraceFitter:
    @mock.patch('banzai_nres.utils.trace_utils.AllTraceFitter._bad_shift')
    @mock.patch('banzai_nres.utils.trace_utils.AllTraceFitter._repeated_fit')
    @mock.patch('banzai_nres.utils.trace_utils.AllTraceFitter._beyond_edge')
    def test_bad_fit(self, beyond_edge, repeated_fit, bad_shift):
        all_trace_fitter = AllTraceFitter()
        beyond_edge.return_value, repeated_fit.return_value, bad_shift.return_value = True, True, True
        assert all_trace_fitter._last_fit_is_bad(trace=None, image_data=None, direction=None)
        beyond_edge.return_value, repeated_fit.return_value, bad_shift.return_value = False, False, False
        assert not all_trace_fitter._last_fit_is_bad(trace=None, image_data=None, direction=None)
        beyond_edge.return_value, repeated_fit.return_value, bad_shift.return_value = True, False, True
        assert all_trace_fitter._last_fit_is_bad(trace=None, image_data=None, direction=None)

    def test_detecting_repeated_fit(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1, 2], 'centers': [centers, centers+0.1]}
        trace = Trace(data=data)
        assert AllTraceFitter._repeated_fit(trace=trace)

    def test_detecting_bad_shifts(self):
        centers = np.array([1, 2, 3])
        data = {'id': [1], 'centers': np.array([centers])}
        trace = Trace(data=data)
        assert not AllTraceFitter._bad_shift(trace=trace, direction='up')

        trace.data = {'id': [1, 2], 'centers': np.array([centers, centers-0.1])}
        assert AllTraceFitter._bad_shift(trace=trace, direction='up')

        trace.data = {'id': [1, 2], 'centers': np.array([centers, centers+0.1])}
        assert AllTraceFitter._bad_shift(trace=trace, direction='down')

    def test_detecting_beyond_edge(self):
        fake_image_data = np.zeros((3, 3))
        data = {'id': [1], 'centers': [[3.1, 3.1, 3.1]]}
        assert AllTraceFitter._beyond_edge(Trace(data=data), fake_image_data)
        data = {'id': [1], 'centers': [[-0.1, -0.1, -0.1]]}
        assert AllTraceFitter._beyond_edge(Trace(data=data), fake_image_data)
        data = {'id': [1], 'centers': [[2, -0.1, 2]]}
        assert AllTraceFitter._beyond_edge(Trace(data=data), fake_image_data)

    def test_step_through_detector(self):
        # TODO
        assert True

    def test_fit_traces(self):
        # TODO
        assert True

    def test_window_for_next_trace_search(self):
        window, step = 100, 6
        ref_x = 1
        fitter = AllTraceFitter(march_parameters={'window': window, 'step_size': step})
        current_trace_centers = np.zeros(10)
        reference_y = current_trace_centers[ref_x]
        yminmax = fitter._window_for_next_trace_search(current_trace_centers,
                                                       reference_x=ref_x, direction='up')
        assert np.allclose(yminmax, (reference_y + step, reference_y + window + step))
        yminmax = fitter._window_for_next_trace_search(current_trace_centers,
                                                       reference_x=ref_x, direction='down')
        assert np.allclose(yminmax, (reference_y - window - step, reference_y - step))


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
        assert np.allclose(fitter.initial_guess_next_fit, np.array([1, 0, 90]))

    def test_generating_initial_guess_raises_error(self):
        with pytest.raises(Exception):
            SingleTraceFitter(image_data=np.zeros((2, 2)),
                              poly_fit_order=2,
                              start_point=None,
                              second_order_coefficient_guess=90)

    def test_changing_initial_guesses(self):
        coefficients = [np.array([0, 0])]
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'coefficients': coefficients})
        fitter.use_fit_as_initial_guess(-1)
        assert np.allclose(fitter.initial_guess_next_fit, coefficients[-1])
        fitter.initial_guess_next_fit += 1
        assert not np.allclose(fitter.initial_guess_next_fit, coefficients[-1])

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
        fitter.initial_guess_next_fit = single_trace_coefficients
        assert np.allclose(fitter.centers_current_guess(), a_line)

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

    @mock.patch('banzai_nres.utils.trace_utils.SingleTraceFitter._centers_from_coefficients')
    @mock.patch('banzai_nres.utils.trace_utils.optimize.minimize')
    def test_fit_trace(self, minimize, centers):
        coefficients = np.arange(4)
        centers.return_value = None
        minimize.return_value = OptimizeResult(x=coefficients)
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False,
                                              'coefficients': [],
                                              'initial_guess_next_fit': coefficients})

        assert fitter.fit_trace() is None
        assert fitter.coefficients == [coefficients]


class TestMatchFilter:
    """
    Tests for the functions in single trace fitter which are used to find
    the approximate locations of the next trace during a march up the detector.
    """
    @mock.patch('banzai_nres.utils.trace_utils.SingleTraceFitter._flux_across_trace', side_effect=np.max)
    def test_convolve_trace_with_image_data(self, flux_across_trace):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        current_trace_centers = np.arange(10)
        offset = current_trace_centers[2]
        flux_per_trace, coords = fitter._convolve_trace_with_image_data(current_trace_centers,
                                                                        y_min=0, y_max=4, reference_x=2)
        expected_flux_per_trace = (np.max(current_trace_centers) - offset) + np.arange(0, 4)
        assert np.allclose(flux_per_trace, expected_flux_per_trace)

    @mock.patch('banzai_nres.utils.trace_utils.SingleTraceFitter._flux_across_trace')
    @mock.patch('banzai_nres.utils.trace_utils.SingleTraceFitter._convolve_trace_with_image_data')
    def test_match_filter(self, flux_vs_pos, reference_flux):
        fitter = SingleTraceFitter(extraargs={'initialize_fit_objects': False})
        fitter.match_filter_parameters = {'min_peak_spacing': 5, 'neighboring_peak_flux_ratio': 20}
        positive_trace_signal, centroids, x_coords = array_with_two_peaks()
        no_trace_signal = np.random.normal(loc=0, scale=np.max(positive_trace_signal)/100,
                                           size=len(positive_trace_signal))
        reference_flux.return_value = np.max(positive_trace_signal)
        for trace_signal, outcome in zip([positive_trace_signal, no_trace_signal], [centroids[0], None]):
            fitter.initial_guess_next_fit = [0]
            flux_vs_pos.return_value = trace_signal
            next_trace_xy_coordinate = fitter.match_filter(None, y_min=None, y_max=None, reference_x=30)
            if next_trace_xy_coordinate is not None:
                assert np.isclose(next_trace_xy_coordinate[1], outcome, atol=3, rtol=0)
                assert np.isclose(next_trace_xy_coordinate[0], 30)


class TestTraceMaker:
    def test_properties(self):
        assert TraceMaker(FakeContext(settings=banzai_nres.settings.NRESSettings())).calibration_type is 'TRACE'

    def test_trace_fit_does_not_crash_on_blank_frame(self):
        readnoise = 11.0
        order_of_poly_fit = 4
        image = FakeTraceImage(nx=100, ny=100)
        image.readnoise = readnoise
        noisify_image(image)
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        fake_context.db_address = ''
        trace_fitter = TraceMaker(fake_context)
        trace_fitter.order_of_poly_fit = order_of_poly_fit
        trace_fitter.do_stage([image])
        assert True

    @mock.patch('banzai_nres.utils.trace_utils.AllTraceFitter.fit_traces')
    def test_trace_maker(self, fit_traces):
        trace_table_name = 'test'
        data = {'id': [1], 'centers': [np.arange(3)]}
        fit_traces.return_value = Trace(data=data)
        expected_trace = Trace(data=data)
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        trace_maker = TraceMaker(fake_context)
        trace_maker.trace_table_name = trace_table_name
        traces = trace_maker.do_stage(images=[FakeImage()])
        loaded_trace = traces[0][0]
        assert np.allclose(loaded_trace.get_centers(0), expected_trace.get_centers(0))
        assert np.allclose(loaded_trace.get_id(0), expected_trace.get_id(0))

    def test_accuracy_of_trace_fitting(self):
        """
        test type: Mock Integration Test with metrics for how well trace fitting is doing.
        info: This tests trace making via a blind fit.
        WARNING: Because trace fitting is defined with polynomials which are normalized from -1 to 1, if one squeezes
        the x axis of the image, then the traces bend more drastically. Thus it is recommended you do not change the
        size of the FakeTraceImage.
        """
        readnoise = 11.0
        poly_fit_order = 4

        image = FakeTraceImage()
        image.fiber0_lit, image.fiber1_lit, image.fiber2_lit = False, True, True
        image.readnoise = readnoise

        image, trace_centers, second_order_coefficient_guess = fill_image_with_traces(image,
                                                                                      poly_fit_order=poly_fit_order)
        noisify_image(image)
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        images = [image]

        trace_maker = TraceMaker(fake_context)
        trace_maker.order_of_poly_fit = poly_fit_order
        trace_maker.second_order_coefficient_guess = second_order_coefficient_guess
        traces = trace_maker.do_stage(images)
        assert traces[0][0].data['centers'].shape[0] == trace_centers.shape[0]
        difference = traces[0][0].data['centers'] - trace_centers
        logger.debug('median absolute deviation in unit-test trace fitting is {0} pixels'
                     .format(np.median(np.abs(difference - np.median(difference)))))
        logger.debug('standard deviation in unit-test trace fitting is {0} pixels'
                     .format(np.std(difference)))
        logger.debug('worst error (max of true minus found) in unit-test trace fitting is {0} pixels'
                     .format(np.max(np.abs(difference))))
        logger.debug('median error (median of true minus found) in unit-test trace fitting is {0} pixels'
                     .format(np.abs(np.median(difference))))
        assert np.median(np.abs(difference - np.median(difference))) < 1/100
        assert np.abs(np.median(difference)) < 1/100


class TestLoadTrace:
    def test_properties(self):
        assert LoadTrace(FakeContext(settings=banzai_nres.settings.NRESSettings())).calibration_type is 'TRACE'

    @mock.patch('banzai.dbs.get_master_calibration_image', return_value=None)
    def test_load_trace_removes_images_without_calibration(self, mock_get_cal):
        images = [FakeImage(), FakeImage(), FakeImage()]
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        setattr(fake_context, 'db_address', None)
        trace_loader = LoadTrace(fake_context)
        images = trace_loader.do_stage(images=images)
        assert len(images) == 0

    @mock.patch('banzai.dbs.get_master_calibration_image', return_value='path')
    @mock.patch('os.path.exists', return_value=True)
    @mock.patch('astropy.io.fits.open', return_value=None)
    @mock.patch('banzai_nres.utils.trace_utils.Trace.load')
    def test_load_trace(self, mock_load, mock_open, mock_os, mock_get_cal):
        data = {'id': [1], 'centers': [np.arange(3)]}
        expected_trace = Trace(data=data)
        mock_load.return_value = Trace(data=data)
        fake_context = FakeContext(settings=banzai_nres.settings.NRESSettings())
        setattr(fake_context, 'db_address', None)
        trace_loader = LoadTrace(fake_context)
        images = [FakeImage(), FakeImage(), FakeImage()]
        images = trace_loader.do_stage(images=images)
        for image in images:
            assert np.allclose(image.trace.get_centers(0), expected_trace.get_centers(0))
            assert np.allclose(image.trace.get_id(0), expected_trace.get_id(0))
