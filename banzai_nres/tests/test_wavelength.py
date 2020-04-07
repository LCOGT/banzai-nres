import numpy as np
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace
from banzai_nres.wavelength import IdentifyFeatures, WavelengthCalibrate, get_ref_ids_and_fibers
from scipy.ndimage import gaussian_filter
from banzai_nres.frames import EchelleSpectralCCDData, NRESObservationFrame
from banzai import context
import sep
import pytest
import mock


class TestIdentifyFeatures:
    sigma = 1.0
    data = np.zeros((100, 100))
    err = np.ones_like(data)
    xcoords = [10, 30, 50, 70]  # 4 features
    ycoords = np.array(xcoords) + 5  # 4 features
    data[ycoords, xcoords] = 1
    data = gaussian_filter(data, sigma=sigma)
    data /= np.max(data) # make features have peak fluxes of 1

    def test_finds_features(self):
        features = identify_features(self.data, self.err, nsigma=0.5, fwhm=self.sigma)
        assert np.allclose(features['peak'], 1, atol=0.001)
        assert np.allclose(features['pixel'], self.xcoords, atol=0.001)
        assert np.allclose(features['ycentroid'], self.ycoords, atol=0.001)
        assert len(features) == 4

    def test_ignores_features(self):
        features = identify_features(self.data, self.err, nsigma=1.5, fwhm=self.sigma)
        assert len(features) == 0

    def test_extract(self):
        # small test to make sure sep.sum_circle is behaving.
        fluxes, _, _ = sep.sum_circle(self.data, self.xcoords, self.ycoords, 10 * self.sigma, gain=1.0, err=self.err)
        # assert the summed flux is the expected flux for a 2d (unnormalized) gaussian.
        assert np.allclose(fluxes, 2 * np.pi * self.sigma**2, rtol=1E-4)

    @pytest.mark.integration
    def test_do_stage(self):
        blaze_factor = 0.5
        input_context = context.Context({})
        image = NRESObservationFrame([EchelleSpectralCCDData(data=self.data, uncertainty=self.err,
                                                             meta={'OBJECTS': 'tung&tung&none'},
                                       traces=np.ones_like(self.data, dtype=int),
                                       blaze={'blaze': blaze_factor * np.ones((1, self.data.shape[1]), dtype=int)})], 'foo.fits')
        stage = IdentifyFeatures(input_context)
        stage.fwhm, stage.nsigma = self.sigma, 0.5
        image = stage.do_stage(image)
        assert np.allclose(image.features['corrected_flux'], image.features['flux'] / blaze_factor, rtol=1E-4)
        assert np.allclose(image.features['pixel'], self.xcoords, atol=0.001)
        assert np.allclose(image.features['ycentroid'], self.ycoords, atol=0.001)
        assert np.allclose(image.features['id'], 1)

    @pytest.mark.integration
    def test_do_stage_on_empty_features(self):
        input_context = context.Context({})
        image = NRESObservationFrame([EchelleSpectralCCDData(data=self.data, uncertainty=self.err,
                                      meta={'OBJECTS': 'tung&tung&none'},
                                      traces=np.ones_like(self.data, dtype=int),
                                      blaze={'blaze': np.ones_like(self.data, dtype=int)})], 'foo.fits')
        stage = IdentifyFeatures(input_context)
        stage.fwhm, stage.nsigma = self.sigma, 1.5
        image = stage.do_stage(image)
        assert len(image.features) == 0


class TestWavelengthCalibrate:
    def generate_image(self):
        traces = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 0], [0, 0, 0], [2, 2, 2], [3, 3, 3]])
        line_list = np.array([10, 11, 12])
        features = {'pixel': np.array([1, 2, 1, 2]), 'id': np.array([1, 1, 2, 3])}
        image = NRESObservationFrame([EchelleSpectralCCDData(data=np.zeros((2, 2)), uncertainty=np.zeros((2, 2)),
                                      meta={'OBJECTS': 'tung&tung&none'}, features=features,
                                      traces=traces, line_list=line_list)], 'foo.fits')
        return image

    @mock.patch('banzai_nres.wavelength.WavelengthCalibrate.calibrate_fiber',
                side_effect=lambda x, y, z, image, t: image)
    def test_do_stage(self, mock_calibrate_fiber):
        input_context = context.Context({})
        image = WavelengthCalibrate(input_context).do_stage(self.generate_image())
        assert np.allclose(image.features['order'], [0, 0, 0, 1])
        assert np.allclose(image.features['fiber'], [0, 0, 1, 0])

    def test_calibrate_fiber(self):
        # TODO
        assert True


@mock.patch('banzai_nres.wavelength.WavelengthSolution.solve', side_effect=lambda x, y, z: (x, y, z))
def test_recalibrate(mock_solve):
    measured_lines = {'wavelength': [10, 11, 50], 'pixel': [1, 2, 3], 'order': [0, 0, 0]}
    line_list = np.array([10.05, 11.05, 60, 11.5])
    pixel, order, principle_order_number = np.arange(10), np.arange(10), 30
    wavelength_solution = WavelengthCalibrate.recalibrate(measured_lines, line_list, pixel, order, principle_order_number)
    # the mock patch populated wavelength_solution.model_coefficients with the arguments fed to WavelengthSolution.solve
    # could use mock_solve.assert_called_with here, but it is not as straightforward because measured_lines changes.
    measured_lines, wavelengths_to_fit, weights = wavelength_solution.model_coefficients
    assert np.allclose(wavelengths_to_fit, line_list[:3])
    assert np.allclose(weights, [1, 1, 0])


def test_group_features_by_trace():
    traces = np.array([[0, 0, 1, 1], [2, 2, 0, 0]])
    features = {'xcentroid': [0, 2, 0, 2], 'ycentroid': [0, 0, 1, 1]}
    features = group_features_by_trace(features, traces)
    assert np.allclose(features['id'], [0, 1, 2, 0])


def test_get_principle_order_number():
    m0_values, ref_id = np.arange(10, 100), np.arange(20)
    true_m0 = 30
    center_wavelengths = np.random.normal(5000, 10, len(ref_id))/(true_m0 + ref_id)
    assert true_m0 == WavelengthCalibrate.get_principle_order_number(m0_values, center_wavelengths, ref_id)


def test_get_principle_order_number_warns_on_degenerate_m0():
    m0_values, ref_id = 30*np.ones(40), np.arange(20)
    true_m0 = 30
    center_wavelengths = np.random.normal(5000, 50, len(ref_id))/(true_m0 + ref_id)
    WavelengthCalibrate.get_principle_order_number(m0_values, center_wavelengths, ref_id)
    # right now this just tests whether the warning if statement does not cause a crash.
    assert True


def test_get_center_wavelengths():
    trace_ids = np.array([1, 2])
    traces = np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [0, 0, 0], [0, 0, 0], [2, 2, 2]])
    center_wavelengths = WavelengthCalibrate.get_center_wavelengths(traces * 3, traces, trace_ids)
    assert np.allclose(center_wavelengths, [3, 6])


def test_get_ref_ids_and_fibers():
    ref_id, fibers = get_ref_ids_and_fibers(3)
    assert np.allclose(ref_id, [0, 0, 1])
    assert np.allclose(fibers, [0, 1, 0])
    ref_id, fibers = get_ref_ids_and_fibers(6)
    assert np.allclose(ref_id, [0, 0, 1, 1, 2, 2])
    assert np.allclose(fibers, [0, 1, 0, 1, 0, 1])
