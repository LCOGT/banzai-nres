import numpy as np
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace, \
    get_principle_order_number, get_center_wavelengths
from banzai_nres.wavelength import IdentifyFeatures, WavelengthCalibrate, get_ref_ids_and_fibers
from banzai_nres.wavelength import ArcLoader, LineListLoader, ArcStacker
from scipy.ndimage import gaussian_filter
from banzai_nres.frames import EchelleSpectralCCDData, NRESObservationFrame
from banzai_nres.qc.qc_wavelength import AssessWavelengthSolution
from banzai import context
from astropy.table import Table
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
        image.features.sort('pixel')
        assert np.allclose(image.features['corrected_flux'], image.features['flux'] / blaze_factor, rtol=1E-4)
        assert np.allclose(image.features['pixel'], self.xcoords, atol=0.001)
        assert np.allclose(image.features['ycentroid'], self.ycoords, atol=0.001)
        assert np.allclose(image.features['id'], 1)

    @pytest.mark.integration
    def test_do_stage_no_blaze(self):
        input_context = context.Context({})
        image = NRESObservationFrame([EchelleSpectralCCDData(data=self.data, uncertainty=self.err,
                                                             meta={'OBJECTS': 'tung&tung&none'},
                                       traces=np.ones_like(self.data, dtype=int))], 'foo.fits')
        stage = IdentifyFeatures(input_context)
        stage.fwhm, stage.nsigma = self.sigma, 0.5
        image = stage.do_stage(image)
        image.features.sort('pixel')
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
        pixel_positions = np.array([1, 2, 1, 2])
        features = Table({'pixel': pixel_positions, 'id': np.array([1, 1, 2, 3]), 'order': np.array([1, 1, 2, 3]),
                          'wavelength': pixel_positions, 'centroid_err': np.ones_like(pixel_positions)})
        image = NRESObservationFrame([EchelleSpectralCCDData(data=np.zeros((2, 2)), uncertainty=np.zeros((2, 2)),
                                      meta={'OBJECTS': 'tung&tung&none'}, features=features,
                                      traces=traces, line_list=line_list)], 'foo.fits')
        return image

    @mock.patch('banzai_nres.wavelength.WavelengthCalibrate.refine_wavelengths')
    def test_do_stage_with_existing_wavelengths(self, mock_refine_wavelengths):
        # test that feature wavelengths are populated from the old solutions.
        image = self.generate_image()
        image.wavelengths = np.random.random(size=image.traces.shape)
        image.features['xcentroid'], image.features['ycentroid'] = np.array([0, 1, 2, 0]), np.array([2, 0, 1, 0])
        expected_wavelengths = image.wavelengths[image.features['ycentroid'], image.features['xcentroid']]
        image = WavelengthCalibrate(context.Context({})).do_stage(image)
        assert np.allclose(image.features['wavelength'], expected_wavelengths)

    @mock.patch('banzai_nres.wavelength.WavelengthCalibrate.refine_wavelengths')
    @mock.patch('banzai_nres.wavelength.find_feature_wavelengths', return_value=np.arange(4))
    def test_do_stage(self, mock_find_wavelengths, mock_refine_wavelengths):
        image = self.generate_image()
        image.features['id'] = np.ones_like(image.features['pixel'])  # so that only one fiber is populated.
        image = WavelengthCalibrate(context.Context({})).do_stage(image)
        assert np.allclose(image.features['wavelength'], np.arange(4))

    def test_refine_wavelengths(self):
        # TODO
        assert True

    def test_fit_wavelength_model(self):
        # TODO
        assert True

    def test_init_feature_labels(self):
        features = {'id': np.array([1, 2, 3, 3, 4, 5, 6])}
        features = WavelengthCalibrate.init_feature_labels(6, features)
        assert np.allclose(features['fiber'], [0, 1, 0, 0, 1, 0, 1])
        assert np.allclose(features['order'], [0, 0, 1, 1, 1, 2, 2])


@mock.patch('banzai_nres.wavelength.WavelengthSolution.solve', side_effect=lambda x, y, z: (x, y, z))
def test_recalibrate(mock_solve):
    features = {'wavelength': [10, 11, 50], 'pixel': [1, 2, 3], 'order': [0, 0, 1]}
    line_list = np.array([10.05, 11.05, 60, 11.5])
    wavelength_solution = WavelengthCalibrate(context.Context({})).fit_wavelength_model(features, line_list, 30)
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


@mock.patch('banzai_nres.utils.wavelength_utils.get_center_wavelengths')
def test_get_principle_order_number(mock_wavelengths):
    m0_values, ref_ids = np.arange(10, 100), np.arange(20)
    true_m0 = 30
    mock_wavelengths.return_value = np.random.normal(5000, 10, len(ref_ids))/(true_m0 + ref_ids)
    assert true_m0 == get_principle_order_number(m0_values, {'order': ref_ids})


@mock.patch('banzai_nres.utils.wavelength_utils.get_center_wavelengths')
def test_get_principle_order_number_warns_on_degenerate_m0(mock_wavelengths):
    m0_values, ref_ids = 30*np.ones(40), np.arange(20)
    true_m0 = 30
    mock_wavelengths.return_value = np.random.normal(5000, 50, len(ref_ids))/(true_m0 + ref_ids)
    get_principle_order_number(m0_values, {'order': ref_ids})
    # right now this just tests whether the warning if statement does not cause a crash.
    assert True


def test_get_center_wavelengths():
    features = {'pixel': np.array([2, 4, 1, 5, 3]), 'order': np.array([1, 1, 2, 2, 3]),
                'wavelength': np.array([20, 40, 10, 50, 30])}
    center_wavelengths = get_center_wavelengths(features)
    assert np.allclose(center_wavelengths, [30, 30, 30])


def test_get_ref_ids_and_fibers():
    ref_id, fibers, trace_ids = get_ref_ids_and_fibers(3)
    assert np.allclose(ref_id, [0, 0, 1])
    assert np.allclose(fibers, [0, 1, 0])
    ref_id, fibers, trace_ids = get_ref_ids_and_fibers(6)
    assert np.allclose(ref_id, [0, 0, 1, 1, 2, 2])
    assert np.allclose(fibers, [0, 1, 0, 1, 0, 1])


def test_stage_caltypes():
    assert ArcStacker(context.Context({})).calibration_type == 'DOUBLE'
    assert LineListLoader(context.Context({})).calibration_type == 'LINELIST'


class TestLineListLoader:
    stage = LineListLoader(context.Context({}))
    @mock.patch('banzai_nres.wavelength.LineListLoader.on_missing_master_calibration', return_value=None)
    def test_do_stage_aborts(self, fake_miss):
        stage = LineListLoader(context.Context({}))
        stage.LINE_LIST_FILENAME = 'some/bad/path'
        assert stage.do_stage('image') is None

    @mock.patch('numpy.genfromtxt', return_value=np.array([[1, 1]]))
    def test_do_stage(self, fake_load):
        image = type('image', (), {})
        image = self.stage.do_stage(image)
        assert np.allclose(image.line_list, [1])

    def test_apply_master_calibration(self):
        test_image = type('image', (), {})
        assert np.allclose(self.stage.apply_master_calibration(test_image, [1, 2]).line_list, [1, 2])


class TestArcLoader:
    stage = ArcLoader(context.Context({}))

    def test_double_on_missing_master_calibration(self):
        test_image = type('image', (), {'obstype': 'DOUBLE'})
        assert self.stage.on_missing_master_calibration(test_image).obstype == 'DOUBLE'

    @mock.patch('banzai_nres.wavelength.CalibrationUser.on_missing_master_calibration')
    def test_on_missing_master_calibration(self, mock_parent_miss):
        test_image = type('image', (), {'obstype': 'ELSE'})
        assert self.stage.on_missing_master_calibration(test_image) is None

    def test_apply_master_calibration(self):
        master_cal = type('image', (), {'wavelengths': [1, 2], 'filename': 'foo.fits', 'fibers': [0, 1]})
        test_image = type('image', (), {'wavelengths': None, 'meta': {}})
        assert np.allclose(self.stage.apply_master_calibration(test_image, master_cal).wavelengths, [1, 2])
        assert self.stage.calibration_type == 'DOUBLE'
        assert test_image.meta['L1IDARC'][0] == 'foo.fits'


class TestQCChecks:
    test_image = TestWavelengthCalibrate().generate_image()
    input_context = context.Context({})

    def test_qc_do_stage(self):
        image = AssessWavelengthSolution(self.input_context).do_stage(self.test_image)
        assert image is not None

    def test_qc_checks(self):
        Delta_lambda = AssessWavelengthSolution(self.input_context).calculate_delta_lambda(self.test_image,self.test_image.features['wavelength'])
        sigma_Dlambda, good_sigma_Dlambda, raw_chi_squared, good_chi_squared = AssessWavelengthSolution(self.input_context).calculate_1d_metrics(self.test_image,Delta_lambda)
        assert sigma_Dlambda >= good_sigma_Dlambda
        assert raw_chi_squared >= good_chi_squared
        x_diff_Dlambda, order_diff_Dlambda = AssessWavelengthSolution(self.input_context).calculate_2d_metrics(self.test_image,Delta_lambda)
        assert np.any(np.isfinite(x_diff_Dlambda))
        assert np.any(np.isfinite(order_diff_Dlambda))
