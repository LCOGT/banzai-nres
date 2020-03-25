import numpy as np
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace
from banzai_nres.wavelength import IdentifyFeatures
from scipy.ndimage import gaussian_filter
from banzai_nres.frames import EchelleSpectralCCDData
from banzai import context
import sep
import pytest


class TestIdentifyFeatures:
    sigma = 1.0
    data = np.zeros((100, 100))
    err = np.ones_like(data)
    coords = [10, 30, 50, 70]  # 4 features
    data[coords, coords] = 1
    data = gaussian_filter(data, sigma=sigma)
    data /= np.max(data) # make features have peak fluxes of 1

    def test_finds_features(self):
        features = identify_features(self.data, self.err, nsigma=0.5, fwhm=self.sigma)
        assert np.allclose(features['peak'], 1, atol=0.001)
        assert np.allclose(features['pixel'], self.coords, atol=0.001)
        assert np.allclose(features['ycentroid'], self.coords, atol=0.001)
        assert len(features) == 4

    def test_ignores_features(self):
        features = identify_features(self.data, self.err, nsigma=1.5, fwhm=self.sigma)
        assert len(features) == 0

    def test_extract(self):
        # small test to make sure sep.sum_circle is behaving.
        fluxes, _, _ = sep.sum_circle(self.data, self.coords, self.coords, 10 * self.sigma, gain=1.0, err=self.err)
        # assert the summed flux is the expected flux for a 2d (unnormalized) gaussian.
        assert np.allclose(fluxes, 2 * np.pi * self.sigma**2, rtol=1E-4)

    @pytest.mark.integration
    def test_do_stage(self):
        blaze_factor = 0.5
        input_context = context.Context({})
        image = EchelleSpectralCCDData(data=self.data, uncertainty=self.err, meta={'OBJECTS': 'tung&tung&none'},
                                       traces=np.ones_like(self.data, dtype=int),
                                       blaze={'blaze': blaze_factor * np.ones_like(self.data, dtype=int)})
        stage = IdentifyFeatures(input_context)
        stage.fwhm, stage.nsigma = self.sigma, 0.5
        image = stage.do_stage(image)
        assert np.allclose(image.features['corrected_flux'], image.features['flux'] / blaze_factor, rtol=1E-4)
        assert np.allclose(image.features['pixel'], self.coords, atol=0.001)
        assert np.allclose(image.features['ycentroid'], self.coords, atol=0.001)
        assert np.allclose(image.features['id'], 1)

    @pytest.mark.integration
    def test_do_stage_on_empty_features(self):
        input_context = context.Context({})
        image = EchelleSpectralCCDData(data=self.data, uncertainty=self.err, meta={'OBJECTS': 'tung&tung&none'},
                                       traces=np.ones_like(self.data, dtype=int),
                                       blaze={'blaze': np.ones_like(self.data, dtype=int)})
        stage = IdentifyFeatures(input_context)
        stage.fwhm, stage.nsigma = self.sigma, 1.5
        image = stage.do_stage(image)
        assert len(image.features) == 0


def test_group_features_by_trace():
    traces = np.array([[0, 0, 1, 1], [2, 2, 0, 0]])
    features = {'xcentroid': [0, 2, 0, 2], 'ycentroid': [0, 0, 1, 1]}
    features = group_features_by_trace(features, traces)
    assert np.allclose(features['id'], [0, 1, 2, 0])
