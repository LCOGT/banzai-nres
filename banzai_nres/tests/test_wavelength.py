import numpy as np
from scipy.ndimage.morphology import binary_dilation
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace
from scipy.ndimage import gaussian_filter
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
        assert np.allclose(features['pixel'], [10, 30, 50, 70], atol=0.001)
        assert np.allclose(features['ycentroid'], [10, 30, 50, 70], atol=0.001)
        assert len(features) == 4

    def test_ignores_features(self):
        features = identify_features(self.data, self.err, nsigma=1.5, fwhm=self.sigma)
        assert len(features) == 0


def test_group_features_by_trace():
    assert True