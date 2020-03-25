import numpy as np
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace, aperture_extract
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

    def test_aperture_extract(self):
        fluxes = aperture_extract(self.coords, self.coords, self.data, aperture_width=self.sigma * 10)
        # assert the summed flux is the expected flux for a 2d (unnormalized) gaussian.
        assert np.allclose(fluxes, 2 * np.pi * self.sigma**2, rtol=1E-4)


def test_group_features_by_trace():
    traces = np.array([[0, 0, 1, 1], [2, 2, 0, 0]])
    features = {'xcentroid': [0, 2, 0, 2], 'ycentroid': [0, 0, 1, 1]}
    features = group_features_by_trace(features, traces)
    assert np.allclose(features['id'], [0, 1, 2, 0])
