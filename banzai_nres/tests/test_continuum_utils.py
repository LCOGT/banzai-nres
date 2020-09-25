import numpy as np
import pytest
from banzai_nres.utils import continuum_utils as cu
from scipy.ndimage import gaussian_filter1d\


def lorentz_profile(x, center, width):
    return 1 / np.pi * 1 / 2 * width / ((x - center) ** 2 + 1 / 4 * width ** 2)


@pytest.mark.parametrize('random_seed', [2131209899, 21311222, 6913219])
class TestFeatureRejection:
    n_lines = 30
    line_sigma = 2
    size = 800
    broad_line_center = size / 2
    broad_line_width = 40

    def make_spectrum(self, seed):
        np.random.seed(seed)
        flux = np.random.normal(loc=100, scale=1, size=self.size)
        broad_line_width = self.broad_line_width + -4 * (np.random.random() - 0.5)

        sharp_line_locs = np.zeros_like(flux)
        sharp_line_locs[np.random.randint(0, len(sharp_line_locs), size=self.n_lines)] = 60
        sharp_lines = -gaussian_filter1d(sharp_line_locs, sigma=self.line_sigma)
        broad_lines = -lorentz_profile(np.arange(len(flux)), self.broad_line_center, broad_line_width) * 3000
        return flux, sharp_lines, broad_lines, {'sharp': sharp_line_locs, 'broad': self.broad_line_center}

    def run_masking_procedure(self, flux):
        mask = cu.mark_absorption_or_emission_features(flux, 2 * self.line_sigma)
        #mask = np.logical_or(mask, cu.mark_pressure_broadened_features(flux, self.broad_line_width))
        return mask

    def test_null_reject(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        mask = self.run_masking_procedure(flux)
        import matplotlib.pyplot as plt
        plt.plot(np.arange(len(flux)), flux)
        plt.plot(np.arange(len(flux))[~mask], flux[~mask])
        plt.show()
        # test that less than 1% of the spectrum is identified.
        assert np.count_nonzero(mask) < self.size/100

    def test_rejecting_sharp_lines(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        flux = flux + sharp_lines
        mask = self.run_masking_procedure(flux)
        import matplotlib.pyplot as plt
        plt.plot(np.arange(len(flux)), flux)
        plt.plot(np.arange(len(flux))[~mask], flux[~mask])
        plt.show()
        assert True

    def test_rejecting_broad_lines(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        flux = flux + broad_lines
        mask = self.run_masking_procedure(flux)
        import matplotlib.pyplot as plt
        plt.plot(np.arange(len(flux)), flux)
        plt.plot(np.arange(len(flux))[~mask], flux[~mask])
        plt.show()
        assert True

    def test_rejecting_sharp_and_broad_lines(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        flux = flux + sharp_lines + broad_lines
        mask = self.run_masking_procedure(flux)
        import matplotlib.pyplot as plt
        plt.plot(np.arange(len(flux)), flux)
        plt.plot(np.arange(len(flux))[~mask], flux[~mask])
        plt.show()
        assert True



