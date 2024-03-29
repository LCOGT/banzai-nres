import numpy as np
import pytest
from banzai_nres.utils import continuum_utils as cu
from banzai_nres.continuum import ContinuumNormalizer
from scipy.ndimage import gaussian_filter1d
from banzai_nres.fitting import fit_polynomial


SHOW_CONTINUUM_FITTING_ON_TEST_CASES = False


@pytest.mark.parametrize('random_seed', [2131209899, 21311222, 6913219, 75322, 50981234, 9942111])
class TestContinuumFittingUtilities:
    n_lines = 30
    line_sigma = 2
    size = 800
    broad_line_center = size / 2
    broad_line_width = 40
    continuum = np.linspace(100, 200, size)

    def make_spectrum(self, seed, broad_line_width=None):
        np.random.seed(seed)
        flux = np.random.normal(loc=self.continuum, scale=1, size=self.size)
        if broad_line_width is None:
            broad_line_width = self.broad_line_width + -4 * (np.random.random() - 0.5)

        sharp_lines = np.zeros_like(flux)
        line_locs = np.random.randint(0, len(sharp_lines), size=self.n_lines)
        sharp_lines[line_locs] = 60
        sharp_lines = -gaussian_filter1d(sharp_lines, sigma=self.line_sigma)
        broad_lines = -lorentz_profile(np.arange(len(flux)), self.broad_line_center, broad_line_width) * 3000
        return flux, sharp_lines, broad_lines, {'sharp': line_locs.astype(int), 'broad': [int(self.broad_line_center)]}

    def run_masking_procedure(self, flux):
        # testing the basic single iteration masking procedure, utilities that are used in ContinuumNormalizer
        mask = cu.mark_features(flux, sigma=3, detector_resolution=4)
        return mask

    def fit_passes(self, flux, mask, atol=1):
        # testing the basic single iteration fitting procedure, utilities that are used in ContinuumNormalizer
        flux_error = np.sqrt(np.ones_like(flux) ** 2 + np.sqrt(flux) ** 2)  # read plus poisson
        x = np.arange(len(flux))
        best_fit = fit_polynomial(flux, flux_error, x=x, order=3, sigma=3, mask=mask)
        # because the continuum is 100, an atol of 1 represents continuum normalization within 1%
        if SHOW_CONTINUUM_FITTING_ON_TEST_CASES:
            # debug
            import matplotlib.pyplot as plt
            plt.plot(flux, label='flux')
            plt.plot(np.arange(len(flux))[~mask], flux[~mask], label='unmasked flux')
            plt.plot(best_fit(np.arange(len(flux))), label='best fit continuum model')
            plt.title(f'max error {np.max(np.abs(best_fit(x) - self.continuum))}')
            plt.legend(loc='best')
            plt.show()
        return np.allclose(best_fit(x), self.continuum, atol=atol)

    @staticmethod
    def fraction_lines_masked(mask, locations):
        return np.count_nonzero(np.isclose(mask[locations], 1)) / len(locations)

    def test_null_reject(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        mask = self.run_masking_procedure(flux)
        # check that the masking was good enough for the continuum fit to recover the continuum.
        assert self.fit_passes(flux, mask, atol=1)

    def test_rejecting_sharp_lines(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        flux = flux + sharp_lines
        mask = self.run_masking_procedure(flux)
        # check that most the absorption lines have been marked by the mask
        assert self.fraction_lines_masked(mask, locs['sharp']) > .95
        # check that the masking was good enough for the continuum fit to recover the continuum.
        assert self.fit_passes(flux, mask, atol=5)

    def test_rejecting_narrow_broad_lines(self, random_seed):
        # In this test, the algorithm thinks the broad lines are much broader than they actually are in the spectrum.
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed, broad_line_width=15)
        flux = flux + broad_lines
        mask = self.run_masking_procedure(flux)
        # check that all the absorption lines have been marked by the mask
        assert self.fraction_lines_masked(mask, locs['broad']) == 1
        # check that the masking was good enough for the continuum fit to recover the continuum.
        assert self.fit_passes(flux, mask, atol=5)


@pytest.mark.parametrize('random_seed', [2131209899, 21311222, 6913219, 75322, 50981234, 9942111])
class TestContinuumFitting:
    # TODO move this to test_continuum.py
    # similar to TestBasicContinuumFitting but runs using the actual fit procedure in continuum.py
    n_lines = 30
    line_sigma = 2
    size = 800
    broad_line_center = size / 2
    broad_line_width = 40
    continuum = np.linspace(100, 200, size)

    def make_spectrum(self, seed, broad_line_width=None):
        np.random.seed(seed)
        flux = np.random.normal(loc=self.continuum, scale=1, size=self.size)
        if broad_line_width is None:
            broad_line_width = self.broad_line_width + -4 * (np.random.random() - 0.5)

        sharp_lines = np.zeros_like(flux)
        line_locs = np.random.randint(0, len(sharp_lines), size=self.n_lines)
        sharp_lines[line_locs] = 60
        sharp_lines = -gaussian_filter1d(sharp_lines, sigma=self.line_sigma)
        broad_lines = -lorentz_profile(np.arange(len(flux)), self.broad_line_center, broad_line_width) * 3000
        return flux, sharp_lines, broad_lines, {'sharp': line_locs.astype(int), 'broad': [int(self.broad_line_center)]}

    def fit_passes(self, flux, atol=1):
        # testing the basic single iteration fitting procedure, utilities that are used in ContinuumNormalizer
        flux_error = np.sqrt(np.ones_like(flux) ** 2 + np.sqrt(flux) ** 2)  # read plus poisson
        norm_flux, norm_uncertainty = ContinuumNormalizer.normalize(flux, flux_error, np.arange(len(flux)))
        if SHOW_CONTINUUM_FITTING_ON_TEST_CASES:
            # debug
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(flux, label='flux')
            plt.plot(self.continuum, ls='--', lw=4, alpha=1, label='true continuum')
            plt.plot(flux/norm_flux, label='best fit continuum model')
            plt.title(f'max error {np.max(np.abs(flux/norm_flux - self.continuum))}')
            plt.legend(loc='best')

            fig, ax = plt.subplots(1, 2)
            ax[0].errorbar(x=np.arange(len(norm_flux)), y=norm_flux, yerr=norm_uncertainty)
            ax[1].errorbar(x=np.arange(len(norm_flux)), y=flux/self.continuum, yerr=flux_error/self.continuum)
            plt.show()
        return np.allclose(flux/norm_flux, self.continuum, atol=atol)


def lorentz_profile(x, center, width, amplitude=1):
    return amplitude / np.pi * 1 / 2 * width / ((x - center) ** 2 + 1 / 4 * width ** 2)
