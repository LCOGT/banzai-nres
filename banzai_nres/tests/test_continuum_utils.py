import numpy as np
import pytest
from banzai_nres.utils import continuum_utils as cu
from scipy.ndimage import gaussian_filter1d
from banzai_nres.fitting import fit_polynomial, fit_stiff_polynomial


@pytest.mark.parametrize('random_seed', [2131209899, 21311222, 6913219, 75322, 50981234, 9942111])
class TestContinuumFitting:
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
        broad_lines = -cu.lorentz_profile(np.arange(len(flux)), self.broad_line_center, broad_line_width) * 3000
        return flux, sharp_lines, broad_lines, {'sharp': line_locs.astype(int), 'broad': [int(self.broad_line_center)]}

    def run_masking_procedure(self, flux):
        flux_error = np.sqrt(np.ones_like(flux) ** 2 + np.sqrt(flux) ** 2)  # read plus poisson
        mask = cu.mark_features(flux, sigma=3, detector_resolution=4)
        #broad_line_mask = cu.mark_broad_features(flux, flux_error, mask, self.broad_line_width/2, level_to_mask=1E-2)
        #mask = np.logical_or(mask, broad_line_mask)
        return mask

    def fit_passes(self, flux, mask, atol=1):
        flux_error = np.sqrt(np.ones_like(flux) ** 2 + np.sqrt(flux) ** 2)  # read plus poisson
        x = np.arange(len(flux))
        best_fit = fit_polynomial(flux, flux_error, x=x, order=3, sigma=3, mask=mask)
        #best_fit2 = fit_stiff_polynomial(flux, flux_error, x=x, order=3, sigma=3, mask=mask, derivative_stiffness=[0, 1E2, 1E10])
        # because the continuum is 100, an atol of 1 represents continuum normalization within 1%
        if True:
            # debug
            import matplotlib.pyplot as plt
            plt.plot(flux, label='flux')
            plt.plot(np.arange(len(flux))[~mask], flux[~mask], label='unmasked flux')
            #plt.plot(self.continuum, ls='--', lw=4, alpha=1, label='true continuum')
            plt.plot(best_fit(np.arange(len(flux))), label='best fit continuum model')
            #plt.plot(best_fit2(np.arange(len(flux))), label='best fit stiff model')
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

    def test_rejecting_broad_lines(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        flux = flux + broad_lines
        mask = self.run_masking_procedure(flux)
        # check that all the absorption lines have been marked by the mask
        assert self.fraction_lines_masked(mask, locs['broad']) == 1
        # check that the masking was good enough for the continuum fit to recover the continuum.
        assert self.fit_passes(flux, mask, atol=1)

    def test_rejecting_sharp_and_broad_lines(self, random_seed):
        flux, sharp_lines, broad_lines, locs = self.make_spectrum(random_seed)
        flux = flux + sharp_lines + broad_lines
        mask = self.run_masking_procedure(flux)

        assert self.fraction_lines_masked(mask, locs['broad']) == 1
        assert self.fraction_lines_masked(mask, locs['broad']) > 0.95
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
        assert self.fit_passes(flux, mask, atol=1)



