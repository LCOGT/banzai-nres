import numpy as np
from banzai_nres.fitting import fit_polynomial


class TestFitPolynomial:
    x = np.arange(100)
    poly_coeffs = [2E-3, 2E-2, 2E-1, 2]
    continuum = np.polyval(poly_coeffs, x)
    sigma = np.median(continuum)/30
    continuum_error = sigma * np.ones_like(continuum, dtype=float)
    # set the seed for the noise generation
    np.random.seed(18232201)
    noise = np.random.normal(0, sigma, size=len(continuum))
    outliers = np.zeros_like(continuum)
    num_outliers = 5
    outliers[np.random.randint(0, high=len(outliers), size=num_outliers)] = np.median(continuum)*100

    def test_perfect(self):
        best_fit = fit_polynomial(self.continuum, self.continuum_error, x=self.x)
        assert np.allclose(best_fit(self.x), self.continuum)

    def test_perfect_with_outliers(self):
        best_fit = fit_polynomial(self.continuum + self.outliers, self.continuum_error, x=self.x)
        assert np.allclose(best_fit(self.x), self.continuum, atol=0.05)

    def test_with_noise(self):
        best_fit = fit_polynomial(self.continuum + self.noise, self.continuum_error, x=self.x)
        assert np.allclose(best_fit(self.x), self.continuum, atol=self.sigma/2)

    def test_with_noise_and_outliers(self):
        best_fit = fit_polynomial(self.continuum + self.noise + self.outliers, self.continuum_error, x=self.x)
        assert np.allclose(best_fit(self.x), self.continuum, atol=self.sigma/2)
