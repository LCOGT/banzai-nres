import numpy as np
from astropy.table import Table
from banzai_nres.utils.tellurics import generate_telluric_mask


def test_generate_telluric_mask_small_window():
    spectrum = Table({'wavelength': np.array([1, 2, 3, 4, 5, 6])})
    telluric_spectrum = np.array([[1, 2, 3, 4, 5, 6], [1, 1, 1, 0.5, 1, 1]])
    mask = generate_telluric_mask(spectrum, telluric_spectrum, 0.5001, window=0.5)
    assert np.allclose(mask, [0, 0, 0, 1, 0, 0])


def test_generate_telluric_mask_wide_window():
    spectrum = Table({'wavelength': np.array([1, 2, 3, 4, 5, 6])})
    telluric_spectrum = np.array([[1, 2, 3, 4, 5, 6], [1, 1, 1, 0.5, 1, 1]])
    mask = generate_telluric_mask(spectrum, telluric_spectrum, 0.5001, window=1.1)
    assert np.allclose(mask, [0, 0, 1, 1, 1, 0])