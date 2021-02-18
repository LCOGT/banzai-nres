import numpy as np
from xwavecal.utils.wavelength_utils import find_nearest


def generate_telluric_mask(spectrum, telluric_spectrum, cutoff_transmission=.9, window=4/40):
    """
    Note: all wavelengths of spectrum and telluric spectrum should be in Angstroms.
    :param: spectrum: astropy.table.Table. spectrum['wavelength'] should give an nd.array of wavelengths.
    :param: telluric_spectrum: 2d numpy array of shape [2, N]. Row [0] should be the wavelengths, row [1] the
    :param: normalized transmission (i.e. 1 means perfect transmission through the atmosphere, 0 means total absorption)
    :param: cutoff_transmission: float between 0 and 1. Any wavelengths corresponding to transmission less than this value will
    :param: be added to the mask
    :param: window: float. Wavelength window in angstroms. For each masked region, we will expand the mask by this amount
    in wavelengths. This will catch the wings of telluric lines.
    """
    mask = np.zeros(spectrum['wavelength'].shape)
    wavelengths_ignore = np.sort(telluric_spectrum[0][telluric_spectrum[1] < cutoff_transmission])
    if len(wavelengths_ignore) == 1:
        # force length 2 so that find_nearest does not fail.
        wavelengths_ignore = np.array([wavelengths_ignore[0], wavelengths_ignore[0]])
    if len(wavelengths_ignore) >= 2:
        # only try to generate the mask if there are any wavelengths we need to ignore.
        mask = np.isclose(spectrum['wavelength'].data - find_nearest(np.array(spectrum['wavelength'].data), wavelengths_ignore), 0, atol=window)
    return mask
