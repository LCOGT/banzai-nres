from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from banzai_nres.fitting import fit_polynomial
import numpy as np

from scipy.ndimage.morphology import binary_dilation
from scipy.stats import median_absolute_deviation

# Wavelength regions where there are strong Balmer or other absorption lines. This is from CERES (thanks!)
WAVELENGTHS_TO_MASK = np.array([[6755, 6769], [6530, 6600], [4840, 4880], [4320, 4360],
                                [4085, 4120], [3950, 3990], [3880, 3910], [3825, 3850]])


class ContinuumNormalizer(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        for fiber, order in zip(*image.spectrum.fibers_and_orders):
            spectrum = image.spectrum[fiber, order]
            blaze_corrected_flux = spectrum['flux'] / spectrum['blaze']
            blaze_corrected_error = blaze_corrected_flux * np.sqrt((spectrum['uncertainty'] / spectrum['flux']) ** 2.0 + (spectrum['blaze_error'] / spectrum['blaze']) ** 2.0)
            detector_resolution = 4  # pixels
            # identify features quickly
            first_derivative = np.hstack([[0], blaze_corrected_flux[1:] - blaze_corrected_flux[:-1]])
            mad = median_absolute_deviation(first_derivative)
            features = np.abs(first_derivative) > 1.5 * mad
            # Mask features
            mask = np.zeros_like(blaze_corrected_flux, dtype=int)
            mask[features] = 1
            # Mask the prohibited wavelength regions.
            for mask_region in WAVELENGTHS_TO_MASK:
                mask[np.logical_and(spectrum['wavelength'] >= min(mask_region), spectrum['wavelength'] <= max(mask_region))] = 1
            # Mask the immediate +- pixels around each absorption feature as well.
            # 2 binary dilations work too, but we will binary dilate up to the detector resolution.
            mask = binary_dilation(mask, iterations=int(detector_resolution))
            # do the fit.
            best_fit = fit_polynomial(blaze_corrected_flux, blaze_corrected_error, x=spectrum['wavelength'],
                                      order=3, sigma=5, mask=mask)
            image.spectrum[fiber, order, 'normflux'] = blaze_corrected_flux / best_fit(spectrum['wavelength'])
            # Technically, this should take into account the fit uncertainty using e.g. MCMC but the covariance it
            # adds is complicated enough to track that we do not attempt that here.
            image.spectrum[fiber, order, 'normuncertainty'] = blaze_corrected_error / best_fit(spectrum['wavelength'])
        return image
