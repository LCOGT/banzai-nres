from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from banzai_nres.fitting import fit_polynomial
import numpy as np

from banzai_nres.utils.continuum_utils import mark_absorption_or_emission_features

# Wavelength regions where there are strong Balmer or other absorption lines. This is from CERES :
# https://ui.adsabs.harvard.edu/link_gateway/2017PASP..129c4002B/doi:10.1088/1538-3873/aa5455

WAVELENGTHS_TO_MASK = np.array([[6755, 6769], [6530, 6600], [4840, 4880], [4320, 4360],
                                [4085, 4120], [3950, 3990], [3880, 3910], [3825, 3850]])


class ContinuumNormalizer(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        for fiber, order in zip(*image.spectrum.fibers_and_orders):
            if np.isclose(fiber, image.science_fiber):
                spectrum = image.spectrum[fiber, order]
                blaze_corrected_flux = spectrum['flux'] / spectrum['blaze']
                blaze_corrected_error = blaze_corrected_flux * np.sqrt((spectrum['uncertainty'] / spectrum['flux']) ** 2.0 + (spectrum['blaze_error'] / spectrum['blaze']) ** 2.0)
                detector_resolution = 4  # pixels
                mask = mark_absorption_or_emission_features(blaze_corrected_flux, int(detector_resolution))
                # Mask the prohibited wavelength regions. Consider masking before mark_absorption_or_emission_features
                for mask_region in WAVELENGTHS_TO_MASK:
                    mask[np.logical_and(spectrum['wavelength'] >= min(mask_region), spectrum['wavelength'] <= max(mask_region))] = 1
                best_fit = fit_polynomial(blaze_corrected_flux, blaze_corrected_error, x=spectrum['wavelength'],
                                          order=3, sigma=5, mask=mask)
                image.spectrum[fiber, order, 'normflux'] = blaze_corrected_flux / best_fit(spectrum['wavelength'])
                # Technically, this should take into account the fit uncertainty using e.g. MCMC but the covariance it
                # adds is complicated enough to track that we do not attempt that here.
                image.spectrum[fiber, order, 'normuncertainty'] = blaze_corrected_error / best_fit(spectrum['wavelength'])
        return image
