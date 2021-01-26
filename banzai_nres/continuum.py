from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from banzai_nres.fitting import fit_polynomial
import numpy as np

from banzai_nres.utils.continuum_utils import mark_features

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
                blaze_corrected_uncertainty = blaze_corrected_flux * np.sqrt(
                    (spectrum['uncertainty'] / spectrum['flux']) ** 2.0 + (
                                spectrum['blaze_error'] / spectrum['blaze']) ** 2.0)
                # fit away any remaining residuals
                norm_flux, norm_uncertainty, mask = self.normalize(blaze_corrected_flux, blaze_corrected_uncertainty, spectrum['wavelength'])

                image.spectrum[fiber, order, 'normflux'] = norm_flux
                image.spectrum[fiber, order, 'normuncertainty'] = norm_uncertainty
        return image

    @staticmethod
    def normalize(norm_flux_init, norm_uncertainty_init, wavelength):
        # first guess at the normalized flux is just the input (e.g. the blaze corrected flux)
        norm_flux = 1. * norm_flux_init
        norm_uncertainty = 1. * norm_uncertainty_init

        #mask = np.zeros_like(norm_flux)
        #for iteration in range(3): # could do an iterative procedure here.
        mask = mark_features(norm_flux, sigma=3, detector_resolution=4)
        # Mask the prohibited wavelength regions. Consider masking before mark_absorption_or_emission_features
        for mask_region in WAVELENGTHS_TO_MASK:
            mask[np.logical_and(wavelength >= min(mask_region), wavelength <= max(mask_region))] = 1
        best_fit = fit_polynomial(norm_flux, norm_uncertainty, x=wavelength,
                                  order=3, sigma=1, mask=mask)
        norm_flux /= best_fit(wavelength)
        # Technically, the normalized uncertainty should take into account the fit uncertainty using e.g.
        # MCMC but the covariance it adds is complicated enough to track that we do not attempt that here.
        norm_uncertainty /= best_fit(wavelength)
        return norm_flux, norm_uncertainty, mask


class MaskBlueHookRegion(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        mask = image.spectrum.mask
        mask[:, :1000] = 1  # mask the first pixels of every 1d spectrum to remove the blue hook.
        mask[:, -300:] = 1  # mask the last pixels of every 1d spectrum to remove any red edge effects
        image.spectrum.mask = mask
        return image