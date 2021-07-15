from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from scipy.signal import medfilt
from scipy import interpolate
import pkg_resources
import numpy as np
from banzai_nres.utils.continuum_utils import mark_features
from banzai_nres.utils.tellurics import generate_telluric_mask

# Wavelength regions where there are strong Balmer or other absorption lines. This is from CERES :
# https://ui.adsabs.harvard.edu/link_gateway/2017PASP..129c4002B/doi:10.1088/1538-3873/aa5455

WAVELENGTHS_TO_MASK = np.array([[6755, 6769], [6530, 6600], [4840, 4880], [4320, 4360],
                                [4085, 4120], [3950, 3990], [3880, 3910], [3825, 3850]])


class ContinuumNormalizer(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        for fiber, order in zip(*image.spectrum.fibers_and_orders):
            if np.isclose(fiber, image.science_fiber):
                # TODO, because masked parts of the spectrum are thrown out, any masked regions
                #  in the middle of the spectrum will cause steps in the spectrum according to the continuum fitter.
                #   the continuum fitter here should be able to know that. See Note B below.
                spectrum = image.spectrum[fiber, order]
                blaze_corrected_flux = spectrum['flux'] / spectrum['blaze']
                blaze_corrected_uncertainty = blaze_corrected_flux * np.sqrt(
                    (spectrum['uncertainty'] / spectrum['flux']) ** 2.0 + (
                                spectrum['blaze_error'] / spectrum['blaze']) ** 2.0)
                # normalize away any remaining residuals from the blaze division.
                norm_flux, norm_uncertainty = self.normalize(blaze_corrected_flux, blaze_corrected_uncertainty,
                                                             spectrum['wavelength'])

                image.spectrum[fiber, order, 'normflux'] = norm_flux
                image.spectrum[fiber, order, 'normuncertainty'] = norm_uncertainty
        return image

    def normalize(self, flux, uncertainty, wavelength):
        continuum_model = self.get_continuum_model(flux, wavelength)
        norm_flux = flux / continuum_model
        norm_uncertainty = uncertainty / continuum_model
        return norm_flux, norm_uncertainty

    @staticmethod
    def get_continuum_model(norm_flux_init, wavelength, window=201):
        # first guess at the normalized flux is just the input (e.g. the blaze corrected flux)
        norm_flux = 1. * norm_flux_init
        x = np.arange(len(norm_flux))

        mask = np.zeros_like(norm_flux, dtype=bool)
        # Mask the prohibited wavelength regions. Consider masking before mark_absorption_or_emission_features
        for mask_region in WAVELENGTHS_TO_MASK:
            mask[np.logical_and(wavelength >= min(mask_region), wavelength <= max(mask_region))] = 1
        mask[np.logical_not(mask)] = mark_features(norm_flux[np.logical_not(mask)], sigma=3, detector_resolution=4)
        continuum_model = np.copy(norm_flux)
        # interpolate the masked regions
        continuum_model[mask] = interpolate.interp1d(x[np.logical_not(mask)], norm_flux[np.logical_not(mask)],
                                                     kind='nearest', bounds_error=False, fill_value='extrapolate')(
            x[mask])
        # TODO NOTE B: x[..] here I think should be wavelength, so that masked portions are automatically skipped over
        #  and not interpreted as a step in the spectrum.
        # smooth the model
        continuum_model = medfilt(continuum_model, window)
        return continuum_model


class MaskBlueHookRegion(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        mask = image.spectrum.mask
        #mask[:, :500] += 8  # mask the first pixels of every 1d spectrum to remove the blue hook.
        #mask[:, -50:] += 8  # mask the last pixels of every 1d spectrum to remove any red edge effects
        image.spectrum.mask = mask
        return image


class MaskTellurics(Stage):
    # NOTE: This stage should come after continuum normalizer so that this does not create gaps
    # in the spectrum that we try to fit with a smooth model in continuum normalizer. Also regions affected by
    # tellurics still tell us about the continuum -- we don't want to liberally cut those regions out when
    # fitting the continuum. So I think we always want MaskTellurics to come after continuum fitting.
    TELLURIC_FILENAME = pkg_resources.resource_filename('banzai_nres', 'data/telluric_spectrum_50percent_humidity.dat')

    def do_stage(self, image) -> NRESObservationFrame:
        telluric_spectrum = np.genfromtxt(self.TELLURIC_FILENAME)
        mask = generate_telluric_mask(image.spectrum.table, telluric_spectrum)
        image.spectrum.mask[mask] |= 16
        return image
