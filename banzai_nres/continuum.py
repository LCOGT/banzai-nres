from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from banzai_nres.fitting import fit_polynomial


class ContinuumNormalizer(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        for fiber, order in zip(image.spectrum.orders):
            spectrum = image.spectrum[fiber, order]
            normalized_flux = spectrum['flux'] / spectrum['blaze']
            normalized_error = spectrum['uncertainty'] / spectrum['blaze']
            best_fit = fit_polynomial(normalized_flux, normalized_error, x=spectrum['wavelength'])
            image.spectrum[fiber, order, 'normflux'] = normalized_flux / best_fit(spectrum['wavelength'])
            image.spectrum[fiber, order, 'normuncertainty'] = normalized_error / best_fit(spectrum['wavelength'])
        return image
