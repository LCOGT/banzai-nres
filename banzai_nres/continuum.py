from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from banzai_nres.fitting import fit_polynomial
import numpy as np


class ContinuumNormalizer(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        for fiber, order in zip(*image.spectrum.fibers_and_orders):
            spectrum = image.spectrum[fiber, order]
            normalized_flux = spectrum['flux'] / spectrum['blaze']
            normalized_error = normalized_flux * np.sqrt((spectrum['uncertainty'] / spectrum['flux']) ** 2.0 + (spectrum['blaze_error'] / spectrum['blaze']) ** 2.0)
            # set the weights of
            eighty_fifth_percentile = np.percentile(spectrum['flux'], 85)
            cont_points = (spectrum['flux'] > eighty_fifth_percentile).astype(int)
            initial_weights = (normalized_error ** -2.0) * cont_points
            best_fit = fit_polynomial(normalized_flux, normalized_error, x=spectrum['wavelength'], initial_weights=initial_weights)
            image.spectrum[fiber, order, 'normflux'] = normalized_flux / best_fit(spectrum['wavelength'])
            # Technically, this should take into account the fit uncertainty using e.g. MCMC but the covariance it
            # adds is complicated enough to track that we do not attempt that here.
            image.spectrum[fiber, order, 'normuncertainty'] = normalized_error / best_fit(spectrum['wavelength'])
        return image
