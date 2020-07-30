from banzai_nres.frames import NRESObservationFrame
from banzai.stages import Stage
from banzai_nres.fitting import fit_polynomial
import numpy as np

from scipy.ndimage.morphology import binary_dilation
from scipy.signal import find_peaks


class ContinuumNormalizer(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        for fiber, order in zip(*image.spectrum.fibers_and_orders):
            spectrum = image.spectrum[fiber, order]
            normalized_flux = spectrum['flux'] / spectrum['blaze']
            normalized_error = normalized_flux * np.sqrt((spectrum['uncertainty'] / spectrum['flux']) ** 2.0 + (spectrum['blaze_error'] / spectrum['blaze']) ** 2.0)
            # identify absorption features. negative of Normalized flux runs from -1 to 0 ideally. Deepest lines at -1.
            detector_resolution = 4  # pixels
            peaks, properties = find_peaks(-1*normalized_flux, height=-0.75, distance=2*detector_resolution,
                                           prominence=0.2, width=[detector_resolution, 200])
            # TODO use signal.peak_widths (or maybe find_peaks has these in properties) to
            #  mask pixels later on using the peak widths, not the detector resolution.
            # Mask absorption features
            mask = np.zeros_like(normalized_flux, dtype=int)
            mask[peaks] = 1
            # Mask the immediate +- two pixels around each absorption feature as well.
            mask = binary_dilation(mask, iterations=int(detector_resolution/2))
            # do the fit.
            best_fit = fit_polynomial(normalized_flux, normalized_error, x=spectrum['wavelength'],
                                      order=3, sigma=2, mask=mask)
            image.spectrum[fiber, order, 'normflux'] = normalized_flux / best_fit(spectrum['wavelength'])
            # Technically, this should take into account the fit uncertainty using e.g. MCMC but the covariance it
            # adds is complicated enough to track that we do not attempt that here.
            image.spectrum[fiber, order, 'normuncertainty'] = normalized_error / best_fit(spectrum['wavelength'])
        return image
