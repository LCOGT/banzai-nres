import numpy as np
from astropy.io import fits
from banzai import context
from banzai_nres.qc.qc_science import CalculateScienceFrameMetrics
from banzai_nres.qc.qc_science import get_snr
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.frames import EchelleSpectralCCDData, Spectrum1D


class TestCalculateScienceFrameMetrics:
    input_context = context.Context({})
    test_wavelengths = np.linspace(5100.0, 5200.0, 4096)
    test_flux = -0.001 * test_wavelengths**2 + 10.3 * test_wavelengths
    test_uncertainty = np.sqrt(test_flux)
    snr_order = 90
    spectrum = []
    header = fits.Header({'OBJECTS': 'test&none&none'})
    for i in range(5):
        order = snr_order + i
        row = {'wavelength': test_wavelengths, 'flux': test_flux, 'uncertainty': test_uncertainty,
               'fiber': 0, 'order': snr_order}
        spectrum.append(row)
    test_image = NRESObservationFrame([EchelleSpectralCCDData(np.zeros((1, 1)), meta=header,
                                                              spectrum=Spectrum1D(spectrum))], 'test.fits')

    def test_do_stage_does_not_crash(self):
        image = CalculateScienceFrameMetrics(self.input_context).do_stage(self.test_image)
        assert image is not None

    def test_snr_calculation(self):
        snr = get_snr(self.test_image, self.snr_order)
        assert np.isclose(snr, np.max(self.test_flux/self.test_uncertainty), rtol=0.1)