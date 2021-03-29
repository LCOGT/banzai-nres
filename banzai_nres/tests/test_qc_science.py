import numpy as np

from banzai import context
from banzai.nres.qc.qc_science import CalculateScienceFrameMetrics
from banzai.nres.qc.qc_science import get_snr
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.frames import EchelleSpectralCCDData, Spectrum1D

class TestCalculateScienceFrameMetrics:
    nput_context = context.Context({})
    test_wavelengths = np.linspace(5100.0,5200.0,4096)
    test_flux =  -0.001 * test_wavelengths**2 + 10.3 * test_wavelengths
    test_uncertainty = np.sqrt(test_flux)
    snr_order = 90
    spectrum = {'wavelength': test_wavelengths, 'flux': test_flux, 'uncertainty': test_uncertainty, 'fiber': 0, 'order': snr_order}
    image = NRESObservationFrame([EchelleSpectralCCDData(np.zeros((1, 1)), spectrum=Spectrum1D(spectrum))], science_fiber=0,
                                 'test.fits')
    

    def test_do_stage_does_not_crash(self):
        image = CalculateScienceFrameMetrics(self.input_context).do_stage(self.test_image)
        assert image is not None

    def test_snr_calculation(self):
        snr = get_snr(image,snr_order)
        assert np.isclose(snr, np.max(test_flux/test_uncertainty), rtol=0.1)