from banzai_nres.continuum import ContinuumNormalizer
import mock
import numpy as np
from banzai_nres.frames import NRESObservationFrame, EchelleSpectralCCDData, Spectrum1D
from types import SimpleNamespace
from astropy.table import Table


@mock.patch('banzai_nres.continuum.ContinuumNormalizer.normalize')
def test_do_stage(mock_normalize):
    # tests that the normalized spectrum are assigned correctly.
    expected = 1
    # make a single order spectrum
    flux = np.arange(1, 101) * expected
    x = np.arange(1, 101)
    spectrum = Table({'wavelength': [x, 2*x], 'flux': [flux, 2 * flux], 'blaze': [flux, 2 * flux],
                      'blaze_error': [flux, 2 * flux], 'uncertainty': [flux, 2*flux], 'fiber': [0, 0], 'order': [1, 2]})
    image = NRESObservationFrame([EchelleSpectralCCDData(np.zeros((1, 1)), meta={'OBJECTS': 'tung&tung&none'},
                                                         spectrum=Spectrum1D(spectrum))],
                                 'test.fits')
    # make it so that ContinuumNormalizer.normalize just returns ones.
    mock_normalize.return_value = (expected * np.ones_like(flux), expected * np.ones_like(flux))

    # Run the normalizer code
    stage = ContinuumNormalizer(SimpleNamespace(db_address='foo'))
    image = stage.do_stage(image)
    for order in [1, 2]:
        assert np.allclose(image.spectrum[0, order]['normflux'], expected)
        assert np.allclose(image.spectrum[0, order]['normuncertainty'], expected)


def test_do_stage_does_not_fit_non_science_fiber():
    # make a single order spectrum
    flux = np.arange(1, 101)*2
    x = np.arange(1, 101)
    spectrum = Table({'wavelength': [x], 'flux': [flux], 'blaze': [flux],
                      'blaze_error': [flux], 'uncertainty': [flux], 'fiber': [1], 'order': [1]})
    image = NRESObservationFrame([EchelleSpectralCCDData(np.zeros((1, 1)), meta={'OBJECTS': 'tung&tung&none'},
                                                         spectrum=Spectrum1D(spectrum))],
                                 'test.fits')

    # Run the normalizer code
    stage = ContinuumNormalizer(SimpleNamespace(db_address='foo'))
    image = stage.do_stage(image)
    assert 'normflux' not in image.spectrum._table.colnames
    assert 'normuncertainty' not in image.spectrum._table.colnames
