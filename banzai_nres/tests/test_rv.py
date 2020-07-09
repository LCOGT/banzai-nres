from banzai_nres.rv import cross_correlate, c, RVCalculator, barycentric_correction
import numpy as np
import mock
from astropy.io import fits
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.frames import EchelleSpectralCCDData, Spectrum1D
from types import SimpleNamespace
from astropy.table import Table


def gaussian(x, mu, sigma):
    return (2.0 * np.pi) ** -0.5 / sigma * np.exp(-0.5 * (x - mu) ** 2 / sigma ** 2)


def test_cross_correlate():
    test_wavelengths = np.arange(5000.0, 6000.0, 0.01)
    flux = np.ones(test_wavelengths.shape) * 1000

    read_noise = 10.0
    # add in fake lines
    for i in range(30):
        central_wavelength = np.random.uniform() * 1000.0 + 5000.0
        width = 0.1  # angstroms
        sigma = width / 2 / np.sqrt(2.0 * np.log(2.0))
        flux += gaussian(test_wavelengths, central_wavelength, sigma) * 10000.0 * np.random.uniform()

    true_v = 1.0
    noisy_flux = np.random.poisson(flux).astype(float)
    noisy_flux += np.random.normal(0.0, read_noise, size=len(flux))

    uncertainty = np.sqrt(flux + read_noise**2.0)

    # in km/s
    velocity_steps = np.arange(-5.0, 5.0, 0.1)
    ccf = cross_correlate(velocity_steps, test_wavelengths * (1.0 + true_v / c), noisy_flux,
                          uncertainty, test_wavelengths, flux)
    assert np.argmax(ccf) == np.argmin(np.abs(velocity_steps - true_v))


@mock.patch('banzai.dbs.get_site')
@mock.patch('banzai_nres.rv.fits.open')
def test_rv(mock_fits, mock_db):
    # Make fake data from orders 74-101
    test_wavelengths = np.arange(4500.0, 6500.0, 0.01)
    flux = np.ones(test_wavelengths.shape) * 1000

    read_noise = 10.0
    # add in fake lines
    for i in range(1000):
        central_wavelength = np.random.uniform() * 2000.0 + 4500.0
        width = 0.1  # angstroms
        sigma = width / 2 / np.sqrt(2.0 * np.log(2.0))
        flux += gaussian(test_wavelengths, central_wavelength, sigma) * 10000.0 * np.random.uniform()

    noisy_flux = np.random.poisson(flux).astype(float)
    noisy_flux += np.random.normal(0.0, read_noise, size=len(flux))
    uncertainty = np.sqrt(flux + read_noise**2.0)

    mock_fits.return_value = fits.HDUList([fits.PrimaryHDU(),
                                           fits.BinTableHDU(Table({'wavelength': test_wavelengths, 'flux': flux}))])
    true_v = 1.205  # km / s

    # Make the fake image
    # Set the site to the north pole, and the ra and dec to the ecliptic pole. The time we chose is just
    # to minimize the rv correction
    header = fits.Header({'DATE-OBS': '2020-09-12T00:00:00.000000',  'RA': 18.0, 'DEC': 66.55911605,
                          'OBJECTS': 'foo&none&none', 'EXPTIME': 1200.0})
    site_info = {'longitude': 0.0, 'latitude': 90.0, 'elevation': 0.0}

    # Each NRES order is ~70 Angstroms wide
    spectrum = []
    for i in range(26):
        order = slice(i * 7000, (i + 1) * 7000, 1)
        row = {'wavelength': test_wavelengths[order] * (1.0 + true_v / c),
               'flux': noisy_flux[order], 'uncertainty': uncertainty[order], 'fiber': 0, 'order': i + 75}
        spectrum.append(row)

    image = NRESObservationFrame([EchelleSpectralCCDData(np.zeros((1, 1)), meta=header, spectrum=Spectrum1D(spectrum))],
                                 'test.fits')
    image.instrument = SimpleNamespace()
    image.instrument.site = 'npt'
    mock_db.return_value = SimpleNamespace(**site_info)

    # Run the RV code
    stage = RVCalculator(SimpleNamespace(db_address='foo'))
    image = stage.do_stage(image)

    # Assert that the true_v + rv_correction is what is in the header within 5 m/s
    assert np.abs(true_v * 1000.0 + image.meta['BARYCORR'] - image.meta['RV']) < 5.0


def test_bc_correction():
    # Celestial coordinates of the North Ecliptic Pole (decimal hours, degrees)
    ra, dec = 18.0, 66.55911605
    # From the north pole
    site_info = {'longitude': 0.0, 'latitude': 90.0, 'elevation': 0.0}
    time, exptime = '2020-09-12T00:00:00.000000', 0.0
    bc_rv, bjd_tdb = barycentric_correction(time, exptime, ra, dec, SimpleNamespace(**site_info))
    assert np.abs(bc_rv) < 10.0  # Check RV correction is within 10 m/s of 0
