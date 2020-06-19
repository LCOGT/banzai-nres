from banzai_nres.rv import cross_correlate, c, barycentric_correction
import numpy as np



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


def test_rv():
    # Make fake data from orders 75-100
    # Run the RV code on the solstice
    pass

def test_bc_correction():
    #Celestial coordinates of the North Ecliptic Pole
    ra, dec = 18.0, 66.+33./60.+38.55/3600.
    site = 'elp'
    time, exptime = '2019-07-10T07:27:13.647', 0.0
    bc_rv, bjd_tdb = barycentric_correction(time, exptime, ra, dec, site)
    #assert np.close(np.abs(bc_rv),0.0,100) #Check RV correction is within 100 m/s of 0
    pass