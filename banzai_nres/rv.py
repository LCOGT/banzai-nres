from banzai.stages import Stage
from banzai.frames import ObservationFrame
import pkg_resources
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy import constants

# Speed of light in km/s
c = constants.c.to('km / s').value


def cross_correlate(velocities, wavelength, flux, flux_uncertainty, template_wavelength, template_flux):
    """
    :param velocities: in km / s
    :param wavelength:
    :param flux:
    :param flux_uncertainty:
    :param template_wavelength:
    :param template_flux:
    :return:

    This is consistent with what we did in the satellite pipeline and
    Zackay & Ofek, 2017, ApJ, 836, 187
    """
    x_cor = []
    # for steps in 1 km/s from -2000 to +2000 km/s
    # calculate the variance ahead of time
    variance = flux_uncertainty * flux_uncertainty

    for v in velocities:
        doppler = 1.0 / (1.0 + v / c)
        kernel = np.interp(doppler * wavelength, template_wavelength, template_flux)
        correlation = kernel * flux
        correlation /= variance

        normalization = kernel * kernel / variance

        # This gives us the signal to noise ratio in terms of N-sigma
        snr = correlation.sum() / np.sqrt(normalization.sum())
        x_cor.append(snr)

    return np.array(x_cor)


class RVCalculator(Stage):
    TEMPLATE_FILENAME = pkg_resources.resource_filename('banzai_nres', 'data/g2v_template.fits')

    def do_stage(self, image) -> ObservationFrame:
        # Load in the template
        template_hdu = fits.open(self.TEMPLATE_FILENAME)
        ccfs = []
        for i in image.fibers[image.fibers['fiber'] == image.science_fiber]['trace']:
            # for steps in 1 km/s from -2000 to +2000 km/s
            velocities = np.arange(-2000, 2001, 1)
            # calculate the variance ahead of time
            x_cor = cross_correlate(velocities, image.spectrum[i]['wavelength'], image.spectrum[i]['flux'],
                                    image.spectrum[i]['uncertainty'],
                                    template_hdu[1].data['wavelength'], template_hdu[1].data['flux'])

            # take the peak
            best_v = velocities[np.argmax(x_cor)]
            # for steps of 1m/s around that peak take +- 2 km/s and repeat the process
            velocities = np.arange(best_v - 2, best_v + 2 + 1e-4, 1e-3)
            x_cor = cross_correlate(velocities, image.spectrum[i]['wavelength'], image.spectrum[i]['flux'],
                                    image.spectrum[i]['uncertainty'],
                                    template_hdu[1].data['wavelength'], template_hdu[1].data['flux'])

            # Save the 1 m/s cross correlation function
            ccfs.append({'v': velocities, 'xcor': x_cor})

        rvs_per_order = [ccf['v'][np.argmax(ccf['xcor'])] for ccf in ccfs]
        # Save the peak v in the header (converting to m/s)
        image.header['RV'] = (np.mean(rvs_per_order) * 1000, 'Radial Velocity [m/s]')
        image.header['RVERR'] = (np.std(rvs_per_order) * 1000, 'Radial Velocity Uncertainty [m/s]')
        image.ccf = Table(ccfs)
        return image
