from banzai.stages import Stage
from banzai.frames import ObservationFrame
import pkg_resources
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy import constants
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris
import astropy.units as u

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


def barycentric_correction(time, exptime, ra, dec, site):
    # TODO: figure out something better than hard-coding the NRES site coordinates
    if site == 'lsc': lat, lon, height = '-30:10:2.64', '-70:48:17.28', 2198.
    if site == 'elp': lat, lon, height = '30:40:12', '-104:01:12', 2070.
    if site == 'cpt': lat, lon, height = '-32:22:48', '20:48:36', 1460.
    if site == 'tlv': lat, lon, height = '30:35:45', '34:45:48', 875.
    # The following is largely taken from the astropy example page
    site_location = EarthLocation.from_geodetic(lat=lat, lon=lon, height=height*u.m)
    sky_coordinates = SkyCoord(ra=ra*u.hourangle, dec=dec*u.deg)
    # The time given in the NRES header is the exposure start time; correct this to the midpoint
    # Eventually, we will want to use the flux-weighted mid-point calculated from emeter data
    obs_time = Time(time, location=site_location) + exptime/2.*u.s
    barycorr_rv = sky_coordinates.radial_velocity_correction(obstime=obs_time)
    bjd_tdb = obs_time.tdb + obs_time.light_travel_time(sky_coordinates)
    return barycorr_rv.to(u.m/u.s).value, bjd_tdb.to_value('jd')


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
        # Calculate the peak v (converting to m/s) in the spectrograph frame
        rv_measured = (np.mean(rvs_per_order)) * 1000
        # Compute the barycentric RV and observing time corrections
        rv_correction, bjd_tdb = barycentric_correction(image.header['DATE-OBS'],image.header['EXPTIME'],
                                                image.header['RA'],image.header['DEC'],image.header['SITEID'])
        # Correct the RV per Wright & Eastman (2014) and save in the header
        image.header['RV'] = rv_measured + rv_correction + rv_measured * rv_correction / c, 'Radial Velocity [m/s]'
        # The following assumes that the uncertainty on the barycentric correction is negligible w.r.t. that
        # on the RV measured from the CCF, which should generally be true
        image.header['RVERR'] = (np.std(rvs_per_order) * 1000, 'Radial Velocity Uncertainty [m/s]')
        image.header['TCORR'] = bjd_tdb, 'Exposure Mid-Time (Barycentric Julian Date)'
        image.header['TCORVERN'] = 'astropy.time.light_travel_time', 'Time correction code version'
        image.header['TCORCOMP'] = 'ROMER, CLOCK', 'Time corrections done'
        image.header['TCOREPOS'] = 'ERFA', 'Source of Earth position'
        image.header['TCORSYST'] = 'BJD_TDB ', 'Ref. frame_timesystem of TCORR column'
        image.header['PLEPHEM'] = solar_system_ephemeris.get, 'Source of planetary ephemerides'
        image.ccf = Table(ccfs)
        return image
