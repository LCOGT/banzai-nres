from banzai.stages import Stage
from banzai.frames import ObservationFrame
import pkg_resources
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy import constants
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris
from astropy import units
from banzai.utils import stats
import logging


logger = logging.getLogger('banzai')

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
    # The following is largely taken from the astropy example page
    site_location = EarthLocation.from_geodetic(lat=site.latitude*units.deg, lon=site.longitude*units.deg,
                                                height=site.elevation*units.m)
    sky_coordinates = SkyCoord(ra=ra, dec=dec, unit=(units.hourangle, units.deg))
    # The time given in the NRES header is the exposure start time; correct this to the midpoint
    # Eventually, we will want to use the flux-weighted mid-point calculated from emeter data
    obs_time = Time(time, location=site_location) + exptime / 2.0 * units.s
    barycorr_rv = sky_coordinates.radial_velocity_correction(obstime=obs_time)
    bjd_tdb = obs_time.tdb + obs_time.light_travel_time(sky_coordinates)
    return barycorr_rv.to('m/s').value, bjd_tdb.to_value('jd')


def cross_correlate_over_traces(image, orders_to_use, velocities, template):
    ccfs = []

    for i in orders_to_use:
        logger.info(f'Cross correlating for order: {i}')
        order = image.spectrum[np.logical_and(image.spectrum['fiber'] == image.science_fiber, image.spectrum['order'] == i)][0]
        # calculate the variance ahead of time
        x_cor = cross_correlate(velocities, order['wavelength'], order['flux'], order['uncertainty'],
                                template['wavelength'], template['flux'])
        ccfs.append({'order': i, 'v': velocities, 'xcor': x_cor})
    return Table(ccfs)


class RVCalculator(Stage):
    TEMPLATE_FILENAME = pkg_resources.resource_filename('banzai_nres', 'data/g2v_template.fits')

    def do_stage(self, image) -> ObservationFrame:
        # Load in the template
        template_hdu = fits.open(self.TEMPLATE_FILENAME)

        template = {'wavelength': template_hdu[1].data['wavelength'], 'flux': template_hdu[1].data['flux']}
        # This is mostly arbitrary. Just pick orders near the center of the detector
        orders_to_use = np.arange(75, 101, 1)

        # for steps in 1 km/s from -2000 to +2000 km/s
        velocities = np.arange(-2000, 2001, 1)

        coarse_ccfs = cross_correlate_over_traces(image, orders_to_use, velocities, template)

        # take the peak
        velocity_peaks = np.array([velocities[np.argmax(ccf['xcor'])] for ccf in coarse_ccfs])
        best_v = stats.sigma_clipped_mean(velocity_peaks, 3.0)
        velocities = np.arange(best_v - 2, best_v + 2 + 1e-4, 1e-3)

        final_ccfs = cross_correlate_over_traces(image, orders_to_use, velocities, template)

        rvs_per_order = np.array([ccf['v'][np.argmax(ccf['xcor'])] for ccf in final_ccfs])
        # Calculate the peak v (converting to m/s) in the spectrograph frame
        rv_measured = stats.sigma_clipped_mean(rvs_per_order, 3.0) * 1000

        # TODO: Add ra and dec to observation frame object
        # Compute the barycentric RV and observing time corrections
        rv_correction, bjd_tdb = barycentric_correction(image.dateobs, image.exptime,
                                                        image.meta['RA'], image.meta['DEC'], image.instrument.site)
        # Correct the RV per Wright & Eastman (2014) and save in the header
        rv = rv_measured + rv_correction + rv_measured * rv_correction / c
        image.meta['RV'] = rv, 'Radial Velocity in Barycentric Frame [m/s]'
        # The following assumes that the uncertainty on the barycentric correction is negligible w.r.t. that
        # on the RV measured from the CCF, which should generally be true
        image.meta['RVERR'] = stats.robust_standard_deviation(rvs_per_order) * 1000, 'Radial Velocity Uncertainty [m/s]'
        image.meta['BARYCORR'] = rv_correction, 'Barycentric Correction Applied to the RV [m/s]'
        image.meta['TCORR'] = bjd_tdb, 'Exposure Mid-Time (Barycentric Julian Date)'
        image.meta['TCORVERN'] = 'astropy.time.light_travel_time', 'Time correction code version'
        #TODO: verify that the following are in fact the time corrections done by astropy;
        #this isn't completely obvious from the online documentation
        image.meta['TCORCOMP'] = 'ROMER, CLOCK', 'Time corrections done'
        image.meta['TCOREPOS'] = 'ERFA', 'Source of Earth position'
        image.meta['TCORSYST'] = 'BJD_TDB ', 'Ref. frame_timesystem of TCORR column'
        image.meta['PLEPHEM'] = solar_system_ephemeris.get(), 'Source of planetary ephemerides'
        image.ccf = final_ccfs
        return image
