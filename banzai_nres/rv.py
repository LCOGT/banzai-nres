from banzai.stages import Stage
from banzai.frames import ObservationFrame
from banzai import dbs
import numpy as np
from astropy.table import Table
from astropy import constants
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris
from astropy import units
from banzai.utils import stats
import logging
from banzai_nres.fitting import fit_polynomial
from banzai_nres import phoenix
from banzai_nres.utils.continuum_utils import mark_features

logger = logging.getLogger('banzai')

# Speed of light in km/s
c = constants.c.to('km / s').value


# wavelength regions selected by eye by Mirek on 09 16 2020 wherever there are large (H alpha like, i.e. broad wings)
# absorption features in the G2V pheonix template (banzai nres name phoenix-05700-p4.5-p0.0-p0.0.fits).
# This shares many regions with the region masks in banzai_nres.continuum.WAVELENGTHS_TO_MASK
PHOENIX_WAVELENGTHS_TO_MASK = np.array([[9215, 9250], [9000, 9035], [8725, 8760], [8800, 8815], [8855, 8875], [8525, 8560], [8650, 8675],
                                        [8480, 8510], [6530, 6600], [5880, 5910], [5260, 5280], [5320, 5340], [5160, 5190], [4880, 4840],
                                        [4380, 4390], [4310, 4360], [4090, 4115], [4040, 4055], [4060, 4070], [3950, 3980], [3925, 3945]])


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
    sky_coordinates = SkyCoord(ra=ra, dec=dec, unit=(units.deg, units.deg))
    # The time given in the NRES header is the exposure start time; correct this to the midpoint
    # Eventually, we will want to use the flux-weighted mid-point calculated from emeter data
    obs_time = Time(time, location=site_location) + exptime / 2.0 * units.s
    barycorr_rv = sky_coordinates.radial_velocity_correction(obstime=obs_time)
    bjd_tdb = obs_time.tdb + obs_time.light_travel_time(sky_coordinates)
    return barycorr_rv.to('m/s').value, bjd_tdb.to_value('jd')


def cross_correlate_over_traces(image, orders_to_use, velocities, template):
    ccfs = []

    for i in orders_to_use:
        logger.info(f'Cross correlating for order: {i}', image=image)
        order = image.spectrum[image.science_fiber, i]

        # Only pass in the given wavelength range +- 1 Angstrom to boost performance
        relevant_region = np.logical_and(template['wavelength'] >= np.min(order['wavelength']) - 1.0,
                                         template['wavelength'] <= np.max(order['wavelength']) + 1.0)
        template_to_fit = {'wavelength': template['wavelength'][relevant_region], 'flux': template['flux'][relevant_region]}
        # Set the model S/N = 1000 -- from looking by eye at
        # scatter in the model and from systematic uncertainties in the model
        template_error = 1e-3 * template_to_fit['flux']
        mask = mark_features(template_to_fit['flux'], sigma=3, detector_resolution=40)
        # Mask the prohibited wavelength regions.
        for mask_region in PHOENIX_WAVELENGTHS_TO_MASK:
            mask[np.logical_and(template_to_fit['wavelength'] >= min(mask_region),
                                template_to_fit['wavelength'] <= max(mask_region))] = 1
        continuum_model = fit_polynomial(template_to_fit['flux'], template_error, x=template_to_fit['wavelength'],
                                         order=3, mask=mask)
        normalized_template = {'wavelength': template_to_fit['wavelength'],
                               'flux': template_to_fit['flux'] / continuum_model(template_to_fit['wavelength'])}
        # NOTE PHOENIX WAVELENGTHS ARE IN VACUUM
        # NRES WAVELENGTHS ARE TIED TO WHATEVER LINE LIST WAS USED (e.g. nres wavelengths will be in air if ThAr atlas air
        # was used, and they will be in vacuum if ThAr_atlas_ESO_vacuum.txt was used.).
        # AS OF Aug 27 2020, NRES WAVELENGTHS ARE IN VACUUM BECAUSE ThAr_atlas_ESO_vacuum.txt IS THE LINE LIST USED.
        x_cor = cross_correlate(velocities, order['wavelength'], order['normflux'], order['normuncertainty'],
                                normalized_template['wavelength'], normalized_template['flux'])
        ccfs.append({'order': i, 'v': velocities, 'xcor': x_cor})
    return Table(ccfs)


def calculate_rv(image, orders_to_use, template):
    # for steps in 1 km/s from -2000 to +2000 km/s
    coarse_velocities = np.arange(-2000, 2001, 1)

    coarse_ccfs = cross_correlate_over_traces(image, orders_to_use, coarse_velocities, template)

    # take the peak
    velocity_peaks = np.array([coarse_velocities[np.argmax(ccf['xcor'])] for ccf in coarse_ccfs])
    best_v = stats.sigma_clipped_mean(velocity_peaks, 3.0)
    velocities = np.arange(best_v - 2, best_v + 2 + 1e-4, 1e-3)

    ccfs = cross_correlate_over_traces(image, orders_to_use, velocities, template)

    rvs_per_order = np.array([ccf['v'][np.argmax(ccf['xcor'])] for ccf in ccfs])
    # Calculate the peak v (converting to m/s) in the spectrograph frame
    rv = stats.sigma_clipped_mean(rvs_per_order, 3.0) * 1000
    rv_err = stats.robust_standard_deviation(rvs_per_order) * 1000
    return rv, rv_err, coarse_ccfs, ccfs


class RVCalculator(Stage):
    def do_stage(self, image) -> ObservationFrame:
        if image.classification is None:
            logger.warning('No classification to use for an RV template', image=image)
            return image
        phoenix_loader = phoenix.PhoenixModelLoader(self.runtime_context.db_address)
        template = phoenix_loader.load(image.classification)
        # Pick orders near the center of the detector that have a high Signal to noise and are free of tellurics.
        orders_to_use = np.arange(self.runtime_context.MIN_ORDER_TO_CORRELATE,
                                  self.runtime_context.MAX_ORDER_TO_CORRELATE, 1)

        rv_measured, rv_err, coarse_ccfs, ccfs = calculate_rv(image, orders_to_use, template)

        # Compute the barycentric RV and observing time corrections
        rv_correction, bjd_tdb = barycentric_correction(image.dateobs, image.exptime,
                                                        image.ra, image.dec,
                                                        dbs.get_site(image.instrument.site, self.runtime_context.db_address))
        # Correct the RV per Wright & Eastman (2014) and save in the header
        rv = rv_measured + rv_correction + rv_measured * rv_correction / c
        image.meta['RV'] = rv, 'Radial Velocity in Barycentric Frame [m/s]'
        # The following assumes that the uncertainty on the barycentric correction is negligible w.r.t. that
        # on the RV measured from the CCF, which should generally be true
        image.meta['RVERR'] = rv_err, 'Radial Velocity Uncertainty [m/s]'
        image.meta['BARYCORR'] = rv_correction, 'Barycentric Correction Applied to the RV [m/s]'
        image.meta['TCORR'] = bjd_tdb, 'Exposure Mid-Time (Barycentric Julian Date)'
        image.meta['TCORVERN'] = 'astropy.time.light_travel_time', 'Time correction code version'
        # TODO: verify that the following are in fact the time corrections done by astropy;
        # this isn't completely obvious from the online documentation
        image.meta['TCORCOMP'] = 'ROMER, CLOCK', 'Time corrections done'
        image.meta['TCOREPOS'] = 'ERFA', 'Source of Earth position'
        image.meta['TCORSYST'] = 'BJD_TDB ', 'Ref. frame_timesystem of TCORR column'
        image.meta['PLEPHEM'] = solar_system_ephemeris.get(), 'Source of planetary ephemerides'

        # save the fine + coarse ccfs together
        # sort the coarse and fine ccf's
        coarse_ccfs.sort('order')
        ccfs.sort('order')
        # remove entries from the coarse ccf's that fall within the velocity range of the fine ccf's
        ccf_range = (np.min(ccfs['v'][0]), np.max(ccfs['v'][0]))
        entries_to_keep = np.logical_or(coarse_ccfs['v'][0] < min(ccf_range), coarse_ccfs['v'][0] > max(ccf_range))
        coarse_ccfs['v'] = coarse_ccfs['v'][:, entries_to_keep]
        coarse_ccfs['xcor'] = coarse_ccfs['xcor'][:, entries_to_keep]
        #
        final_ccfs = Table({'order': ccfs['order'],
                            'v': np.hstack([coarse_ccfs['v'], ccfs['v']]),
                            'xcor': np.hstack([coarse_ccfs['xcor'], ccfs['xcor']])})
        # sorting xcor by the velocity grid so that things are in order
        sort_array = np.argsort(final_ccfs['v'][0])
        final_ccfs['xcor'] = final_ccfs['xcor'][:, sort_array]
        final_ccfs['v'] = final_ccfs['v'][:, sort_array]
        image.ccf = final_ccfs
        return image
