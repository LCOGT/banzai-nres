from banzai.stages import Stage
from banzai_nres.frames import NRESObservationFrame
from banzai_nres import dbs
import warnings
from astropy import units
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
simbad = Simbad()
simbad.add_votable_fields('pmra', 'pmdec', 'fe_h')
Gaia.ROW_LIMIT = 200


def get_initial_guess(ra, dec, pm_ra, pm_dec):
    """

    :param ra: hour angle
    :param dec: deg
    :param pm_ra: arcsec/yr, Note this needs to include the cos dec term like it is in simbad
    :param pm_dec: arcsec/yr
    :return:
    """

    # Assume that the equinox and input epoch are both j2000.
    # Gaia uses an equinox of 2000, but epoch of 2015.5 for the proper motion
    coordinate = SkyCoord(ra=ra, dec=dec, unit=(units.hourangle, units.deg),
                          frame='icrs', pm_ra_cosdec=pm_ra * units.arcsec / units.year,
                          pm_dec=pm_dec * units.arcsec / units.year, equinox='j2000',
                          obstime=Time(2000.0, format='decimalyear'))
    transformed_coordinate = coordinate.apply_space_motion(new_obstime=Time(2015.5, format='decimalyear'))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # 10 arcseconds should be a large enough radius to capture bright objects.
        results = Gaia.query_object(coordinate=transformed_coordinate, radius=10.0 * units.arcsec)

    # Filter out objects fainter than r=12
    results = results[results['phot_rp_mean_mag'] < 12.0]
    # If nothing in Gaia fall back to simbad. This should only be for stars that are brighter than mag = 3
    if len(results) == 0:
        results = simbad.query_region(coordinate, radius='0d0m10s')
        results.rename_column('Fe_H_log_g', 'log_g')
        # If there are still no results, then abort
        if len(results) == 0:
            return None
        else:
            return results[0]['RA', 'DEC', 'PMRA', 'PMDEC', 'Fe_H_Teff', 'log_g']

    else:
        return results[0]['ra', 'dec', 'pmra', 'pmdec', 'teff_val', 'lum_val']


class StellarClassifier(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        initial_guess = get_initial_guess(image.ra, image.dec, image.pm_ra, image.pm_dec)
        if initial_guess is None:
            return image

        if 'log_g' in initial_guess.colnames:
            ra, dec, pm_ra, pm_dec, T_effective, log_g = initial_guess
            # Assume solar alpha abundance and metallicity to start
            image.classification = dbs.get_closest_phoenix_models(self.runtime_context.db_address, T_effective, log_g)
        # Assume HR diagram
        else:
            ra, dec, pm_ra, pm_dec, T_effective, luminosity = initial_guess
            image.classification = dbs.get_closest_HR_phoenix_models(self.runtime_context.db_address, T_effective, luminosity)
        # Update the ra and dec to the catalog coordinates as those are basically always better than a user enters
        # manually.
        image.ra, image.dec = ra, dec
        image.pm_ra, image.pm_dec = pm_ra, pm_dec
        # TODO: For each param: Fix the other params, get the N closest models and save the results
        return image
