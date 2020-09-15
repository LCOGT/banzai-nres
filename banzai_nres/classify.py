from banzai.stages import Stage
from banzai_nres.frames import NRESObservationFrame
from banzai_nres import dbs
import warnings
from astropy import units
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroquery import gaia, simbad


def find_object_in_catalog(image, db_address):
    """
    Find the object in external catalogs. Update the ra and dec if found. Also add an initial classification if found.
    :return:
    """

    # Assume that the equinox and input epoch are both j2000.
    # Gaia uses an equinox of 2000, but epoch of 2015.5 for the proper motion
    coordinate = SkyCoord(ra=image.ra, dec=image.dec, unit=(units.deg, units.deg),
                          frame='icrs', pm_ra_cosdec=image.pm_ra * units.mas / units.year,
                          pm_dec=image.pm_dec * units.mas / units.year, equinox='j2000',
                          obstime=Time(2000.0, format='decimalyear'))
    transformed_coordinate = coordinate.apply_space_motion(new_obstime=Time(2015.5, format='decimalyear'))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # 10 arcseconds should be a large enough radius to capture bright objects.
        gaia_connection = gaia.GaiaClass()
        gaia_connection.ROW_LIMIT = 200
        results = gaia_connection.query_object(coordinate=transformed_coordinate, radius=10.0 * units.arcsec)

    # Filter out objects fainter than r=12
    results = results[results['phot_rp_mean_mag'] < 12.0]
    if len(results) > 0:
        image.classification = dbs.get_closest_HR_phoenix_models(db_address, results[0]['teff_val'],
                                                                 results[0]['lum_val'])
        # Update the ra and dec to the catalog coordinates as those are basically always better than a user enters
        # manually.
        image.ra, image.dec = results[0]['ra'], results[0]['dec']
        image.pm_ra, image.pm_dec = results[0]['pmra'], results[0]['pmdec']
    # If nothing in Gaia fall back to simbad. This should only be for stars that are brighter than mag = 3
    else:
        simbad_connection = simbad.Simbad()
        simbad_connection.add_votable_fields('pmra', 'pmdec', 'fe_h')
        results = simbad_connection.query_region(coordinate, radius='0d0m10s')
        if results:
            image.classification = dbs.get_closest_phoenix_models(db_address, results[0]['Fe_H_Teff'],
                                                                  results[0]['Fe_H_log_g'])
            # Update the ra and dec to the catalog coordinates as those are basically always better than a user enters
            # manually.
            image.ra, image.dec = results[0]['RA'], results[0]['DEC']
            image.pm_ra, image.pm_dec = results[0]['PMRA'], results[0]['PMDEC']
        # If there are still no results, then do nothing


class StellarClassifier(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        find_object_in_catalog(image, self.runtime_context.db_address)

        # TODO: For each param: Fix the other params, get the N closest models and save the results

        if image.classification is None:
            image.meta['CLASSIFY'] = 0, 'Was this spectrum classified'
        else:
            image.meta['CLASSIFY'] = 1, 'Was this spectrum classified'
            image.meta['TEFF'] = image.classification.T_effective
            image.meta['LOG_G'] = image.classification.log_g
            image.meta['FE_H'] = image.classification.metallicity
            image.meta['ALPHA'] = image.classification.alpha

        return image
