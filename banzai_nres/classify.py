from banzai.stages import Stage
from banzai_nres.frames import NRESObservationFrame
from banzai_nres import dbs
from banzai.logs import get_logger
import warnings
from astropy import units
from astropy.time import Time
from astropy.coordinates import SkyCoord
from banzai_nres.rv import cross_correlate_over_traces, calculate_rv
from banzai_nres import phoenix
import numpy as np
from banzai.utils import import_utils
import astroquery.exceptions


logger = get_logger()


def find_object_in_catalog(image, db_address, gaia_class, simbad_class):
    """
    Find the object in external catalogs. Update the ra and dec if found. Also add an initial classification if found.
    :return:
    """
    # Assume that the equinox and input epoch are both j2000.
    coordinate = SkyCoord(ra=image.ra, dec=image.dec, unit=(units.deg, units.deg),
                          frame='icrs', pm_ra_cosdec=image.pm_ra * units.mas / units.year,
                          pm_dec=image.pm_dec * units.mas / units.year, equinox='j2000',
                          obstime=Time(2000.0, format='decimalyear'))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # 10 arcseconds should be a large enough radius to capture bright objects.
        gaia = import_utils.import_attribute(gaia_class)
        gaia_connection = gaia()
        gaia_connection.ROW_LIMIT = 200
        # Gaia uses an equinox of 2000, but epoch of 2016 for the proper motion
        transformed_coordinate = coordinate.apply_space_motion(new_obstime=Time(2016.0, format='decimalyear'))

        # Use a radius of 2.6 arcseconds which is the fiber size
        results = gaia_connection.query_object(coordinate=transformed_coordinate, radius=2.6 * units.arcsec)

    # Filter out objects fainter than r=12 and brighter than r = 5.
    # There is at least one case (gamma cas) that is in gaia but does not have a complete catalog record like proper
    # motions and effective temperatures.
    results = results[np.logical_and(results['phot_rp_mean_mag'] < 12.0, results['phot_rp_mean_mag'] > 5.0)]
    if len(results) > 0:
        # Take the brightest object in the fiber
        brightest = np.ma.argmin(results['phot_rp_mean_mag'])
        # If the entire set of results has masked magnitudes, move on to
        # simbad.
        if not np.ma.is_masked(results['phot_rp_mean_mag'][brightest]):
            image.classification = dbs.get_closest_HR_phoenix_models(
                db_address,
                results[brightest]['teff_gspphot'],
                results[brightest]['logg_gspphot']
            )
            # Update the ra and dec to the catalog coordinates as those are basically always better than a user enters
            # manually.
            # We now have to rewind the coordinates from 2016 to 2000 since
            # that is our convention
            if results[brightest]['pmra'] is not np.ma.masked:
                image.pm_ra, image.pm_dec = results[brightest]['pmra'], results[brightest]['pmdec']
            gaia_coordinate = SkyCoord(
                ra=results[brightest]['ra'],
                dec=results[brightest]['dec'],
                unit=(units.deg, units.deg),
                frame='icrs',
                pm_ra_cosdec=image.pm_ra * units.mas / units.year,
                pm_dec=image.pm_dec * units.mas / units.year, equinox='j2000',
                obstime=Time(2016.0, format='decimalyear')
            )
            epoch2000 = Time(2000.0, format='decimalyear')
            transformed_gaia_coordinate = gaia_coordinate.apply_space_motion(new_obstime=epoch2000)
            image.ra, image.dec = transformed_gaia_coordinate.ra.deg, transformed_gaia_coordinate.dec.deg

    # If nothing in Gaia fall back to simbad. This should only be for stars that are brighter than mag = 3
    else:
        # IMPORTANT NOTE:
        # During e2e tests we do not import astroquery.simbad.Simbad. We import a mocked simbad call
        # which can be found in banzai_nres.tests.utils.MockSimbad . This returns a simbad table that is
        # truncated. If you add a new votable field, you will need to add it to the mocked table as well.
        simbad = import_utils.import_attribute(simbad_class)
        simbad_connection = simbad()
        simbad_connection.add_votable_fields('pmra', 'pmdec', 'mesfe_h', 'otype', 'V', 'r')
        try:
            # Use a radius of 2.6 arcseconds which is the fiber size
            results = simbad_connection.query_region(coordinate, radius='0d0m2.6s')
        except astroquery.exceptions.TableParseError:
            response = simbad_connection.last_response.content
            logger.error(f'Error querying SIMBAD. Response from SIMBAD: {response}', image=image)
            results = []
        if results:
            results = remove_planets_from_simbad(results)
            if len(results) == 0:
                logger.warning('All objects in SIMBAD query are planets. Will not classify.', image=image)
                return
            # Get the brightest object in the fiber
            brightest = np.ma.argmin([row['V'] for row in results])
            if results[brightest]['V'] is np.ma.masked:
                # Try r band if there no V band
                brightest = np.ma.argmin([row['r'] for row in results])
                if results[brightest]['r'] is np.ma.masked:
                    logger.warning(
                        'All objects in SIMBAD query are missing V and r magnitudes. '
                        'Cannot determine which is brightest. Will not classify.',
                        image=image
                    )
                    return
            results = results[brightest]  # get the brightest source.

            phoenix_models = dbs.get_closest_phoenix_models(
                db_address,
                results['mesfe_h.teff'],
                results['mesfe_h.log_g']
            )
            image.classification = phoenix_models[0]
            # note that we always assume the proper motions are in mas/yr... which they should be.
            # And that they include cos(dec)
            if results['pmra'] is not np.ma.masked:
                image.pm_ra, image.pm_dec = results['pmra'], results['pmdec']
            # Update the ra and dec to the catalog coordinates as those will be consistent across observations.
            # The new version of astroquery returns ra and dec in degrees, so we can use them directly.
            coord = SkyCoord(results['ra'], results['dec'], unit=(units.deg, units.deg))
            image.ra, image.dec = coord.ra.deg, coord.dec.deg
        # If there are still no results, then do nothing


def remove_planets_from_simbad(results):
    # Remove planets. See https://simbad.u-strasbg.fr/simbad/sim-display?data=otypes for otype designations.
    return results[['Pl' not in row['otype'] for row in results]]


class StellarClassifier(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        # Short circuit if the header coordinates are malformed
        if any(np.isnan([image.ra, image.dec, image.pm_ra, image.pm_dec])):
            if len(image.meta['RADESYS'].strip()) == 0:
                # Assume that we are looking at a non-sidereal object
                pass
            else:
                logger.error('RA or Dec or Proper Motion is malformed in header.', image=image)
            return image
        closest_previous_classification = dbs.get_closest_existing_classification(self.runtime_context.db_address,
                                                                                  image.ra, image.dec)
        if closest_previous_classification is not None:
            previous_coordinate = SkyCoord(closest_previous_classification.ra, closest_previous_classification.dec,
                                           unit=(units.deg, units.deg))
            this_coordinate = SkyCoord(image.ra, image.dec, unit=(units.deg, units.deg))

            # Short circuit if the object is already classified
            # We choose 2.6 arcseconds as the don't reclassify cutoff radius as it is the fiber size
            if this_coordinate.separation(previous_coordinate) < 2.6 * units.arcsec:
                image.classification = dbs.get_phoenix_model_by_id(closest_previous_classification.phoenix_id,
                                                                   self.runtime_context.db_address,)
                image.meta['CLASSIFY'] = 0, 'Was this spectrum classified'
                return image

        # The object was not previously classified, check Gaia and/or SIMBAD
        find_object_in_catalog(image, self.runtime_context.db_address,
                               self.runtime_context.GAIA_CLASS, self.runtime_context.SIMBAD_CLASS)

        # Short circuit if the object is not a star
        if image.classification is None:
            image.meta['CLASSIFY'] = 0, 'Was this spectrum classified'
            logger.warning('Target not found in SIMBAD or Gaia. Assuming it is not a star. Will not classify',
                           image=image)
            return image
        orders_to_use = np.arange(self.runtime_context.MIN_ORDER_TO_CORRELATE,
                                  self.runtime_context.MAX_ORDER_TO_CORRELATE, 1)
        phoenix_loader = phoenix.PhoenixModelLoader(self.runtime_context)
        template = phoenix_loader.load(image.classification)
        rv = calculate_rv(image, orders_to_use, template)[0]
        best_metric = np.sum(cross_correlate_over_traces(image, orders_to_use, np.array([rv.value]) * rv.unit,
                                                         template)['xcor'])
        # For each param: Fix the other params, get the N closest models and save the results
        physical_parameters = ['T_effective', 'log_g', 'metallicity', 'alpha']
        n_steps = [11, 5, 5, 5]
        for parameter, n in zip(physical_parameters, n_steps):
            models_to_test = dbs.get_closest_phoenix_models(self.runtime_context.db_address,
                                                            image.classification.T_effective,
                                                            image.classification.log_g,
                                                            image.classification.metallicity,
                                                            image.classification.alpha,
                                                            fixed=[i for i in physical_parameters if i != parameter],
                                                            n=n)
            for model_to_test in models_to_test:
                template = phoenix_loader.load(model_to_test)
                metric = np.sum(cross_correlate_over_traces(image, orders_to_use, np.array([rv.value]) * rv.unit,
                                                            template)['xcor'])
                if metric > best_metric:
                    image.classification = model_to_test
                    best_metric = metric
        if image.classification is None:
            image.meta['CLASSIFY'] = 0, 'Was this spectrum classified'
        else:
            image.meta['CLASSIFY'] = 1, 'Was this spectrum classified'
            dbs.save_classification(self.runtime_context.db_address, image)

        return image
