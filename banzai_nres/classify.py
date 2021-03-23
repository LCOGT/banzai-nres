from banzai.stages import Stage
from banzai_nres.frames import NRESObservationFrame
from banzai_nres import dbs
import warnings
from astropy import units
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import constants
from banzai_nres.rv import cross_correlate_over_traces, calculate_rv
from banzai_nres import phoenix
import numpy as np
from banzai.utils import import_utils


def find_object_in_catalog(image, db_address, gaia_class, simbad_class):
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
        gaia = import_utils.import_attribute(gaia_class)
        gaia_connection = gaia()
        gaia_connection.ROW_LIMIT = 200
        results = gaia_connection.query_object(coordinate=transformed_coordinate, radius=10.0 * units.arcsec)

    # Filter out objects fainter than r=12
    results = results[results['phot_rp_mean_mag'] < 12.0]
    if len(results) > 0:
        # convert the luminosity from the LSun units that Gaia provides to cgs units
        results[0]['lum_val'] *= constants.L_sun.to('erg / s').value
        image.classification = dbs.get_closest_HR_phoenix_models(db_address, results[0]['teff_val'],
                                                                 results[0]['lum_val'])
        # Update the ra and dec to the catalog coordinates as those are basically always better than a user enters
        # manually.
        image.ra, image.dec = results[0]['ra'], results[0]['dec']
        image.pm_ra, image.pm_dec = results[0]['pmra'], results[0]['pmdec']
    # If nothing in Gaia fall back to simbad. This should only be for stars that are brighter than mag = 3
    else:
        simbad = import_utils.import_attribute(simbad_class)
        simbad_connection = simbad()
        simbad_connection.add_votable_fields('pmra', 'pmdec', 'fe_h')
        results = simbad_connection.query_region(coordinate, radius='0d0m10s')
        if results:
            image.classification = dbs.get_closest_phoenix_models(db_address, results[0]['Fe_H_Teff'],
                                                                  results[0]['Fe_H_log_g'])[0]
            # Update the ra and dec to the catalog coordinates as those are basically always better than a user enters
            # manually.
            image.ra, image.dec = results[0]['RA'], results[0]['DEC']
            image.pm_ra, image.pm_dec = results[0]['PMRA'], results[0]['PMDEC']
        # If there are still no results, then do nothing


class StellarClassifier(Stage):
    def do_stage(self, image) -> NRESObservationFrame:
        find_object_in_catalog(image, self.runtime_context.db_address,
                               self.runtime_context.GAIA_CLASS, self.runtime_context.SIMBAD_CLASS)

        closest_previous_classification = dbs.get_closest_existing_classification(self.runtime_context.db_address,
                                                                                  image.ra, image.dec)

        if closest_previous_classification is not None:
            previous_coordinate = SkyCoord(closest_previous_classification.ra, closest_previous_classification.dec,
                                           unit=(units.deg, units.deg))
            this_coordinate = SkyCoord(image.ra, image.dec, unit=(units.deg, units.deg))

            # Short circuit if the object is already classified
            # We choose 2.6 arcseconds as the don't reclassify cutoff radius as it is the fiber size
            if this_coordinate.separation(previous_coordinate) < 2.6 * units.arcsec:
                image.classification = closest_previous_classification.phoenix_model
                image.meta['CLASSIFY'] = 0, 'Was this spectrum classified'
                return image


        orders_to_use = np.arange(self.runtime_context.MIN_ORDER_TO_CORRELATE,
                                  self.runtime_context.MAX_ORDER_TO_CORRELATE, 1)
        phoenix_loader = phoenix.PhoenixModelLoader(self.runtime_context.db_address)
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
        image.meta['TEFF'] = image.classification.T_effective, 'Estimated stellar effective temperature [K]'
        image.meta['LOGG'] = image.classification.log_g, 'Estimated stellar surface gravity [cgs]'
        image.meta['FEH'] = image.classification.metallicity, 'Estimated stellar metallicity [dex]'
        image.meta['ALPHA'] = image.classification.alpha, 'Estimated stellar alpha abundance [dex]'

        return image
