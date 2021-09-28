from astropy.table import Table
from banzai_nres.classify import parse_simbad_coordinates, restrict_simbad_results_to_stellar_only
from banzai_nres.classify import convert_simbad_coordinates_to_degrees, get_closest_source
from astropy import units
from astropy.coordinates import SkyCoord
import pkg_resources
import numpy as np

SIMBAD_RESPONSE_FILENAME = pkg_resources.resource_filename('banzai_nres.tests', 'data/mock_simbad_response.ecsv')


def test_parse_simbad_coordinates():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    assert type(results['RA'].unit) is units.UnrecognizedUnit
    results = parse_simbad_coordinates(results)
    assert results['RA'].unit == units.hourangle
    assert results['DEC'].unit == units.deg


def test_restrict_simbad_results_to_stellar_only():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    results = restrict_simbad_results_to_stellar_only(results)
    assert len(results) == 1
    assert results['MAIN_ID'] == '* tau Cet'


def test_restrict_simbad_results_to_stellar_only_improper():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    results['OTYPE'][results['MAIN_ID'] == '* tau Cet'] = 'SN*'
    results = restrict_simbad_results_to_stellar_only(results)
    assert len(results) == 0


def test_convert_simbad_coordinates_to_degrees():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    results = parse_simbad_coordinates(results)
    results = convert_simbad_coordinates_to_degrees(results)
    results = results[results['MAIN_ID'] == '* tau Cet']
    assert np.isclose(results[0]['RA'], 26.017014166666666)
    assert np.isclose(results[0]['DEC'], -15.937479444444445)


def test_convert_simbad_coordinates_to_degrees_improper():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    results['RA'].unit = units.UnrecognizedUnit('Incomprehensible Unit')
    results = convert_simbad_coordinates_to_degrees(results)
    assert results is None


def test_get_closest_source():
    results = Table({'NAME': ['wrong', 'right'],
                     'RA': [0, 1], 'DEC': [0, 0]},
                    units=('', units.deg, units.deg))
    coordinate = SkyCoord(0.9, 0, unit=(units.deg, units.deg))
    results = get_closest_source(results, coordinate)
    assert results['NAME'] == 'right'


def test_integration_on_no_match():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    results['OTYPE'][results['MAIN_ID'] == '* tau Cet'] = 'SN*'
    results = restrict_simbad_results_to_stellar_only(results)
    results = parse_simbad_coordinates(results)
    results = convert_simbad_coordinates_to_degrees(results)
    assert results is None


def test_integration_on_match():
    results = Table.read(SIMBAD_RESPONSE_FILENAME)
    results = restrict_simbad_results_to_stellar_only(results)
    results = parse_simbad_coordinates(results)
    results = convert_simbad_coordinates_to_degrees(results)
    results = get_closest_source(results, SkyCoord("01 44 04.0834", "-15 56 14.926", unit=(units.hourangle, units.deg)))
    assert results['MAIN_ID'] == '* tau Cet'
