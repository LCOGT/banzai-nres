from astropy.table import Table
from banzai_nres.classify import remove_planets_from_simbad


def test_remove_planets_from_simbad():
    results = Table({'MAIN_ID': ['First', 'Second'], 'otype': ['Star', 'Planet']})
    results = remove_planets_from_simbad(results)
    assert len(results) == 1
    assert results[0]['MAIN_ID'] == 'First'
