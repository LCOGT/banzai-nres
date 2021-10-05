from astropy.table import Table
import numpy as np


class MockGaiaCatalog:
    def query_object(self, *args, **kwargs):
        return Table({'phot_rp_mean_mag': [21.0]})


class MockSimbad:
    def add_votable_fields(self, *args, **kwargs):
        pass

    def query_region(self, *args, **kwarg):
        data = {
            'MAIN_ID': ['* alf CMi', '* alf CMi B'],
            'RA': ['07 39 18.1195', '07 39 17.8800'],
            'DEC': ['+05 13 29.955', '+05 13 26.800'],
            'COO_BIBCODE': ['2007A&A...474..653V', '2015ApJS..219...19L'],
            'PMRA': [-714.59, -709.0],
            'PMDEC': [-1036.8, -1024.0],
            'Fe_H_Teff': [6474, 7870],
            'Fe_H_log_g': [3.964200019836426, 8.079999923706055],
            'Fe_H_Fe_H': [0.009999999776482582, np.nan],
            'Fe_H_flag': ['', ''],
            'Fe_H_bibcode': ['2019A&A...624A..19B', '2015ApJS..219...19L'],
            'OTYPE': ['SB*', 'WD*'],
            'SCRIPT_NUMBER_ID': [1, 1],
            }
        return Table(data)
