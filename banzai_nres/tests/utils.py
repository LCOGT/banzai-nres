from astropy.table import Table
import numpy as np


class MockGaiaCatalog:
    def query_object(self, *args, **kwargs):
        return Table({'phot_rp_mean_mag': [21.0]})


class MockSimbad:
    def add_votable_fields(self, *args, **kwargs):
        pass

    def query_region(self, *args, **kwarg):
        data = {'RA': ['07 39 18.1195', '07 39 17.8800'],
                'DEC': ['+05 13 29.955', '+05 13 26.800'],
                'PMRA': [-714.59,  -709.],
                'PMDEC': [-1036.8, -1024.],
                'Fe_H_Teff': [6474, 7870],
                'Fe_H_log_g': [3.9642, 8.08],
                'Fe_H_Fe_H': [0.01, np.nan]}
        return Table(data)
