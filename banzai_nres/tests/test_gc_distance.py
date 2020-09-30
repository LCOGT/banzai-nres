from banzai_nres.dbs import great_circle_distance
from astropy.coordinates import SkyCoord
from astropy import units
import math
import numpy as np


def test_gc_distance():
    ra1, dec1 = 150.0, 25.0
    ra2, dec2 = 100.0, 10.0
    coord1 = SkyCoord(ra1, dec1, unit=(units.deg, units.deg))
    coord2 = SkyCoord(ra2, dec2, unit=(units.deg, units.deg))
    expected = coord1.separation(coord2).deg

    actual = great_circle_distance(ra1, dec1, ra2, dec2, math)
    np.testing.assert_allclose(actual, expected)
