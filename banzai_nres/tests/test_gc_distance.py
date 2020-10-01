from banzai_nres.dbs import cos_great_circle_distance
from astropy.coordinates import SkyCoord
from astropy import units
import numpy as np


def test_gc_distance():
    ra1, dec1 = 150.0, 25.0
    ra2, dec2 = 100.0, 10.0
    coord1 = SkyCoord(ra1, dec1, unit=(units.deg, units.deg))
    coord2 = SkyCoord(ra2, dec2, unit=(units.deg, units.deg))
    expected = np.cos(np.deg2rad(coord1.separation(coord2).deg))

    actual = cos_great_circle_distance(np.sin(np.deg2rad(ra1)), np.cos(np.deg2rad(ra1)),
                                       np.sin(np.deg2rad(dec1)), np.cos(np.deg2rad(dec1)),
                                       np.sin(np.deg2rad(ra2)), np.cos(np.deg2rad(ra2)),
                                       np.sin(np.deg2rad(dec2)), np.cos(np.deg2rad(dec2)))
    np.testing.assert_allclose(actual, expected)
