import numpy as np


def find_nearest(array, value):
    return array[(np.abs(array - value)).argmin()]