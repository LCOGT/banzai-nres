import numpy as np
from photutils import DAOStarFinder
from astropy.table import Table

import logging

logger = logging.getLogger('banzai')


def identify_features(data, err, mask=None, nsigma=2., fwhm=6.0, **kwargs):
    """
    :param data: 2d ndarray. The data to identify features (e.g. ThAr emission lines) on.
    :param err: 2d ndarray. The standard 1-sigma error for data.
    :param mask: 2b ndarray, the same shape as data. Masked values take a value of 1.
    :param nsigma: minimum multiple above the read noise (calculated via np.min(err)) for a feature to be cataloged.
    :param fwhm: The full-width half-maximum (FWHM) of the major axis of the Gaussian kernel in units of pixels.
                 Features with fwhm below this will not be counted.
    :param kwargs: extra key word arguments to pass to photutils.DAOStarFinder, e.g. theta .
    :return: features. Table.
             catalog of features (e.g. emission lines). features['xcentroid'][j] gives the horizontal pixel position
             of the jth feature. features['ycentroid'][j] gives the y pixel (vertical) position of the jth feature.
    """
    daofind = DAOStarFinder(fwhm=fwhm, threshold=nsigma, exclude_border=True, **kwargs)
    features = daofind(data / err, mask=mask)
    if features is None:
        features = Table({'xcentroid': [], 'ycentroid': [], 'flux': []})
    features['pixel'] = features['xcentroid']  # because xwavecal uses 'pixel' as the coordinate key.
    return features


def index_of_refraction_Ciddor1996(vacuum_wavelength):
    """
    From Ciddor 1996, Also (copied from Ciddor 1996): See T.-O. Husser et al.: A new extensive library of PHOENIX stellar
    atmospheres and synthetic spectra: https://www.aanda.org/articles/aa/pdf/2013/05/aa19058-12.pdf. Equation 9 and 10.
    :param vacuum_wavelength: wavelength in vacuum
    :return: index of refraction at that wavelength.
    Comments:
    Standard air: dry air at 15 °C, 101.325 kPa and with 450 ppm CO2 content.
    wavelength_vacuum: true
    temperature: 15 °C
    pressure: 101325 Pa
    """
    sig2 = (1E4/vacuum_wavelength)**2
    return 1 + 0.05792105/(238.0185 - sig2) + 0.00167917/(57.362 - sig2)


def index_of_refraction_Edlen(vacuum_wavelength):
    """
    The original 1966 Edlen Equation (Bengt Edlén 1966 Metrologia 2 71, https://doi.org/10.1088/0026-1394/2/2/002)
    :param vacuum_wavelength: wavelength in vacuum
    :return: index of refraction at that wavelength.
    """
    sig2 = (1E4/vacuum_wavelength)**2
    return 1 + 1E-8 * (8342.13 + 2406030/(130 - sig2) + 15997/(38.9 - sig2))


def group_features_by_trace(features, traces):
    """
    :return: features.
             where features['id'][j] gives the id number of the trace that the feature falls in.
    """
    features['id'] = traces[np.array(features['ycentroid'], dtype=int), np.array(features['xcentroid'], dtype=int)]
    return features


def get_center_wavelengths(features):
    ref_ids = np.sort(np.array(list(set(features['order']))))
    wavelengths = [np.average(features['wavelength'][features['order'] == o]) for o in ref_ids]
    return wavelengths


def get_principle_order_number(m0_values, features):
    """
    Finds the principle order number m0. Selects the m0 such that the function y(i) = (m0 + i) * central_wavelengths
    has the smallest slope. I.e. this selects the m0 that allows constant/(m0+i) to best fit central_wavelengths.
    This is exactly what CERES does. See equation 3 of Brahm et al. 2016.

    :param m0_values: ndarray of integers. 1d.
    :param features: dict.
           dictionary of ndarrays of pixel and order positions of measured spectral features (e.g. emission lines)
           and their wavelengths.
           Example:
               measured_lines = {'pixel': np.array([1, 2.5, 6.1]), 'order': np.array([1, 1, 2]),
                                 'wavelength': np.array([4000, 5001, 5005)}
               If the principle order number is 52, then these measured_lines represents 3 spectral features,
               with (pixel, diffraction order) coordinates of (1, 53), (2.5, 53), (6.1, 54), and wavelengths
               4000, 5001, and 5005 Angstroms, respectively.
               respectively. The wavelength solution will calibrate each fiber separately.
    :return: m0: int.
    The principle order number for the fiber from which central_wavelengths were taken.
    This is the true order index of the the trace that corresponds to ref_id[0].
    """
    center_wavelengths = get_center_wavelengths(features)
    ref_ids = np.sort(np.array(list(set(features['order']))))
    slopes = []
    for m0 in m0_values:
        # note: replacing np.ptp with some outlier resistant measure of the scatter would be more robust.
        slopes.append(np.ptp(center_wavelengths * (m0 + ref_ids)))

    if np.count_nonzero(np.isclose(slopes, np.min(slopes), rtol=0.01)) > 1:
        logger.warning('Two or more viable principle order numbers for this fiber! The m0 recovered from the '
                       'wavelength solution could be wrong!')
    return m0_values[np.argmin(slopes)]
