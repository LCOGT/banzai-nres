import numpy as np
from photutils import DAOStarFinder
from astropy.table import Table


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
    daofind = DAOStarFinder(fwhm=fwhm, threshold=np.nanmin(err) * nsigma, exclude_border=True, **kwargs)
    features = daofind(data, mask=mask)
    if features is None:
        features = Table({'xcentroid': [], 'ycentroid': [], 'flux': []})
    features['pixel'] = features['xcentroid']  # because xwavecal uses 'pixel' as the coordinate key.

    # TODO DAOStarFind uses flux = peak/threshold. We want the total flux in the line. So features['flux'] is not
    #  yet appropriate. Expect the overlap fitting to fail.
    #features['flux'] = sum_flux(features['xcentroid'], ['ycentroid'], data, mask, fwhm=fwhm)
    return features


def group_features_by_trace(features, traces):
    """
    :return: features.
             where features['id'][j] gives the id number of the trace that the feature falls in.
    """
    features['id'] = traces[features['ycentroid'].astype(int), features['xcentroid'].astype(int)]
    return features

