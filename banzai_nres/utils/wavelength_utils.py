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
    return features


def group_features_by_trace(features, traces):
    """
    :return: features.
             where features['id'][j] gives the id number of the trace that the feature falls in.
    """
    features['id'] = traces[np.array(features['ycentroid'], dtype=int), np.array(features['xcentroid'], dtype=int)]
    return features


def aperture_extract(x, y, data, aperture_width, mask=None):
    """
    Extract all flux within the aperture of (width, height) = (aperture_width, aperture_width)
    around each point given by x and y.
    :param x: ndarray, column indices of each extraction point in data.
    :param y: ndarray, row indices of each extraction point in data.
    :param data: 2d ndarray. data where (x[i], y[i]) is the coordinate of the ith source of interest.
    :param mask: 2d ndarray. Same shape as data. Masked values take a value of 1.
    :param aperture_width: int. The width (and height) of the box around each point (x[i], y[i]) where flux will be summed.
    :return: flux. ndarray. Same shape as x and y. Summed flux around each point x and y, in order of the points.
    """
    if mask is None:
        mask = np.zeros_like(data)
    # ignore masked points in the sum
    masked_data = np.zeros_like(data, dtype=float)
    masked_data[np.isclose(mask, 0)] = data[np.isclose(mask, 0)]
    # get the set of coordinates to extract
    hw = max(int(aperture_width/2), 1)
    xoffsets, yoffsets = np.meshgrid(np.arange(-hw, hw+1), np.arange(-hw, hw+1))
    #
    fluxes = np.zeros_like(x, dtype=float)
    x_int, y_int = np.array(x, dtype=int), np.array(y, dtype=int)
    for deltax, deltay in zip(xoffsets.ravel(), yoffsets.ravel()):
        fluxes += data[y_int + deltay, x_int + deltax]
    return fluxes
