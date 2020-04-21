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


def get_center_wavelengths(wavelengths, traces, trace_ids):
    # get the center wavelength from each trace in trace_ids.
    cc = wavelengths.shape[1] // 2
    center_wavelengths = [wavelengths[:, cc][traces[:, cc] == i].flatten()[0] for i in trace_ids]
    return center_wavelengths


def get_principle_order_number(m0_values, wavelengths, traces, trace_ids, ref_ids):
    """
    Finds the principle order number m0. Selects the m0 such that the function y(i) = (m0 + i) * central_wavelengths
    has the smallest slope. I.e. this selects the m0 that allows constant/(m0+i) to best fit central_wavelengths.
    This is exactly what CERES does. See equation 3 of Brahm et al. 2016.

    :param m0_values: ndarray of integers. 1d.
    :param wavelengths: ndarray of wavelengths. The wavelengths at every pixel in the image.
    :param traces:
    :param trace_ids:
    :param ref_ids: ndarray of integers. The reference id's of the traces from which central_wavelengths came.
    Same shape as central_wavelengths. i.e. central_wavelengths[0] is the center wavelength of the trace with reference
    id ref_ids[0]
    :return: m0: int.
    The principle order number for the fiber from which central_wavelengths were taken.
    This is the true order index of the the trace that corresponds to ref_id[0].
    """
    center_wavelengths = get_center_wavelengths(wavelengths, traces, trace_ids)
    slopes = []
    for m0 in m0_values:
        # note: replacing np.ptp with some outlier resistant measure of the scatter would be more robust.
        slopes.append(np.ptp(center_wavelengths * (m0 + ref_ids)))

    if np.count_nonzero(np.isclose(slopes, np.min(slopes), rtol=0.01)) > 1:
        logger.warning('Two or more viable principle order numbers for this fiber! The m0 recovered from the '
                       'wavelength solution could be wrong!')
    return m0_values[np.argmin(slopes)]