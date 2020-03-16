import numpy as np
from scipy.interpolate import SmoothBivariateSpline, UnivariateSpline
from statsmodels.robust.norms import HuberT


def fit_smooth_spline(y, error, mask=None, n_iter=5, x=None, sigma=5):
    if mask is None:
        mask = np.zeros(y.shape, dtype=bool)
    y_to_fit = y[np.logical_not(mask)]
    error_to_fit = error[np.logical_not(mask)]
    # Note that from the docs, the weights should be about 1 over the standard deviation for the smoothing to
    # work well
    weights = 1.0 / error_to_fit
    huber_t = HuberT(t=sigma)
    if len(y.shape) == 2 and x is None:
        spline = SmoothBivariateSpline
        x1d, y1d = np.arange(y.shape[0], dtype=np.float), np.arange(y.shape[1], dtype=np.float)
        x2d, y2d = np.meshgrid(x1d, y1d)
        x_to_fit = x2d.T[np.logical_not(mask)], y2d.T[np.logical_not(mask)]
        kwargs = {'grid': False}

    elif len(y.shape) == 2 or len(x) == 2:
        spline = SmoothBivariateSpline
        x_to_fit = (x[0][np.logical_not(mask)], x[1][np.logical_not(mask)])
        kwargs = {'grid': False}

    elif x is not None:
        spline = UnivariateSpline
        x_to_fit = (x[np.logical_not(mask)], )
        kwargs = {}
    else:
        x_to_fit = (x[0][np.logical_not(mask)],)
        spline = UnivariateSpline
        kwargs = {}

    for i in range(n_iter):
        best_fit = spline(*x_to_fit, y_to_fit, w=weights)

        # Iteratively reweight to reject outliers
        weights = huber_t.weights(np.abs(y_to_fit - best_fit(*x_to_fit, **kwargs)) / error_to_fit) / error_to_fit
    return best_fit
