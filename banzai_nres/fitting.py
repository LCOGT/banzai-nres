import numpy as np
from scipy.interpolate import SmoothBivariateSpline, UnivariateSpline
from statsmodels.robust.norms import HuberT


def fit_smooth_spline(y, error, mask=None, n_iter=5, x=None, sigma=5):
    if mask is None:
        mask = np.zeros(y.shape, dtype=bool)
    y_to_fit = y[~mask]
    error_to_fit = error[~mask]
    # Note that from the docs, the weights should be about 1 over the standard deviation for the smoothing to
    # work well
    weights = 1.0 / error_to_fit
    huber_t = HuberT(t=sigma)

    if len(y.shape) == 2 and x is None:
        spline = SmoothBivariateSpline
        x1d, y1d = np.arange(y.shape[0], dtype=np.float), np.arange(y.shape[1], dtype=np.float)
        x2d, y2d = np.meshgrid(x1d, y1d)
        x_to_fit = x2d.T[~mask], y2d.T[~mask]
        x_eval = x1d, y1d
    elif len(y.shape) == 2:
        spline = SmoothBivariateSpline
        x_to_fit = (x[0][~mask], x[1][~mask])
        x_eval = x
    elif x is not None:
        spline = UnivariateSpline
        x_to_fit = (x[~mask], )
        x_eval = (x, )
    else:
        spline = UnivariateSpline
        x_eval = (np.arange(len(y), dtype=np.float), )
        x_to_fit = (x[0][~mask],)
        spline = UnivariateSpline

    for i in range(n_iter):
        best_fit = spline(*x_to_fit, y_to_fit, w=weights)

        # Iteratively reweight to reject outliers
        weights = huber_t.weights(np.abs(y_to_fit - best_fit(*x_eval)[~mask]) / error_to_fit) / error_to_fit
    return best_fit
