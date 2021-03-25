import numpy as np
from statsmodels.robust.norms import HuberT
from astropy.modeling import fitting, polynomial


def fit_polynomial(y, error, mask=None, n_iter=5, x=None, sigma=5, order=3, domain=None):
    if mask is None:
        mask = np.zeros(y.shape, dtype=bool)

    if x is None:
        x = np.arange(len(y), dtype=float)

    y_to_fit = y[np.logical_not(mask)]
    error_to_fit = error[np.logical_not(mask)]
    x_to_fit = x[np.logical_not(mask)]
    # start with inverse variance weights
    weights = error_to_fit ** -2.0
    weight_class = HuberT(t=sigma)
    # instantiate fitter
    fitter = fitting.LinearLSQFitter()
    model = polynomial.Legendre1D(order, domain=domain)
    for i in range(n_iter):
        best_fit = fitter(model, x_to_fit, y_to_fit, weights=weights)
        # Iteratively reweight to reject outliers
        # (residual/error) < sigma then weight is the weight_class given weight * 1/error^2
        # (residual/error) > sigma then weight is: 0 if weight_class is TukeyBiweight, and some fall off if HuberT
        weights = weight_class.weights(np.abs(y_to_fit - best_fit(x_to_fit)) / error_to_fit) / error_to_fit ** 2.0
    return best_fit
