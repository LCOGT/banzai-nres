import numpy as np
from statsmodels.robust.norms import HuberT
from scipy.optimize import minimize
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
    #weight_class = TukeyBiweight(c=sigma)
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


def stiff_cost_function(model_params, x_to_fit, y_to_fit, weights, model_class: np.polynomial.polynomial.Polynomial,
                        derivative_stiffness):
    model = model_class(model_params)
    cost = (model(x_to_fit) - y_to_fit)**2/weights**2
    # Note the below might just be same as multiplying the parameters by x**2 etc. i.e.
    # x**2 times the acceleration term.
    for i in range(len(model_params) - 1):
        cost += derivative_stiffness[i] * (model.deriv(i+1)(x_to_fit)) ** 2 / weights ** 2
    return np.sum(cost)


def fit_stiff_polynomial(y, error, mask=None, n_iter=5, x=None, sigma=5, order=3, derivative_stiffness=None):
    if derivative_stiffness is None:
        derivative_stiffness = np.zeros(order)
    if mask is None:
        mask = np.zeros(y.shape, dtype=bool)
    if x is None:
        x = np.arange(len(y), dtype=float)
    y_to_fit = y[np.logical_not(mask)]
    error_to_fit = error[np.logical_not(mask)]
    x_to_fit = x[np.logical_not(mask)]
    # start with constant polynomial
    x0 = np.zeros(order + 1, dtype=float)
    x0[0] = np.median(y_to_fit)
    # start with inverse variance weights
    weights = error_to_fit ** -2.0
    weight_class = HuberT(t=sigma)
    model_class = np.polynomial.polynomial.Polynomial
    for i in range(n_iter):
        best_fit_params = minimize(stiff_cost_function, x0=x0, args=(x_to_fit, y_to_fit, weights, model_class,
                                                                     derivative_stiffness),
                                   method='Powell').x
        best_fit = lambda xxx: model_class(best_fit_params)(xxx)
        # Iteratively reweight to reject outliers
        # (residual/error) < sigma then weight is the weight_class given weight * 1/error^2
        # (residual/error) > sigma then weight is: 0 if weight_class is TukeyBiweight, and some fall off if HuberT
        weights = weight_class.weights(np.abs(y_to_fit - best_fit(x_to_fit)) / error_to_fit) / error_to_fit ** 2.0
    return best_fit
