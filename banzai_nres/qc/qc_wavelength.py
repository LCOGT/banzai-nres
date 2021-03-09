import numpy as np

from banzai.stages import Stage
from banzai.utils import qc
from xwavecal.utils.wavelength_utils import find_nearest

from scipy.stats import binned_statistic


class AssessWavelengthSolution(Stage):
    """
    Calculate the dispersion in the wavelength solution and assess whether it is photon-limited or not
    @author:mjohnson
    """

    def __init__(self, runtime_context):
        super(AssessWavelengthSolution, self).__init__(runtime_context)

    def do_stage(self, image):
        lab_lines = find_nearest(image.features['wavelength'], np.sort(image.line_list))
        Delta_lambda = self.calculate_delta_lambda(image, lab_lines)
        sigma_Dlambda, good_sigma_Dlambda, raw_chi_squared, good_chi_squared = self.calculate_1d_metrics(image,
                                                                                                         Delta_lambda)
        x_diff_Dlambda, order_diff_Dlambda = self.calculate_2d_metrics(image, Delta_lambda)
        qc_results = {'sigma_Dlambda': sigma_Dlambda, 'good_sigma_Dlambda': good_sigma_Dlambda,
                      'raw_chi_squared': raw_chi_squared, 'good_chi_squared': good_chi_squared,
                      'Delta_lambda': Delta_lambda, 'x_diff_Dlambda': x_diff_Dlambda,
                      'order_diff_Dlambda': order_diff_Dlambda}
        qc.save_qc_results(self.runtime_context, qc_results, image)
        return image

    def calculate_delta_lambda(self, image, lab_lines):
        measured_wavelengths = image.features['wavelength']
        Delta_lambda = measured_wavelengths - lab_lines
        return Delta_lambda

    def calculate_1d_metrics(self, image, Delta_lambda):
        sigma_Dlambda = np.std(Delta_lambda)
        low_scatter_lines = np.isclose(Delta_lambda, 0, atol=0.1)
        good_sigma_Dlambda = np.std(Delta_lambda[low_scatter_lines])
        feature_centroid_uncertainty = image.features['centroid_err']
        raw_chi_squared = np.sum((Delta_lambda/feature_centroid_uncertainty) ** 2) / len(Delta_lambda)
        good_chi_squared = (Delta_lambda[low_scatter_lines] / feature_centroid_uncertainty[low_scatter_lines]) ** 2
        good_chi_squared /= len(Delta_lambda[low_scatter_lines])
        good_chi_squared = np.sum(good_chi_squared)
        return sigma_Dlambda, good_sigma_Dlambda, raw_chi_squared, good_chi_squared

    def calculate_2d_metrics(self, image, Delta_lambda):
        x, order = image.features['pixel'], image.features['order']
        bins = 20
        x_diff_Dlambda = binned_statistic(x, Delta_lambda, statistic='mean', bins=bins,
                                          range=([np.min(x) - 1, np.max(x) + 1]))
        order_diff_Dlambda = binned_statistic(order, Delta_lambda, statistic='mean', bins=bins,
                                              range=([np.min(order) - 1, np.max(order) + 1]))
        return x_diff_Dlambda.statistic, order_diff_Dlambda.statistic
