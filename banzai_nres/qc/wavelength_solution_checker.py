import numpy as np

from banzai.stages import Stage
from banzai.calibrations import CalibrationStacker, CalibrationUser
from banzai_nres.wavelength import LineListLoader
from xwavecal.utils.wavelength_utils import find_nearest

import logging

class AssessWavelengthSolution(Stage):
    """
    Calculate the dispersion in the wavelength solution and assess whether it is photon-limited or not
    @author:mjohnson
    """

    def __init__(self, runtime_context):
        super(AssessWavelengthSolution, self).__init__(runtime_context)
    
    def do_stage(self,image):
        line_list = image.line_list
        lab_lines = find_nearest(image.features['wavelength'], np.sort(line_list))
        self.calculate_delta_lambda(image,lab_lines,Delta_lambda)
        self.calculate_1d_metrics(image,Delta_lambda,sigma_Dlambda,good_sigma_Dlambda,raw_chi_squared,good_chi_squared,)
        self.calculate_2d_metrics(image,Delta_lambda,x_diff_Dlambda, order_diff_Dlambda)
        qc_results = {'sigma_Dlambda':sigma_Dlambda, 'good_sigma_Dlambda':good_sigma_Dlambda, 
                        'raw_chi_squared':raw_chi_squared, 'good_chi_squared':good_chi_squared, 
                        'Delta_lambda':Delta_lambda, 'x_diff_Dlambda':x_diff_Dlambda, 
                        'order_diff_Dlambda':order_diff_Dlambda}
        qc.save_qc_results(self.runtime_context, qc_results, image)

    def calculate_delta_lambda(self,image,lab_lines):
        measured_wavelengths = image.features['wavelength']
        Delta_lambda = measured_wavelengths - lab_lines
        return Delta_lambda

    def calculate_1d_metrics(self,image,Delta_lambda):
        sigma_Dlambda = np.std(Delta_lambda)
        low_scatter_lines = np.isclose(Delta_lambda,0,atol=0.1)
        good_sigma_Dlambda = np.std(Delta_lambda[low_scatter_lines])
        feature_centroid_uncertainty = image.features['centroid_err']
        raw_chi_squared = np.sum((Delta_lambda/feature_centroid_uncertainty)**2)/len(Delta_lambda)
        good_chi_squared = np.sum((Delta_lambda[low_scatter_lines]/feature_centroid_uncertainty[low_scatter_lines])**2)/len(Delta_lambda[low_scatter_lines])
        return sigma_Dlambda, good_sigma_Dlambda, raw_chi_squared, good_chi_squared

    def calculate_2d_metrics(self,image,Delta_lambda):
        x, order = image.features['pixel'], image.features['order']
        bins = 20
        histogramx, bins_x = np.histogram(x, bins=bins, range=([np.min(x)-1,np.max(x)+1]))
        histogramo, bins_order = np.histogram(order, bins=bins, range=([np.min(order)-1,np.max(order)+1]))
        x_indices, order_indices = np.digitize(x, bins_x[1:]), np.digitize(order, bins_order[1:])
        x_diff_Dlambda, order_diff_Dlambda = np.zeros_like(histogramx), np.zeros_like(histogramo)
        for i in range (0,bins):
            if histogramx[i] != 0: x_diff_Dlambda[i] = np.mean(Delta_lambda[x_indices == i])
            if histogramo[i] != 0: order_diff_Dlambda[i] = np.mean(Delta_lambda[order_indices == i])
        return x_diff_Dlambda, order_diff_Dlambda


        #TO DO:
        #clean up code