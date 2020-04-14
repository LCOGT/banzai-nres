import numpy as np

import logging

from banzai.stages import Stage
from banzai.calibrations import CalibrationStacker, CalibrationUser

from banzai_nres.wavelength import LineListLoader

class AssessWavelengthSolution(Stage):
    """
    Calculate the dispersion in the wavelength solution and assess whether it is photon-limited or not
    @author:mjohnson
    """

    def __init__(self,runtime_context):
        super(AssessWavelengthSolution, self).__init__(runtime_context)
    
    def do_stage(self,image):
        line_list = LineListLoader.do_stage(image)
        self.calculate_dispersion(image,line_list,raw_dispersion,good_dispersion,raw_chi_squared,good_chi_squared)
        qc_results = {'raw_dispersion':raw_dispersion, 'good_dispersion':good_dispersion, 'raw_chi_squared':raw_chi_squared, 'good_chi_squared':good_chi_squared}
        qc.save_qc_results(self.runtime_context, qc_results, image)

    def calculate_dispersion(self,image,line_list):
        measured_wavelengths = image.features['wavelength']
        #actual_wavelengths = image.wcs.measured_lines['wavelength']
        difference = measured_wavelengths - line_list
        raw_dispersion = np.std(difference)
        low_scatter_lines = np.isclose(difference,0,atol=0.1)
        good_dispersion = np.std(difference[low_scatter_lines])
        self.calculate_centroid_errors(image,feature_centroid_uncertainty)
        raw_chi_squared = np.sum((difference/feature_centroid_uncertainty)**2)/len(difference)
        good_chi_squared = np.sum((difference[low_scatter_lines]/feature_centroid_uncertainty[low_scatter_lines])**2)/len(difference[low_scatter_lines])
        return raw_dispersion, good_dispersion, raw_chi_squared, good_chi_squared

    def calculate_centroid_errors(self,image):
        feature_width = 2.0 #pixels, guess for now -- assuming all features have same width
        feature_flux = image.features['flux']*np.sqrt(2.*np.pi)/feature_width
        feature_centroid_uncertainty = 0.412 * feature_flux * 2. * feature_width
        return feature_centroid_uncertainty