import numpy as np
from banzai.stages import Stage
from banzai.utils import qc
from astropy import constants
from astropy import units
from xwavecal.utils.wavelength_utils import find_nearest
from banzai.utils.stats import robust_standard_deviation
import logging


logger = logging.getLogger('banzai')


class AssessWavelengthSolution(Stage):

    def __init__(self, runtime_context):
        super(AssessWavelengthSolution, self).__init__(runtime_context)

    def do_stage(self, image):
        lab_lines = find_nearest(image.features['wavelength'], np.sort(image.line_list))
        delta_lambda = image.features['wavelength'] - lab_lines

        sigma_delta_lambda = robust_standard_deviation(delta_lambda)
        low_scatter_lines = np.abs(delta_lambda) < 3. * sigma_delta_lambda

        matched_sigma_delta_lambda = robust_standard_deviation(delta_lambda[low_scatter_lines])
        num_detected_lines = len(image.features['wavelength'])
        num_matched_lines = np.count_nonzero(low_scatter_lines)

        feature_centroid_uncertainty = image.features['centroid_err']

        reduced_chi2 = get_reduced_chi_squared(delta_lambda[low_scatter_lines],
                                               feature_centroid_uncertainty[low_scatter_lines])
        velocity_precision = get_velocity_precision(image.features['wavelength'][low_scatter_lines],
                                                    lab_lines[low_scatter_lines], num_matched_lines)

        if num_matched_lines == 0:  # get rid of nans in the matched statistics if we have zero matched lines.
            matched_sigma_delta_lambda, reduced_chi2, velocity_precision = 0, 0, 0 * units.meter/units.second

        # opensearch keys don't have to be the same as the fits headers
        qc_results = {'SIGLAM': np.round(matched_sigma_delta_lambda, 4),
                      'RVPRECSN': np.round(velocity_precision.to(units.meter/units.second).value, 4),
                      'WAVRCHI2': np.round(reduced_chi2, 4),
                      'NLINEDET': num_detected_lines,
                      'NLINES': num_matched_lines}
        qc_description = {'SIGLAM': 'wavecal residuals [Angstroms]',
                          'RVPRECSN': 'wavecal precision [m/s]',
                          'WAVRCHI2': 'reduced chisquared goodness of wavecal fit',
                          'NLINEDET': 'Number of lines found on detector',
                          'NLINES': 'Number of matched lines'}
        qc.save_qc_results(self.runtime_context, qc_results, image)
        # saving the results to the image header
        for key in qc_results.keys():
            image.meta[key] = (qc_results[key], qc_description[key])

        logger.info(f'wavecal precision (m/s) = {qc_results["RVPRECSN"]}', image=image)
        if qc_results['RVPRECSN'] > 10 or qc_results['RVPRECSN'] < 3:
            logger.warning(f' Final calibration precision is outside the expected range '
                           f'wavecal precision (m/s) = '
                           f'{qc_results["RVPRECSN"]}', image=image)
        return image


def get_reduced_chi_squared(values, uncertainty):
    return np.sum((values / uncertainty)**2) / len(values)


def get_velocity_precision(image_lines, lab_lines, num_matched_lines):
    """
    calculating metrics in velocity space (easily understood by users) del lambda/ lambda * c = delta v.
    then divide delta v by square root of the number of lines, giving the error on the mean of the residuals.
    """
    delta_lambda = image_lines - lab_lines
    dlam_overlam = delta_lambda / lab_lines
    velocity_precision = robust_standard_deviation(dlam_overlam) / np.sqrt(num_matched_lines) * constants.c
    return velocity_precision
