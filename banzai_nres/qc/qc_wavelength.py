import numpy as np

from banzai.stages import Stage
from banzai.utils import qc
from astropy import constants as const
from astropy import units
from xwavecal.utils.wavelength_utils import find_nearest
import logging

from scipy.stats import binned_statistic


logger = logging.getLogger('banzai')


class AssessWavelengthSolution(Stage):
    """
    Calculate the dispersion in the wavelength solution
    @author:mjohnson
    """

    def __init__(self, runtime_context):
        super(AssessWavelengthSolution, self).__init__(runtime_context)

    def do_stage(self, image):
        Dlambda_match_threshold = 0.1  # Angstroms
        lab_lines = find_nearest(image.features['wavelength'], np.sort(image.line_list))
        Dlambda = self.calculate_delta_lambda(image, lab_lines)
        result = self.calculate_1d_metrics(image, Dlambda, lab_lines, Dlambda_match_threshold)
        sigma_Dlambda, matched_sigma_Dlambda, chi2, matched_chi2, num_matched_lines, velocity_precision = result
        #x_diff_Dlambda, order_diff_Dlambda = self.calculate_2d_metrics(image, Dlambda)
        # TODO add the number of overlaps matched during wavelength calibration to this metric.
        #  the number of fit overlaps may be saved in the features header, check this.
        qc_results = {'sigma_Dlambda': sigma_Dlambda,
                      'matched_sigma_Dlambda': matched_sigma_Dlambda,
                      'wavecal_precision(m/s)': velocity_precision.to(units.meter/units.second).value,
                      'chi_squared': chi2,
                      'matched_chi_squared': matched_chi2,
                      'num_matched_lines': num_matched_lines}
        qc_description = {'sigma_Dlambda': 'Wavecal statistic. Standard deviation of the distribution of '
                                           'lab minus measured wavelength residuals',
                          'matched_sigma_Dlambda': 'Wavecal statistic. sigma_Dlambda but only using lines that agree with a '
                                                   f'lab line to within {Dlambda_match_threshold} Angstroms.',
                          'wavecal_precision(m/s)': 'matched_sigma_Dlambda/sqrt(N) in velocity space '
                                                    '(i.e. standard deviation of delta lambda/lambda * c divided by '
                                                    'the square root of the number of matched lines) '
                                                    'for matched lines only. I.e. this is the formal error on the mean '
                                                    'of the residuals distribution (in velocity units). '
                                                    'Units of meter per second.',
                          'chi_squared': 'Wavecal statistic. Formal chi-squared statistic of the wavelength residuals '
                                         'using their formal errors.',
                          'matched_chi_squared': 'Wavecal statistic. chi_squared but only using lines that agree with a '
                                                 f'lab line to within {Dlambda_match_threshold} Angstroms.',
                          'num_matched_lines': f'Wavecal statistic. Number of measured lines that agree with a '
                          f'lab line to within {Dlambda_match_threshold} Angstroms.'}
        qc.save_qc_results(self.runtime_context, qc_results, image)
        # saving the results to the image header
        for key in qc_results.keys():
            image.meta[key] = (qc_results[key], qc_description[key])

        # print the most easily understood metric to the log
        logger.info(f'wavecal precision (m/s) = {qc_results["wavecal_precision(m/s)"]}', image=image)
        if qc_results['wavecal_precision(m/s)'] > 15 or qc_results['wavecal_precision(m/s)'] < 3:
            logger.warning(f' Final calibration precision is too large or small. '
                           f'wavecal precision (m/s) = '
                           f'{qc_results["wavecal_precision(m/s)"]}', image=image)
        return image

    def calculate_delta_lambda(self, image, lab_lines):
        """
        :param image:
        :param lab_lines: list of laboratory wavelengths
        :return: Dlambda: ndarray. The residual between the measured wavelength and the laboratory wavelength

        """
        measured_wavelengths = image.features['wavelength']
        Dlambda = measured_wavelengths - lab_lines
        return Dlambda

    def calculate_1d_metrics(self, image, Delta_lambda, lab_lines, Dlambda_match_threshold=0.1):
        sigma_Dlambda = np.std(Delta_lambda)
        low_scatter_lines = np.isclose(Delta_lambda, 0, atol=Dlambda_match_threshold)
        num_matched_lines = np.count_nonzero(low_scatter_lines)
        matched_sigma_Dlambda = np.std(Delta_lambda[low_scatter_lines])
        feature_centroid_uncertainty = image.features['centroid_err']
        chi2 = np.sum((Delta_lambda/feature_centroid_uncertainty)**2)/len(Delta_lambda)
        matched_chi2 = np.sum((Delta_lambda[low_scatter_lines]/feature_centroid_uncertainty[low_scatter_lines])**2)/len(Delta_lambda[low_scatter_lines])
        # calculating metrics in velocity space (easily understood by users) del lambda/ lambda * c = delta v.
        # then divide delta v by square root of the number of lines, giving the error on the mean of the residuals.
        velocity_precision = np.std((Delta_lambda / lab_lines * const.c)[low_scatter_lines]) / np.sqrt(num_matched_lines)
        if num_matched_lines == 0:  # get rid of nans in the matched statistics if we have zero matched lines.
            matched_sigma_Dlambda, matched_chi2, velocity_precision = 0, 0, 0 * units.meter/units.second
        return sigma_Dlambda, matched_sigma_Dlambda, chi2, matched_chi2, num_matched_lines, velocity_precision

    def calculate_2d_metrics(self, image, Delta_lambda):
        """
        :param image:
        :param Delta_lambda:
        :return: Dlambda_vs_x, Dlambda_vs_order: ndarray, ndarray

        Dlambda_vs_x: The wavelength residual binned by pixel, into 20 bins.
        Dlambda_vs_order: The wavelength residual binned by order number, into 20 bins (so roughly 3 orders per bin).
        """
        x, order = image.features['pixel'], image.features['order']
        bins = 20
        Dlambda_vs_x = binned_statistic(x, Delta_lambda, statistic='mean',
                                         bins=bins, range=([np.min(x)-1, np.max(x)+1]))
        Dlambda_vs_order = binned_statistic(order, Delta_lambda, statistic='mean',
                                            bins=bins, range=([np.min(order)-1, np.max(order)+1]))
        return Dlambda_vs_x.statistic, Dlambda_vs_order.statistic
