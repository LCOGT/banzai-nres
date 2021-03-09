import numpy as np
import mock
import pytest

from banzai_nres.qc.qc_wavelength import AssessWavelengthSolution
from banzai_nres.tests.test_wavelength import TestWavelengthCalibrate
from banzai import context


class TestAssessWavelengthSolution:
    test_image = TestWavelengthCalibrate().generate_image()
    input_context = context.Context({})
    test_qc_results = {}

    def test_do_stage_does_not_crash(self):
        image = AssessWavelengthSolution(self.input_context).do_stage(self.test_image)
        assert image is not None


    @pytest.mark.integration
    @mock.patch('banzai_nres.qc.qc_wavelength.qc.save_qc_results')
    def test_do_stage_posts_to_elastic_search(self, fake_post):
        # define a test function to mock the saving of the quality control metrics

        def save_qc_results_locally(runtime_context, qc_results, image):
            self.test_qc_results = qc_results
            return None
        fake_post.side_effect = save_qc_results_locally
        # run the stage
        AssessWavelengthSolution(self.input_context).do_stage(self.test_image)
        assert len(self.test_qc_results) == 5
        assert np.isfinite(self.test_qc_results['PRECISN'])

    @mock.patch('banzai_nres.qc.qc_wavelength.qc.save_qc_results')
    def test_do_stage_savesqc_toheader(self, fake_post):
        # for now we just test that one of the results (the most important one) was saved to the header.
        image = AssessWavelengthSolution(self.input_context).do_stage(self.test_image)
        assert np.isfinite(image.meta['PRECISN'][0])  # check the calculated wavelength precision
        assert len(image.meta['PRECISN'][1]) > 0  # description string is not empty

    def test_quality_metrics(self):
        lab_lines = self.test_image.features['wavelength']
        Delta_lambda = AssessWavelengthSolution(self.input_context).calculate_delta_lambda(self.test_image, lab_lines)
        result = AssessWavelengthSolution(self.input_context).calculate_1d_metrics(self.test_image, Delta_lambda,
                                                                                   lab_lines)
        sigma_Dlambda, good_sigma_Dlambda, raw_chi_squared, good_chi_squared, num_matched_lines, _ = result
        assert sigma_Dlambda >= good_sigma_Dlambda
        assert raw_chi_squared >= good_chi_squared
