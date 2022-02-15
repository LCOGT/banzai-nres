import numpy as np
import mock

from banzai_nres.qc.qc_wavelength import get_velocity_precision
from banzai_nres.qc.qc_wavelength import AssessWavelengthSolution
from banzai_nres.tests.test_wavelength import TestWavelengthCalibrate
from banzai import context
from xwavecal.utils.wavelength_utils import find_nearest
from banzai.utils.stats import robust_standard_deviation
from astropy import constants
from astropy import units


class TestAssessWavelengthSolution:
    test_image = TestWavelengthCalibrate().generate_image()
    input_context = context.Context({})
    test_qc_results = {}

    def test_do_stage_does_not_crash(self):
        image = AssessWavelengthSolution(self.input_context).do_stage(self.test_image)
        assert image is not None

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
        assert np.isfinite(self.test_qc_results['RVPRECSN'])

    @mock.patch('banzai_nres.qc.qc_wavelength.qc.save_qc_results')
    def test_do_stage_savesqc_toheader(self, fake_post):
        # for now we just test that one of the results (the most important one) was saved to the header.
        image = AssessWavelengthSolution(self.input_context).do_stage(self.test_image)
        assert np.isfinite(image.meta['RVPRECSN'][0])  # check the calculated wavelength precision
        assert len(image.meta['RVPRECSN'][1]) > 0  # description string is not empty

    def test_velocity_precision(self):
        np.random.seed(87213483)
        # make a mock line list
        nlines, expected_precision = 1000, 10 * units.m / units.s
        lab_lines = np.linspace(4000, 5000, nlines)
        # We assume the precision will go like sqrt(n)
        scatter_per_line = expected_precision * np.sqrt(nlines)
        # Our precision is in velocity space so delta lambda = lambda * delta v / c
        features = np.random.normal(lab_lines, scale=lab_lines * scatter_per_line / constants.c)
        velocity_precision = get_velocity_precision(features.value, lab_lines, nlines)
        # For a 1000 lines we expect sigma to be ~3% so we double that here in the tests.
        assert np.isclose(velocity_precision, expected_precision, rtol=6.e-2)

    def test_line_matching(self):
        nlines, wavelength_scatter = 100, 0.1
        mock_lines = np.linspace(4000, 5000, nlines)
        line_list = np.random.permutation(mock_lines)
        features = mock_lines + np.random.randn(nlines) * wavelength_scatter
        lab_lines = find_nearest(features, np.sort(line_list))
        sigma_delta_lambda = robust_standard_deviation(features - lab_lines)
        assert np.isclose(sigma_delta_lambda, wavelength_scatter, rtol=5.e-2)
