import numpy as np

from banzai_nres.utils.extract_utils import get_region
from banzai_nres.frames import NRESObservationFrame, EchelleSpectralCCDData
from banzai_nres.extract import WeightedExtract, GetOptimalExtractionWeights
from banzai import context
import pytest


def test_get_region():
    traces = np.zeros((4, 4))
    traces[[2, 3], :] = 1
    trace_yx_pos = get_region(traces == 1)
    assert np.allclose(traces[trace_yx_pos], 1)


class TestExtract:
    @pytest.mark.integration
    def test_unit_weights_extraction(self):
        image = two_order_image()
        expected_extracted_flux = np.max(image.data) * 3
        expected_extracted_uncertainty = np.sqrt(3) * np.max(image.data)
        input_context = context.Context({})

        stage = WeightedExtract(input_context)
        output_image = stage.do_stage(image)
        spectrum = output_image.spectrum

        assert np.allclose(spectrum['flux'][0], expected_extracted_flux)
        assert np.allclose(spectrum['uncertainty'][0], expected_extracted_uncertainty)
        assert np.allclose(spectrum['id'][0], 1)

        assert np.allclose(spectrum['flux'][1][1:-1], expected_extracted_flux)
        assert np.allclose(spectrum['flux'][1][[0, -1]], 0)
        assert np.allclose(spectrum['uncertainty'][1:-1], expected_extracted_uncertainty)
        assert np.allclose(spectrum['uncertainty'][1][[0, -1]], 0)
        assert np.allclose(spectrum['id'][1], 2)

        assert np.allclose(spectrum['pixel'], np.arange(image.data.shape[1]))

    # TODO write a test which tests that: the optimal extraction weights are correct, and tests that test whether:
    #  1. the variance of the optimal extracted spectrum is the same as the box extracted spectrum in the Poisson noise regime.
    #  2. the variance of the optimal extracted spectrum is ~1.69 times the box extracted spectrum in the read noise regime


class TestGetWeights:
    def test_optimal_weights_zero_on_zero_profile(self):
        image = two_order_image()
        image.profile = np.zeros_like(image.traces, dtype=float)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        output_image = stage.do_stage(image)
        assert np.allclose(output_image.weights, 0)

    def test_optimal_weights_on_box_profile(self):
        image = two_order_image()
        image.profile = np.ones_like(image.traces, dtype=float)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        output_image = stage.do_stage(image)
        assert np.allclose(output_image.weights[~np.isclose(image.traces, 0)], 1/3)


def two_order_image():
    # generate 2 flat traces.
    traces = np.zeros((60, 20))
    traces[[10, 11, 12], :] = 1
    # the second trace that does not span the image entirely
    traces[[50, 51, 53], 1:-1] = 2
    # generate test data with zero noise
    data = np.ones_like(traces, dtype=float)
    data[~np.isclose(traces, 0)] = 100.
    uncertainty = 1. * data
    image = NRESObservationFrame([EchelleSpectralCCDData(data=data, uncertainty=uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    image.traces = traces
    return image