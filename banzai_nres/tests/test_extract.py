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

    @pytest.mark.integration
    def test_weights_in_poisson_regime(self):
        trace_width=20
        image = five_hundred_square_image(1000,10,trace_width)
        image2 = NRESObservationFrame([EchelleSpectralCCDData(data=image.data, uncertainty=image.uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
        image2.traces = image.traces
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        image = stage.do_stage(image)
        image2.weights = np.ones_like(image2.data)/trace_width #correct for normalization
        stage2 = WeightedExtract(input_context)
        optimal_image = stage2.do_stage(image)
        box_image = stage2.do_stage(image2)
        assert np.allclose(optimal_image.spectrum['uncertainty'], box_image.spectrum['uncertainty'],0.1)

    #note the following test does not yet work as intended
    @pytest.mark.integration
    def test_weights_in_readnoise_regime(self):
        trace_width=20
        image = five_hundred_square_image(1,10,trace_width)
        image2 = NRESObservationFrame([EchelleSpectralCCDData(data=image.data, uncertainty=image.uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
        image2.traces = image.traces
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        image = stage.do_stage(image)
        image2.weights = np.ones_like(image2.data)/trace_width #correct for normalization
        stage2 = WeightedExtract(input_context)
        optimal_image = stage2.do_stage(image)
        box_image = stage2.do_stage(image2)
        import pdb; pdb.set_trace()
        assert np.allclose(optimal_image.spectrum['uncertainty']**2, box_image.spectrum['uncertainty']**2,0.1)

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

    def test_optimal_weights_on_one_profile(self):
        image = two_order_image()
        image.profile = np.ones_like(image.traces, dtype=float)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        output_image = stage.do_stage(image)
        assert np.allclose(output_image.weights[~np.isclose(image.traces, 0)], 1/3)

    def test_optimal_weights_on_box_profile(self):
        profile = np.zeros((5,5))
        variance = np.ones_like(profile, dtype=float)
        mask = np.zeros_like(profile, dtype=float)
        profile[1:4,:], variance[1:4,:] = 1., 1.
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        weights = stage.weights(profile,variance,mask)
        #check that the weights in the order are 1/width of the order:
        assert np.allclose(weights[np.isclose(profile,1)],1./3.)
        #check that the weights of the area with zero profile are zero:
        assert np.allclose(weights[np.isclose(profile,0)],0)
        


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

def five_hundred_square_image(maxflux,number_traces,trace_width):
    traces = np.zeros((500,500))
    for i in range (0,number_traces): traces[50*i:50*i+trace_width,:]=i
    data = np.ones_like(traces, dtype=float)
    data[~np.isclose(traces, 0)] = maxflux
    data+=np.random.randn(500,500)+10.
    data+=np.random.poisson(data)
    uncertainty = 10. + np.sqrt(data)
    image = NRESObservationFrame([EchelleSpectralCCDData(data=data, uncertainty=uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    image.traces = traces
    return image
