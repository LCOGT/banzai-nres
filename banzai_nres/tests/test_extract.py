import numpy as np

from banzai_nres.utils.trace_utils import get_trace_region
from banzai_nres.frames import NRESObservationFrame, EchelleSpectralCCDData
from banzai_nres.extract import WeightedExtract, GetOptimalExtractionWeights
from banzai import context
import pytest


def test_get_region():
    traces = np.zeros((4, 4))
    traces[[2, 3], :] = 1
    trace_yx_pos = get_trace_region(traces == 1)
    assert np.allclose(traces[trace_yx_pos], 1)


class TestExtract:
    def test_rejects_on_no_weights(self):
        con = context.Context({})
        assert GetOptimalExtractionWeights(con).do_stage(two_order_image()) is None

    @pytest.mark.integration
    def test_unit_weights_extraction(self):
        image = two_order_image()
        image.weights = np.ones_like(image.data)
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


    @pytest.mark.integration
    def test_weights_in_poisson_regime(self):
        trace_width, number_traces = 20, 10
        image = five_hundred_square_image(50000, number_traces, trace_width)
        image2 = NRESObservationFrame([EchelleSpectralCCDData(data=image.data, uncertainty=image.uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
        image2.traces = image.traces
        image2.profile = np.ones_like(image.data)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        image = stage.do_stage(image)
        image2.weights = np.ones_like(image2.data)
        stage2 = WeightedExtract(input_context)
        optimal_image = stage2.do_stage(image)
        box_image = stage2.do_stage(image2)
        assert np.allclose(optimal_image.spectrum['flux'], box_image.spectrum['flux'], rtol=0.05)
        assert np.allclose(optimal_image.spectrum['uncertainty'], box_image.spectrum['uncertainty'], rtol=0.05)

    @pytest.mark.integration
    def test_weights_in_readnoise_regime(self):
        trace_width, number_traces = 20, 10
        image = five_hundred_square_image(100, number_traces, trace_width, read_noise=100)
        image2 = NRESObservationFrame([EchelleSpectralCCDData(data=image.data, uncertainty=image.uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
        image2.traces = image.traces
        image2.profile = np.ones_like(image.data)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        image = stage.do_stage(image)
        image2.weights = np.ones_like(image2.data)
        stage2 = WeightedExtract(input_context)
        optimal_image = stage2.do_stage(image)
        box_image = stage2.do_stage(image2)
        optimal_median_sn = np.median(optimal_image.spectrum['flux'] / optimal_image.spectrum['uncertainty'])
        box_median_sn = np.median(box_image.spectrum['flux'] / box_image.spectrum['uncertainty'])
        assert optimal_median_sn > 1.45 * box_median_sn


class TestGetWeights:
    def test_rejects_on_no_profile(self):
        con = context.Context({})
        assert GetOptimalExtractionWeights(con).do_stage(two_order_image()) is None

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

def five_hundred_square_image(maxflux,number_traces,trace_width, read_noise=10):
    traces = np.zeros((500,500))
    data = np.ones_like(traces, dtype=float)
    profile = np.zeros_like(traces, dtype=float)
    ix = np.arange(0,trace_width)
    for i in range (0,number_traces): 
        traces[50*i:50*i+trace_width,:]=i
        for j in range (0,trace_width):
            data[50*i+j,:]+=maxflux*np.exp((-1.)*(ix[j]-trace_width/2.)**2/(trace_width/6.)**2)
        for j in range (0,trace_width):
            profile[50*i+j,:]=data[50*i+j,:]/np.sum(data[50*i:50*i+trace_width,0])

    data += np.random.poisson(data)
    data += np.random.normal(0.0, read_noise, size=data.shape)
    uncertainty = np.sqrt(data + read_noise ** 2)

    image = NRESObservationFrame([EchelleSpectralCCDData(data=data, uncertainty=uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    image.traces, image.profile = traces, profile
    return image
