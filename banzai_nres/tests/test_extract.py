import numpy as np

from banzai.data import CCDData
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.extract import WeightedExtract, GetOptimalExtractionWeights
from banzai import context


class TestExtract:
    def test_rejects_on_no_weights(self):
        con = context.Context({})
        assert WeightedExtract(con).do_stage(two_order_image()) is None

    def test_rejects_on_no_wavelengths(self):
        con = context.Context({})
        image = type('image', (), {'weights': 'notnone', 'wavelengths': None})
        assert WeightedExtract(con).do_stage(image) is None

    def test_unit_weights_extraction(self):
        image = two_order_image()
        image.weights = np.ones_like(image.data)
        expected_extracted_flux = np.max(image.data) * 3
        expected_extracted_wavelength = np.max(image.wavelengths)
        expected_extracted_uncertainty = np.sqrt(3) * np.max(image.data)
        input_context = context.Context({})

        stage = WeightedExtract(input_context)
        output_image = stage.do_stage(image)
        spectrum = output_image.spectrum

        assert np.allclose(spectrum[0, 0]['flux'], expected_extracted_flux)
        assert np.allclose(spectrum[0, 0]['wavelength'], expected_extracted_wavelength)
        assert np.allclose(spectrum[0, 0]['uncertainty'], expected_extracted_uncertainty)
        assert np.allclose(spectrum[0, 0]['id'], 1)

        assert np.allclose(spectrum[1, 1]['flux'][1:-1], expected_extracted_flux)
        assert np.allclose(spectrum[1, 1]['wavelength'][1:-1], expected_extracted_wavelength)
        assert len(spectrum[1, 1]['flux']) == image.traces.shape[1] - 2
        assert np.allclose(spectrum[1, 1]['uncertainty'][1:-1], expected_extracted_uncertainty)
        assert len(spectrum[1, 1]['uncertainty']) == image.traces.shape[1] - 2
        assert np.allclose(spectrum[1, 1]['id'], 2)

    def test_extract_in_poisson_regime(self):
        trace_width, number_traces = 20, 10
        seed = 1408235915
        image = five_hundred_square_image(50000, number_traces, trace_width, seed=seed)
        expected_extracted_wavelength = np.max(image.wavelengths)

        image2 = five_hundred_square_image(50000, number_traces, trace_width, seed=seed)
        image2.profile = np.ones_like(image.data)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        image = stage.do_stage(image)
        image2.weights = np.ones_like(image2.data)
        stage2 = WeightedExtract(input_context)
        optimal_image = stage2.do_stage(image)
        box_image = stage2.do_stage(image2)
        for i in range(1, number_traces + 1):
            assert np.allclose(optimal_image.spectrum[i, i]['flux'], box_image.spectrum[i, i]['flux'], rtol=0.05)
            assert np.allclose(optimal_image.spectrum[i, i]['wavelength'], expected_extracted_wavelength)
            assert np.allclose(optimal_image.spectrum[i, i]['uncertainty'], box_image.spectrum[i, i]['uncertainty'],
                               rtol=0.05)

    def test_extract_in_readnoise_regime(self):
        trace_width, number_traces = 20, 10
        seed = 192074123
        image = five_hundred_square_image(100, number_traces, trace_width, read_noise=100, seed=seed)
        image2 = five_hundred_square_image(100, number_traces, trace_width, read_noise=100, seed=seed)
        image2.profile = np.ones_like(image.data)
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        image = stage.do_stage(image)
        image2.weights = np.ones_like(image2.data)
        stage2 = WeightedExtract(input_context)
        optimal_image = stage2.do_stage(image)
        box_image = stage2.do_stage(image2)
        for i in range(1, number_traces + 1):
            optimal_median_sn = np.median(
                optimal_image.spectrum[i, i]['flux'] / optimal_image.spectrum[i, i]['uncertainty'])
            box_median_sn = np.median(box_image.spectrum[i, i]['flux'] / box_image.spectrum[i, i]['uncertainty'])
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
        assert np.allclose(output_image.weights[~np.isclose(image.traces, 0)], 1 / 3)

    def test_optimal_weights_on_box_profile(self):
        profile = np.zeros((5, 5))
        variance = np.ones_like(profile, dtype=float)
        mask = np.zeros_like(profile, dtype=float)
        profile[1:4, :], variance[1:4, :] = 1., 1.
        input_context = context.Context({})

        stage = GetOptimalExtractionWeights(input_context)
        weights = stage.weights(profile, variance, mask)
        # check that the weights in the order are 1/width of the order:
        assert np.allclose(weights[np.isclose(profile, 1)], 1. / 3.)
        # check that the weights of the area with zero profile are zero:
        assert np.allclose(weights[np.isclose(profile, 0)], 0)


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
    wavelengths = np.ones_like(traces) * 5  # dummy wavelengths image that has values distinct from flux and traces.
    image = NRESObservationFrame([CCDData(data=data, uncertainty=uncertainty, meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    image.wavelengths = wavelengths
    image.traces = traces
    image.fibers = {'fiber': np.arange(2), 'order': np.arange(2)}
    image.blaze = {'id': np.arange(2), 'blaze': [np.arange(20), np.arange(20)],
                   'blaze_error': [np.arange(20), np.arange(20)]}
    return image


def five_hundred_square_image(maxflux, number_traces, trace_width, read_noise=10, seed=None):
    traces = np.zeros((500, 500))
    data = np.ones_like(traces, dtype=float)
    profile = np.zeros_like(traces, dtype=float)
    ix = np.arange(trace_width)
    for i in range(1, number_traces + 1):
        traces[40 * i:40 * i + trace_width, :] = i
        for j in range(0, trace_width):
            data[40 * i + j, :] += maxflux * np.exp((-1.) * (ix[j] - trace_width / 2.) ** 2 / (trace_width / 6.) ** 2)
        for j in range(0, trace_width):
            profile[40 * i + j, :] = data[40 * i + j, :] / np.sum(data[40 * i: 40 * i + trace_width, 0])
    np.random.seed(seed=seed)
    data += np.random.poisson(data)
    data += np.random.normal(0.0, read_noise, size=data.shape)
    uncertainty = np.sqrt(data + read_noise ** 2)
    wavelengths = np.ones_like(traces) * 5  # dummy wavelengths image that has values distinct from flux and traces.
    image = NRESObservationFrame([CCDData(data=data, uncertainty=uncertainty, meta={'OBJECTS': 'tung&tung&none'})],
                                 'foo.fits')
    image.traces = traces,
    image.profile = profile
    image.wavelengths = wavelengths
    image.blaze = {'id': np.arange(number_traces) + 1,
                   'blaze': [np.ones(traces.shape[1]) for i in range(number_traces)],
                   'blaze_error': [np.ones(traces.shape[1]) for i in range(number_traces)]}
    image.fibers = {'fiber': np.arange(number_traces) + 1, 'order': np.arange(number_traces) + 1}
    return image
