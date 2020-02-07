import numpy as np
from banzai_nres.traces import find_y_center, TraceInitializer, TraceRefiner
from banzai_nres.images import NRESObservationFrame
from banzai_nres.images import EchelleSpectralCCDData
from banzai import context


def gaussian(y, mu, sigma, a=1.0):
    return a * np.exp(-(y - mu)**2 / 2.0 / sigma**2)


def test_centroid_int_weights():
    nx, ny = 103, 107
    fake_data = np.zeros((ny, nx))
    fake_data[49] = 1
    fake_data[50] = 1
    fake_data[51] = 1
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    center_region = np.logical_and(Y > 10, Y < (ny - 10))
    center, center_errors = find_y_center(Y, center_region, fake_data)
    np.testing.assert_allclose(center, 50.0, atol=0.1)


def test_centroid_flux_weights():
    nx, ny = 103, 107
    fake_data = np.zeros((ny, nx))
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    center_region = np.logical_and(Y > 10, Y < (ny - 10))
    fake_data[center_region] += gaussian(Y[center_region], 45., 3.0, a=10.0)
    center, center_errors = find_y_center(Y, center_region, fake_data)
    np.testing.assert_allclose(center, 45.0, atol=0.1)


def test_blind_solve():
    nx, ny = 401, 403
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_centers = []
    y_0s = [100, 200, 300]
    for i in range(3):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y_0s[i])
        test_data += gaussian(y2d, trace_centers[i], 3)

    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=1e-5*np.ones((ny, nx)),
                                       meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    input_context = context.Context({'trace_separation': 10, 'signal_to_noise_tracing_cutoff': 10,
                                     'min_trace_half_width': 50, 'trace_half_width': 6})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])


def test_refining_on_noisy_data():
    nx, ny = 401, 403
    read_noise = 10.0
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_centers = []
    y_0s = [100, 200, 300]
    for i in range(3):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y_0s[i])
        test_data += gaussian(y2d, trace_centers[i], 3, a=10000.0)

    test_data += np.random.poisson(test_data)
    test_data += np.random.normal(0.0, read_noise, size=test_data.shape)
    uncertainty = np.sqrt(test_data + read_noise ** 2.0)
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=uncertainty,
                                                              meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    input_context = context.Context({'trace_separation': 10, 'signal_to_noise_tracing_cutoff': 10,
                                     'min_trace_half_width': 50, 'trace_half_width': 6})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])


def test_blind_solve_realistic_data():
    nx, ny = 401, 403
    read_noise = 10.0
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_centers = []
    y_0s = [100, 200, 300]
    blaze_function = 1 - 1e-5 * (x2d - nx / 2.) ** 2
    for i in range(3):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y_0s[i])
        test_data += gaussian(y2d, trace_centers[i], 3, a=10000.0) * blaze_function

    # Filter out the numerical noise
    test_data[test_data < 1e-15] = 0.0
    test_data += np.random.poisson(test_data)
    test_data += np.random.normal(0.0, read_noise, size=test_data.shape)
    uncertainty = np.sqrt(test_data + read_noise ** 2.0)
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=uncertainty,
                                                              meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    input_context = context.Context({'trace_separation': 10, 'signal_to_noise_tracing_cutoff': 10,
                                     'min_trace_half_width': 50, 'trace_half_width': 6})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])


def test_refine_traces_with_previous_trace():
    nx, ny = 401, 403
    read_noise = 10.0
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_half_width = 6
    trace_centers = []
    y_0s = [100, 200, 300]
    blaze_function = 1 - 1e-5 * (x2d - nx / 2.) ** 2
    input_traces = np.zeros((ny, nx), dtype=np.int)
    for i in range(3):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y_0s[i])
        test_data += gaussian(y2d, trace_centers[i], 2, a=10000.0) * blaze_function
        input_traces[np.abs(y2d - trace_centers[i]) <= trace_half_width] = i + 1

    # Filter out the numerical noise
    test_data[test_data < 1e-15] = 0.0
    test_data += np.random.poisson(test_data)
    test_data += np.random.normal(0.0, read_noise, size=test_data.shape)
    uncertainty = np.sqrt(test_data + read_noise ** 2.0)
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=uncertainty,
                                                              meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = input_traces
    input_context = context.Context({'trace_separation': 10, 'signal_to_noise_tracing_cutoff': 10,
                                     'min_trace_half_width': 50, 'trace_half_width': trace_half_width})

    stage = TraceRefiner(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])


def test_refine_traces_offset_centroid():
    nx, ny = 401, 403
    read_noise = 10.0
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_half_width = 6
    trace_centers = []
    y_0s = [100, 200, 300]
    blaze_function = 1 - 1e-5 * (x2d - nx / 2.) ** 2
    input_traces = np.zeros((ny, nx), dtype=np.int)
    for i in range(3):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y_0s[i])
        test_data += gaussian(y2d, trace_centers[i], 2, a=10000.0) * blaze_function
        input_traces[np.abs(y2d - trace_centers[i]) <= trace_half_width] = i + 1

    # Filter out the numerical noise
    test_data[test_data < 1e-15] = 0.0
    test_data += np.random.poisson(test_data)
    test_data += np.random.normal(0.0, read_noise, size=test_data.shape)
    uncertainty = np.sqrt(test_data + read_noise ** 2.0)
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=uncertainty,
                                                              meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = input_traces
    input_context = context.Context({'trace_separation': 10, 'signal_to_noise_tracing_cutoff': 10,
                                     'min_trace_half_width': 50, 'trace_half_width': trace_half_width})

    stage = TraceRefiner(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center + 1) <= 4])
