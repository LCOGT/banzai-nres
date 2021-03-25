import numpy as np
from banzai_nres.traces import find_y_center, TraceInitializer, TraceRefiner, refine_traces
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.frames import EchelleSpectralCCDData
from banzai import context


def gaussian(y, mu, sigma, a=1.0):
    return a * np.exp(-(y - mu)**2 / 2.0 / sigma**2)


def make_realistic_trace_image(nx=401, ny=403, y0_centers=[100, 200, 300]):
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_centers = []
    for i in range(len(y0_centers)):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y0_centers[i])
        test_data += gaussian(y2d, trace_centers[i], 3)
    return trace_centers, test_data, x2d, y2d


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
    trace_centers, test_data, x2d, y2d = make_realistic_trace_image()
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=1e-5*np.ones_like(test_data),
                                       meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)
    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])


def test_blind_solve_with_bpm():
    trace_centers, test_data, x2d, y2d = make_realistic_trace_image(y0_centers=[100, 200, 300])
    bpm_mask = np.ones_like(test_data, dtype=bool)
    bpm_mask[150:250, :] = 0  # set the pixels around the center trace as good. Leave the other pixels masked.
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=1e-5*np.ones_like(test_data),
                                       meta={'OBJECTS': 'tung&tung&none'}, mask=bpm_mask)], 'foo.fits')
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)
    assert np.count_nonzero(list(set(output_image.traces[:, 200]))) == 1
    # test that the only valid trace is centered correctly
    assert all(output_image.traces[np.abs(y2d - trace_centers[1]) <= 4])


def test_blind_solve_with_edge_clipping_traces():
    trace_centers, test_data, x2d, y2d = make_realistic_trace_image(nx=401, ny=403, y0_centers=[1, 200, 400])
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=1e-5*np.ones_like(test_data),
                                       meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)
    assert np.count_nonzero(list(set(output_image.traces[:, 200]))) == 1
    # test that the only valid trace is centered correctly
    assert all(output_image.traces[np.abs(y2d - trace_centers[1]) <= 4])


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
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})
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
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})
    stage = TraceInitializer(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])


def test_refine_traces_with_previous_trace():
    nx, ny = 401, 403
    num_traces = 3
    read_noise = 10.0
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_half_width = 6
    trace_centers = []
    y_0s = [100, 200, 300]
    blaze_function = 1 - 1e-5 * (x2d - nx / 2.) ** 2
    input_traces = np.zeros((ny, nx), dtype=np.int)
    for i in range(num_traces):
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
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})

    stage = TraceRefiner(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center) <= 4])
    assert np.isclose(num_traces, output_image.num_traces)


def test_refine_traces_offset_centroid():
    nx, ny = 401, 403
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))
    test_data, trace_centers, input_traces = make_simple_traces(nx, ny)
    read_noise = 10.0

    # Filter out the numerical noise
    test_data[test_data < 1e-15] = 0.0
    test_data += np.random.poisson(test_data)
    test_data += np.random.normal(0.0, read_noise, size=test_data.shape)
    uncertainty = np.sqrt(test_data + read_noise ** 2.0)
    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=uncertainty,
                                                              meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = input_traces
    input_context = context.Context({'TRACE_HALF_HEIGHT': 5})

    stage = TraceRefiner(input_context)
    output_image = stage.do_stage(test_image)

    for trace_center in trace_centers:
        # Make sure that the center +- 4 pixels is in the trace image
        assert all(output_image.traces[np.abs(y2d - trace_center + 1) <= 4])


def make_simple_traces(nx=401, ny=403, trace_half_width=6, blaze=True):
    test_data = np.zeros((ny, nx))
    x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))

    trace_centers = []
    y_0s = [int(2*ny/5), int(3*ny/5), int(4*ny/5)]
    blaze_function = 1 - blaze * 1e-5 * (x2d - nx / 2.) ** 2
    input_traces = np.zeros((ny, nx), dtype=np.int)
    for i in range(3):
        trace_centers.append(5e-4 * (np.arange(nx) - nx / 2.) ** 2 + y_0s[i])
        test_data += gaussian(y2d, trace_centers[i], 3, a=10000.0) * blaze_function
        input_traces[np.abs(y2d - trace_centers[i]) <= trace_half_width] = i + 1
    return test_data, trace_centers, input_traces


def test_refine_traces_curved_trace():
    nx, ny = 1001, 1003
    expected_trace = np.zeros((ny, nx), dtype=int)
    trace_half_height = 5
    x0 = 498
    input_y_center = 2e-4 * (np.arange(nx) - x0) ** 2.0 + 0.01 * (np.arange(nx) - x0) + 502.0
    trace = np.round(input_y_center).astype(int)
    for i in np.arange(nx, dtype=int):
        for j in range(-trace_half_height, trace_half_height + 1):
            expected_trace[trace[i] + j, i] = 1
    x2d, y2d = np.meshgrid(np.arange(nx, dtype=int), np.arange(ny, dtype=int))
    sigma = 1.5
    y_center = np.array([input_y_center.copy() for i in range(ny)])
    test_data = 1 / sigma / (2.0 * np.pi) * np.exp(-0.5 * (y2d - y_center) ** 2.0 / sigma ** 2)

    test_image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=np.zeros_like(test_data),
                                                              meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = np.ones_like(expected_trace)
    refine_traces(test_image, weights=test_image.data, trace_half_height=trace_half_height)

    np.testing.assert_equal(test_image.traces, expected_trace)
