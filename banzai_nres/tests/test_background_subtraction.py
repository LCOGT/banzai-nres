from banzai_nres.images import NRESObservationFrame
from banzai.images import CCDData
import numpy as np
from banzai_nres.background import BackgroundSubtractor
from banzai import context
from scipy.stats import anderson


def test_background_subtraction_on_noisy_data():
    nx, ny = 105, 103
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    noise_sigma = 1.0
    input_background = 30 * np.exp(-(X - nx / 2.0)**2/10**2 - Y**2/45**2)
    test_data = input_background + np.random.normal(0.0, noise_sigma, size=input_background.shape)
    test_image = NRESObservationFrame([CCDData(data=test_data, uncertainty=np.ones((ny, nx)) * noise_sigma,
                                      meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = np.zeros((ny, nx))
    input_context = context.Context({'background_smoothing_scale': 30, 'background_spline_order': 1})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    # Make sure our background estimation is good. This is roughly 1.5 sigma which is not bad
    np.testing.assert_allclose(output_image.background.data, input_background, atol=1.5)
    # Check that the remaining image basically looks like noise
    difference = output_image.data - output_image.background
    anderson_test = anderson(difference.ravel(), 'norm')
    assert np.all(anderson_test.statistic < anderson_test.critical_values)


def test_background_subtraction_with_traces():
    nx, ny = 105, 103
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    noise_sigma = 1.0
    input_background = 30 * np.exp(-(X - nx / 2.0)**2/10**2 - Y**2/45**2)
    test_data = input_background + np.random.normal(0.0, noise_sigma, size=input_background.shape)
    test_image = NRESObservationFrame([CCDData(data=test_data, uncertainty=np.ones((ny, nx)) * noise_sigma,
                                      meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = np.zeros((ny, nx))
    test_image.traces[20:30] = 1
    input_context = context.Context({'background_smoothing_scale': 15, 'background_spline_order': 2})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    # Make sure our background estimation is good. This is roughly 1.5 sigma which is not bad
    np.testing.assert_allclose(output_image.background, input_background, atol=1.5)
    # Check that the remaining image basically looks like noise
    difference = output_image.data - output_image.background.data
    anderson_test = anderson(difference.ravel(), 'norm')
    assert np.all(anderson_test.statistic < anderson_test.critical_values)
