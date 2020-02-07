from banzai_nres.images import NRESObservationFrame
from banzai.images import CCDData
import numpy as np
from banzai_nres.background import BackgroundSubtractor
from banzai import context
from scipy.stats import anderson
import pytest


@pytest.fixture
def seed():
    np.random.seed(11248)


def test_background_subtraction_on_noisy_data(seed):
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
    input_context = context.Context({})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    # Make sure our background estimation is good. This is roughly 1.5 sigma which is not bad
    np.testing.assert_allclose(output_image.background, input_background, atol=2.0)
    # Check that the remaining image basically looks like noise
    difference = test_data - output_image.background
    anderson_test = anderson(difference.ravel(), 'norm')
    assert np.all(anderson_test.statistic < anderson_test.critical_values)


def test_background_subtraction_with_traces(seed):
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
    input_context = context.Context({})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    # Make sure our background estimation is good. This is roughly 1.5 sigma which is not bad
    np.testing.assert_allclose(output_image.background, input_background, atol=2.0)
    # Check that the remaining image basically looks like noise
    difference = test_data - output_image.background
    anderson_test = anderson(difference.ravel(), 'norm')
    assert np.all(anderson_test.statistic < anderson_test.critical_values)
