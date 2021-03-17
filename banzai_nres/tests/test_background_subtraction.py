from banzai_nres.frames import NRESObservationFrame
from banzai.data import CCDData
import numpy as np
from banzai_nres.background import BackgroundSubtractor
from banzai import context
import pytest


@pytest.fixture
def seed():
    np.random.seed(11248)


def test_background_subtraction_on_noisy_data(seed):
    nx, ny = 405, 403
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    noise_sigma = 1.0
    input_background = 30 * np.exp(-(X - nx / 2.0)**2/300**2 - (Y - ny / 2.0 - 50.0)**2 / 200**2)
    test_data = input_background + np.random.normal(0.0, noise_sigma, size=input_background.shape)
    test_image = NRESObservationFrame([CCDData(data=test_data, uncertainty=np.ones((ny, nx)) * noise_sigma,
                                      meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.traces = np.zeros((ny, nx))
    input_context = context.Context({})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    # Make sure our background estimation is good. This is roughly 4 counts which is not bad
    # If we fully model the image, we can do better than this, but that becomes computationally prohibitive on a 4kx4k
    # image
    np.testing.assert_allclose(output_image.background, input_background, atol=4.0)


def test_background_subtraction_with_traces(seed):
    nx, ny = 405, 403
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    noise_sigma = 1.0
    input_background = 30 * np.exp(-(X - nx / 2.0)**2/300**2 - (Y - ny / 2.0 - 50.0)**2 / 200**2)
    test_data = input_background + np.random.normal(0.0, noise_sigma, size=input_background.shape)
    test_image = NRESObservationFrame([CCDData(data=test_data, uncertainty=np.ones((ny, nx)) * noise_sigma,
                                      meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')

    test_image.traces = np.zeros((ny, nx))

    for i in range(1, 8):
        test_image.traces[40*i:40*i+10] = i

    input_context = context.Context({})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    # Make sure our background estimation is good. This is roughly 4 counts which is not bad
    # If we fully model the image, we can do better than this, but that becomes computationally prohibitive on a 4kx4k
    # image
    np.testing.assert_allclose(output_image.background, input_background, atol=5.0)
