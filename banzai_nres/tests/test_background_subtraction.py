from banzai_nres.images import NRESObservationFrame
from banzai.images import CCDData
import numpy as np
from banzai_nres.background import BackgroundSubtractor
from banzai import context


def test_background_subtraction_on_noisy_data():
    nx, ny = 105, 103
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    noise_sigma = 1.0
    input_background = 30 * np.exp(-(X - nx / 2.0)**2/10**2 - Y**2/45**2)
    test_data = input_background + np.random.normal(0.0, noise_sigma, size=input_background.shape)
    test_image = NRESObservationFrame([CCDData(data=test_data, uncertainty=np.ones((ny, nx)) * 3.0,
                                      meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.trace = np.zeros((ny, nx))
    input_context = context.Context({'background_smoothing_scale': 15, 'background_spline_order': 2})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    np.testing.assert_allclose(output_image.background.data, input_background, atol=3.0)


def test_background_subtraction_with_traces():
    nx, ny = 105, 103
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    noise_sigma = 1.0
    input_background = 30 * np.exp(-(X - nx / 2.0)**2/10**2 - Y**2/45**2)
    test_data = input_background + np.random.normal(0.0, noise_sigma, size=input_background.shape)
    test_image = NRESObservationFrame([CCDData(data=test_data, uncertainty=np.ones((ny, nx)) * 3.0,
                                      meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    test_image.trace = np.zeros((ny, nx))
    test_image.trace[20:30] = 1
    input_context = context.Context({'background_smoothing_scale': 15, 'background_spline_order': 2})
    stage = BackgroundSubtractor(input_context)
    output_image = stage.do_stage(test_image)
    np.testing.assert_allclose(output_image.background.data, input_background, atol=3.0)

