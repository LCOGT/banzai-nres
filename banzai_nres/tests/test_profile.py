import numpy as np

from banzai_nres.tests.test_traces import make_simple_traces
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.profile import ProfileFitter
from banzai import context


def test_profile_fit_without_noise_without_blaze():
    nx, ny = 401, 403
    input_profile, trace_centers, input_traces = make_simple_traces(nx, ny, blaze=False)
    # Filter out the numerical noise
    input_profile[input_profile < 1e-15] = 0.0
    uncertainty = np.sqrt(input_profile)
    image = NRESObservationFrame([EchelleSpectralCCDData(data=input_profile.copy(), uncertainty=uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    image.traces = input_traces
    input_context = context.Context({})

    stage = ProfileFitter(input_context)
    image = stage.do_stage(image)
    scale_factor = np.max(input_profile) / np.max(image.profile)
    assert np.allclose(input_profile[input_traces != 0], image.profile[input_traces != 0] * scale_factor, rtol=5e-3, atol=1.0)


def test_profile_fit_with_noise_with_blaze():
    nx, ny = 401, 403
    input_profile, trace_centers, input_traces = make_simple_traces(nx, ny, blaze=True)
    read_noise = 10.0

    # Filter out the numerical noise
    input_profile[input_profile < 1e-15] = 0.0
    input_profile += np.random.poisson(input_profile)
    input_profile += np.random.normal(0.0, read_noise, size=input_profile.shape)
    uncertainty = np.sqrt(input_profile + read_noise ** 2.0)
    image = NRESObservationFrame([EchelleSpectralCCDData(data=input_profile.copy(), uncertainty=uncertainty,
                                                         meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
    image.traces = input_traces
    input_context = context.Context({})

    stage = ProfileFitter(input_context)
    image = stage.do_stage(image)
    scale_factor = np.max(input_profile) / np.max(image.profile)
    assert np.allclose(input_profile[input_traces != 0], image.profile[input_traces != 0] * scale_factor, rtol=5e-3, atol=1.0)