import numpy as np

from banzai_nres.tests.test_traces import make_simple_traces
from banzai_nres.frames import NRESObservationFrame
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.profile import ProfileFitter
from banzai import context


class TestProfileFitter:
    def test_fit(self):
        nx, ny = 401, 403
        x2d, y2d = np.meshgrid(np.arange(nx), np.arange(ny))
        test_data, trace_centers, input_traces = make_simple_traces(nx, ny)
        read_noise = 10.0

        # Filter out the numerical noise
        test_data[test_data < 1e-15] = 0.0
        test_data += np.random.poisson(test_data)
        test_data += np.random.normal(0.0, read_noise, size=test_data.shape)
        uncertainty = np.sqrt(test_data + read_noise ** 2.0)
        image = NRESObservationFrame([EchelleSpectralCCDData(data=test_data, uncertainty=uncertainty,
                                                                  meta={'OBJECTS': 'tung&tung&none'})], 'foo.fits')
        image.traces = input_traces
        input_context = context.Context({})

        stage = ProfileFitter(input_context)
        image = stage.do_stage(image)

        import matplotlib.pyplot as plt
        plt.imshow(image.profile)
        plt.show()
        import pdb; pdb.set_trace()