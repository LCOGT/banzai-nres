import pytest
import mock
import banzai_nres.tests.utils as test_utils
from banzai_nres.traces import BlindTraceMaker
from banzai.tests.utils import FakeContext
import numpy as np
from banzai import logs


logger = logs.get_logger(__name__)


@mock.patch('banzai_nres.traces.Image')
def test_blind_trace_maker(mock_images):
    nimages = 3
    # in counts
    readnoise = 11.0

    for x in range(nimages):
        images = [test_utils.FakeTraceImage(readnoise=readnoise)]

        test_utils.make_random_yet_realistic_trace_coefficients(images[0])
        test_utils.fill_image_with_traces(images[0])
        test_utils.add_overscan_and_noisify_image(images[0])

        maker = BlindTraceMaker(FakeContext())
        maker.do_stage(images)

        args, kwargs = mock_images.call_args
        master_trace = kwargs['data']

        difference = test_utils.differences_between_found_and_generated_trace_vals(master_trace, images[0])
        logger.info('test_error in trace fitting is less than %s of a pixel' % np.median(np.abs(difference)))
        assert np.median(np.abs(difference)) < 1/10