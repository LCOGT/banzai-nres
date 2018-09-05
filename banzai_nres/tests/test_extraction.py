from banzai.tests.utils import FakeContext
from banzai import logs

from banzai_nres.tests.utils import append_good_region_info, fill_with_simple_inverse_variances,\
    generate_image_with_two_flat_traces
from banzai_nres.coordinate_transform import MakeTraceCentricCoordinates
from banzai_nres import extraction

import numpy as np

logger = logs.get_logger(__name__)


def undefined_testing_of_extraction_method(extract_method=extraction.ExtractSpectrumVersusPixel):
    real_full_width_half_max = 1.25
    image = generate_image_with_two_flat_traces(order_width=real_full_width_half_max)
    append_good_region_info(image)

    fill_with_simple_inverse_variances(image)
    images = [image]
    # appending coordinate info
    coordinate_stage = MakeTraceCentricCoordinates(FakeContext())
    images = coordinate_stage.do_stage(images)
    fiber_profile_image = generate_image_with_two_flat_traces(order_width=real_full_width_half_max,
                                                              normalized_traces=True, add_noise=False)
    images[0].fiber_profile.normalized_fiber_profile_image = fiber_profile_image.data
    images[0].fiber_profile.median_full_width_half_max = real_full_width_half_max

    extract_stage = extract_method(FakeContext())
    extract_stage.extraction_window_fwhm = 5
    images = extract_stage.do_stage(images)
    extracted_spectra = images[0].spectra.intensity_versus_x_per_order
    for order in range(extracted_spectra.shape[0]):
        normed_order_spectrum = extracted_spectra[order]/np.mean(extracted_spectra[order])
        logger.info('mean of extracted and normalized spectrum with method %s is %s' % (extract_method,
                                                                                        np.mean(normed_order_spectrum)))
        logger.info('the standard deviation is %s' % np.std(normed_order_spectrum))
        assert np.abs(np.mean(normed_order_spectrum) - 1) < 1E-5


def test_box_extraction():
    undefined_testing_of_extraction_method(extract_method=extraction.OptimallyExtractSpectrum)


def test_optimal_extraction():
    undefined_testing_of_extraction_method(extract_method=extraction.BoxExtractSpectrum)
