from banzai.calibrations import CalibrationComparer
from banzai.utils import qc
import numpy as np


class NRESCalibrationComparer(CalibrationComparer):
    SIGNAL_TO_NOISE_THRESHOLD = 10.0

    def is_frame_bad(self, image, master_calibration_image):
        orders = np.unique(image.traces)
        orders = orders[orders > 0]
        order_quality = []
        for order in orders:
            in_order = image.traces == order
            data = image.data[in_order]
            error = image.uncertainty[in_order]
            reference_data = master_calibration_image.data[in_order]
            reference_error = master_calibration_image.data[in_order]

            # Normalize the data using only data above the SNR cutoff
            data_region = data / error > self.SIGNAL_TO_NOISE_THRESHOLD
            data_normalizaton = np.median(data[data_region])

            data /= data_normalizaton
            error /= data_normalizaton

            reference_region = reference_data / reference_error > self.SIGNAL_TO_NOISE_THRESHOLD
            reference_normalization = np.median(reference_data[reference_region])

            reference_data /= reference_normalization
            reference_error /= reference_normalization

            # Compare pixels if either reference or image has a SNR larger than the cutoff
            comparison_region = np.logical_or(data_region, reference_region)

            difference = data[comparison_region] - reference_data[comparison_region]
            difference_error = np.sqrt(error[comparison_region] ** 2.0 + reference_error[comparison_region] ** 2.0)
            outliers = (difference / difference_error) >= self.SIGNAL_TO_NOISE_THRESHOLD
            bad_pixel_fraction = outliers.sum() / float(outliers.size)
            order_quality.append(bad_pixel_fraction > self.ACCEPTABLE_PIXEL_FRACTION)

        bad_pixel_fraction = np.sum(order_quality) / float(len(order_quality))
        frame_is_bad = bad_pixel_fraction >= self.ACCEPTABLE_PIXEL_FRACTION

        qc_results = {"master_comparison.fraction": bad_pixel_fraction,
                      "master_comparison.snr_threshold": self.SIGNAL_TO_NOISE_THRESHOLD,
                      "master_comparison.pixel_threshold": self.ACCEPTABLE_PIXEL_FRACTION,
                      "master_comparison.comparison_master_filename": master_calibration_image.filename,
                      "master_comparison.failed": frame_is_bad}

        qc.save_qc_results(self.runtime_context, qc_results, image)
        return frame_is_bad
