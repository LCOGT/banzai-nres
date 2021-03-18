from banzai.stages import Stage
import numpy as np
from banzai_nres.utils.trace_utils import get_trace_region
from astropy.table import Table


class ProfileFitter(Stage):
    def do_stage(self, image):
        image.profile = np.zeros_like(image.data)
        blazes = np.zeros((image.num_traces, image.data.shape[1]), dtype=float)
        blazes_errors = np.zeros((image.num_traces, image.data.shape[1]), dtype=float)
        trace_ids = range(1, image.num_traces + 1)
        for i, trace_id in enumerate(trace_ids):
            # Extract the pixels from the spectrum extension in that order
            this_trace = get_trace_region(image.traces == trace_id)

            # Evaluate the best fit spline along the center of the trace, normalize, save as BLAZE
            # do a flux weighted mean in the y direction
            fluxes = image.data[this_trace]
            flux_errors = image.uncertainty[this_trace]
            blaze = (fluxes * fluxes).sum(axis=0) / fluxes.sum(axis=0)
            # Uncertainty propagation is a bear
            blaze_errors = np.sqrt(np.sum((np.sqrt(2.0) * fluxes * flux_errors) ** 2.0, axis=0))
            flux_sum_error = np.sqrt((flux_errors * flux_errors).sum(axis=0))
            blaze_errors = blaze * np.sqrt(
                (blaze_errors / (fluxes * fluxes).sum(axis=0)) ** 2.0 + (flux_sum_error / fluxes.sum(axis=0)) ** 2.0)
            # Normalize so the sum of the blaze is 1
            blaze_sum = blaze.sum()
            blaze_sum_error = np.sqrt((blaze * blaze).sum())
            blaze /= blaze_sum
            blaze_errors = blaze * np.sqrt(
                (blaze_errors / (blaze * blaze_sum)) ** 2 + (blaze_sum_error / blaze_sum) ** 2.0)
            # Evaluate the best fit spline/BLAZE save in an empty array, normalize (PROFILE extension)
            image.profile[this_trace] = fluxes / blaze
            image.profile[this_trace] /= image.profile[this_trace].sum(axis=0)

            image.profile[this_trace][fluxes < 0.0] = 0.0
            # Use the 4 bit to represent negative values in the profile
            image.mask[this_trace][fluxes < 0.0] |= 4
            # Divide the original image by best fit spline (Pixel-to-pixel sensitivity)
            # TODO: propagate the division of blaze and profile to the uncertainties of the image!!!
            image.data[this_trace] /= blaze
            image.data[this_trace] /= image.profile[this_trace]
            # Normalize the pixel-to-pixel sensitivity by the median
            image.data[this_trace] /= np.median(this_trace)
            # get the horizontal (x) extent of the trace.
            x_extent = slice(np.min(this_trace[1]), np.max(this_trace[1]) + 1)
            blazes[i, x_extent] = blaze
            blazes_errors[i, x_extent] = blaze_errors
        image.blaze = Table({'id': trace_ids, 'blaze': blazes, 'blaze_error': blazes_errors})

        return image
