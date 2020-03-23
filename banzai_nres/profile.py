from banzai.stages import Stage
import numpy as np
from banzai_nres.utils.trace_utils import get_trace_region
from astropy.table import Table


class ProfileFitter(Stage):
    def do_stage(self, image):
        image.profile = np.zeros_like(image.data)
        blazes = np.zeros((image.num_traces, image.data.shape[1]), dtype=float)
        trace_ids = range(1, image.num_traces + 1)
        for i, trace_id in enumerate(trace_ids):
            # Extract the pixels from the spectrum extension in that order
            this_trace = get_trace_region(image.traces == trace_id)

            # Evaluate the best fit spline along the center of the trace, normalize, save as BLAZE
            # do a flux weighted mean in the y direction
            fluxes = image.data[this_trace]
            blaze = (fluxes * fluxes).sum(axis=0) / fluxes.sum(axis=0)
            # Normalize so the sum of the blaze is 1
            blaze /= blaze.sum()

            # Evaluate the best fit spline/BLAZE save in an empty array, normalize (PROFILE extension)
            image.profile[this_trace] = fluxes / blaze
            image.profile[this_trace] /= image.profile[this_trace].sum(axis=0)

            # Divide the original image by best fit spline (Pixel-to-pixel sensitivity)
            image.data[this_trace] /= blaze
            image.data[this_trace] /= image.profile[this_trace]
            # Normalize the pixel-to-pixel sensitivity by the median
            image.data[this_trace] /= np.median(this_trace)
            x_extent = slice(np.min(this_trace[1]), np.max(this_trace[1]) + 1)  # get the horizontal (x) extent of the trace.
            blazes[i, x_extent] = blaze
        image.blaze = Table({'id': trace_ids, 'blaze': blazes})

        return image
