from banzai.stages import Stage
from banzai_nres import fitting
import numpy as np
from banzai_nres.utils.trace_utils import get_trace_region


class ProfileFitter(Stage):
    def do_stage(self, image):
        image.profile = np.zeros_like(image.data)
        x2d, y2d = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
        for trace_id in range(1, image.num_traces + 1):
            # Extract the pixels from the spectrum extension in that order
            this_trace = get_trace_region(image.traces == trace_id)

            # Fit a smooth B-spline to the pixels
            best_fit = fitting.fit_smooth_spline(image.data[this_trace], image.uncertainty[this_trace],
                                                 mask=image.mask[this_trace], x=(x2d[this_trace], y2d[this_trace]))
            # Evaluate the best fit spline along the center of the trace, normalize, save as BLAZE
            # do a flux weighted mean in the y direction
            fluxes = best_fit(x2d[this_trace], y2d[this_trace], grid=False)
            import pdb; pdb.set_trace()
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
        return image
