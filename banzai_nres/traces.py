"""
traces.py: Stages for finding traces in an NRES frame

Authors
    Curtis McCully
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from banzai.stages import Stage
import logging
from scipy import ndimage
import numpy as np
from banzai_nres.fitting import fit_polynomial

logger = logging.getLogger('banzai')

# Minimum separation between peaks in the image that could be separate traces
MIN_TRACE_SEPARATION = 10

# Cutoff in signal-to-noise to stop following a trace
SIGNAL_TO_NOISE_TRACING_CUTOFF = 50

# Minimum half width of a feature to be considered a trace
MIN_TRACE_HALF_WIDTH = 50

# The final trace will be +- this from the center in the y-direction
TRACE_HALF_HEIGHT = 7


def find_y_center(y, indices, weights=None):
    if weights is None:
        weights = np.ones(y.shape, dtype=y.dtype)
    centers = (y * indices * weights).sum(axis=0) / (indices * weights).sum(axis=0)
    # Errors are sqrt sum of the weights
    errors = 1.0 / np.sqrt((indices * weights).sum(axis=0))
    return centers, errors


def refine_traces(image, weights=None):
    x2d, y2d = np.meshgrid(np.arange(image.traces.shape[1]), np.arange(image.traces.shape[0]))
    # For each label
    for i in range(1, np.max(image.traces) + 1):
        x_stamp = slice(min(x2d[image.traces == i]), max(x2d[image.traces == i]) + 1, 1)
        y_stamp = slice(min(y2d[image.traces == i]), max(y2d[image.traces == i]) + 1, 1)
        y_center, y_center_errors = find_y_center(y2d[y_stamp, x_stamp], (image.traces == i)[y_stamp, x_stamp],
                                                  weights=weights[y_stamp, x_stamp])
        # Refit the centroids to reject cosmic rays etc, but only evaluate where the S/N is good
        x_center = np.arange(min(x2d[image.traces == i]), max(x2d[image.traces == i]) + 1, dtype=np.float)
        logger.info(f'Fitting a polynomial to order {i}', image=image)
        # we chose order 5 based on visually inspecting the residuals between the trace centers and the model fit centers
        # TODO we need to verify that an order 5 polynomial fit is the best thing to do.
        best_fit = fit_polynomial(y_center[x_center.astype(int)], y_center_errors[x_center.astype(int)], x=x_center, order=5)
        # TODO: extrapolate fit and go out to pixels that have a cumulative S/N = 10
        y_center = best_fit(x_center)

        y_center = np.round(y_center).astype(int)
        # Pad y_center with zeros so that it has the same dimension as y2d
        padded_y_center = np.zeros(y2d.shape[1])
        padded_y_center[int(min(x_center)):int(max(x_center) + 1)] = y_center[:]
        pixels_to_label = np.logical_and(np.logical_and(x2d >= min(x_center), x2d <= max(x_center)),
                                         np.abs(y2d - padded_y_center) <= TRACE_HALF_HEIGHT)
        # Reset the previously marked traces to 0. Then mark the newly measured traces.
        image.traces[image.traces == i] = 0
        image.traces[pixels_to_label] = i


class TraceInitializer(Stage):
    def do_stage(self, image):
        if image.traces is None:
            image.traces = self.blind_solve(image)
            refine_traces(image, weights=image.traces > 0)
        return image

    def blind_solve(self, image):
        # Find the peaks of each of the traces using a max filter
        peaks = ndimage.maximum_filter1d(image.data.data,
                                         size=MIN_TRACE_SEPARATION, axis=0)
        signal_to_noise = image.data.data / image.uncertainty > SIGNAL_TO_NOISE_TRACING_CUTOFF
        binary_map = np.logical_and(peaks == image.data.data, signal_to_noise)

        # Dilate the label map to make sure all traces are connected
        binary_map = ndimage.morphology.binary_dilation(binary_map)
        labeled_image, n_labels = ndimage.label(binary_map)
        X, Y = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
        labeled_indices = np.arange(1, n_labels + 1)

        # Find the widths of the traces by finding the min and max x position
        x_maxes = ndimage.labeled_comprehension(X, labeled_image, labeled_indices, np.max, float, None)
        x_mins = ndimage.labeled_comprehension(X, labeled_image, labeled_indices, np.min, float, None)
        # Pick out only features that are wide like traces and span the center
        # Note labeled_indices is one indexed
        true_labels = labeled_indices[np.logical_and(x_maxes > (image.shape[1] // 2 + MIN_TRACE_HALF_WIDTH),
                                                     x_mins < (image.shape[1] // 2 - MIN_TRACE_HALF_WIDTH))]
        # Reset the values that are not actually in traces
        labeled_image[np.logical_not(np.isin(labeled_image, true_labels))] = 0

        # Reindex the labels to start at 1
        for i, label in enumerate(true_labels):
            labeled_image[labeled_image == label] = i + 1
        return labeled_image


# TODO do flux weighted mean residuals across 10, 11, 20,21 30,31,.. 40, ..., 130
#  use get_trace_region to check the residuals and that order 5 is correct.


class TraceRefiner(Stage):
    def do_stage(self, image):
        refine_traces(image, weights=image.data)
        return image
