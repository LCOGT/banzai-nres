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
from numpy.polynomial.polynomial import Polynomial, polyfit

logger = logging.getLogger('banzai')

# Minimum separation between peaks in the image that could be separate traces
MIN_TRACE_SEPARATION = 10

# Cutoff in signal-to-noise to stop following a trace
SIGNAL_TO_NOISE_TRACING_CUTOFF = 10

# Maximum degree of polynomial which describes the y position of the trace as a function of pixel
POLY_DEGREE = 5

# The final trace will be +- this from the center in the y-direction
TRACE_HALF_HEIGHT = 6


def find_y_center(y, indices, weights=None):
    if weights is None:
        weights = 1
    centers = (y * indices * weights).sum(axis=0) / (indices * weights).sum(axis=0)
    # Errors are sqrt sum of the weights squared
    errors = np.sqrt((indices * weights ** 2).sum(axis=0))
    return centers, errors


def refine_traces(traces, image, weights=None):
    x2d, y2d = np.meshgrid(np.arange(traces.shape[1]), np.arange(traces.shape[0]))
    # For each label
    for i in range(1, np.max(traces) + 1):
        y_center, y_center_errors = find_y_center(y2d, traces == i, weights=weights)
        # Refit the centroids to reject cosmic rays etc, but only evaluate where the S/N is good
        x_center = np.arange(min(x2d[traces == i]), max(x2d[traces == i]) + 1, dtype=np.float)
        logger.info(f'Fitting a polynomial to trace {i}', image=image)
        best_fit = Polynomial(polyfit(x_center, y_center[x_center.astype(int)], deg=POLY_DEGREE,
                                      w=1/y_center_errors[x_center.astype(int)]**2))
        y_center = best_fit(x_center)

        # Pad y_center with zeros so that it has the same dimension as y2d
        padded_y_center = np.zeros(y2d.shape[1])
        padded_y_center[int(min(x_center)):int(max(x_center) + 1)] = y_center[:]
        pixels_to_label = np.logical_and(np.logical_and(x2d >= min(x_center), x2d <= max(x_center)),
                                         np.abs(y2d - padded_y_center) <= TRACE_HALF_HEIGHT)
        # Reset the previously marked traces to 0. Then mark the newly measured traces.
        traces[traces == i] = 0
        traces[pixels_to_label] = i
    return traces


class TraceInitializer(Stage):
    # Minimum half width of a feature to be considered a trace
    min_trace_half_width = 500

    def do_stage(self, image):
        if image.traces is None:
            image.traces = self.blind_solve(image)
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
        true_labels = labeled_indices[np.logical_and(x_maxes > (image.shape[1] // 2 + self.min_trace_half_width),
                                                     x_mins < (image.shape[1] // 2 - self.min_trace_half_width))]
        # Reset the values that are not actually in traces
        labeled_image[np.logical_not(np.isin(labeled_image, true_labels))] = 0

        # Reindex the labels to start at 1
        for i, label in enumerate(true_labels):
            labeled_image[labeled_image == label] = i + 1
        return refine_traces(labeled_image, image, weights=labeled_image > 0)


class TraceRefiner(Stage):
    def do_stage(self, image):
        image.traces = refine_traces(image.traces, image, weights=image.data)
        return image
