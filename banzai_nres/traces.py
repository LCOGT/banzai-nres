"""
traces.py: Stages for finding traces in an NRES frame

Authors
    Curtis McCully
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from banzai.stages import Stage
from banzai.logs import get_logger
from scipy import ndimage
import numpy as np
from banzai_nres.fitting import fit_polynomial
from banzai.data import ArrayData

logger = get_logger()

# Minimum separation between peaks in the image that could be separate traces
MIN_TRACE_SEPARATION = 10

# Cutoff in signal-to-noise to stop following a trace
SIGNAL_TO_NOISE_TRACING_CUTOFF = 10

# Minimum half width of a feature to be considered a trace
MIN_TRACE_HALF_WIDTH = 50

# Trace half height in pixels. The final trace will be +- this from the center in the y-direction.
TRACE_HALF_HEIGHT = 5
# NOTE: 5 pixels covers the NRES orders more than sufficiently. Changing this value may break
# the tracing part of the pipeline.


def find_y_center(y, mask, weights):
    centers = (y * mask * weights).sum(axis=0) / (mask * weights).sum(axis=0)
    # Errors are inv sqrt sum of the weights
    inverse_variance = (mask * weights).sum(axis=0)
    errors = np.zeros_like(centers)
    errors[inverse_variance > 0] = 1.0 / np.sqrt(inverse_variance[inverse_variance > 0])
    return centers, errors


def refine_traces(image, weights=None, trace_half_height=5):
    x2d, y2d = np.meshgrid(np.arange(image.traces.shape[1], dtype=int), np.arange(image.traces.shape[0], dtype=int))
    if weights is None:
        weights = np.ones_like(image.data.shape)
    # For each label
    for i in range(1, np.max(image.traces) + 1):
        x_stamp = slice(min(x2d[image.traces == i]), max(x2d[image.traces == i]) + 1, 1)
        y_stamp = slice(min(y2d[image.traces == i]), max(y2d[image.traces == i]) + 1, 1)
        stamp_weights = weights[y_stamp, x_stamp] * (image.mask[y_stamp, x_stamp] == 0)
        y_center, y_center_errors = find_y_center(y2d[y_stamp, x_stamp], (image.traces == i)[y_stamp, x_stamp],
                                                  weights=stamp_weights)
        # Refit the centroids to reject cosmic rays etc, but only evaluate where the S/N is good
        x_center = np.arange(min(x2d[image.traces == i]), max(x2d[image.traces == i]) + 1, dtype=float)
        logger.info(f'Fitting a polynomial to order {i}', image=image)
        # order 5 chosen based on visually inspecting the residuals between the trace centers and the model fit centers
        # TODO we need to verify that an order 5 polynomial fit is the best thing to do.
        best_fit = fit_polynomial(y_center, y_center_errors, mask=y_center_errors == 0, x=x_center, order=5)
        y_center = best_fit(x_center)

        y_center = np.round(y_center).astype(int)
        # Pad y_center with zeros so that it has the same dimension as y2d
        padded_y_center = np.zeros(y2d.shape[1])
        padded_y_center[int(min(x_center)):int(max(x_center) + 1)] = y_center[:]
        pixels_to_label = np.logical_and(np.logical_and(x2d >= min(x_center), x2d <= max(x_center)),
                                         np.abs(y2d - padded_y_center) <= trace_half_height)
        # Reset the previously marked traces to 0. Then mark the newly measured traces.
        image.traces[image.traces == i] = 0
        image.traces[pixels_to_label] = i


class TraceInitializer(Stage):
    def do_stage(self, image):
        if image.traces is None:
            image.add_or_update(ArrayData(self.blind_solve(image, TRACE_HALF_HEIGHT), name='TRACES'))
            refine_traces(image, weights=(image.traces > 0).astype(np.float32),
                          trace_half_height=TRACE_HALF_HEIGHT)
        else:
            image.add_or_update(ArrayData(image.traces, name='TRACES'))
        return image

    @staticmethod
    def blind_solve(image, trace_half_height=5):
        # Find the peaks of each of the traces using a max filter
        peaks = ndimage.maximum_filter1d(image.data.data,
                                         size=MIN_TRACE_SEPARATION, axis=0)
        significant = (image.data.data / image.uncertainty) > SIGNAL_TO_NOISE_TRACING_CUTOFF
        # ignore pixels in the bpm
        significant = np.logical_and(significant, image.mask == 0)
        # identify the traces.
        binary_map = np.logical_and(peaks == image.data.data, significant)

        # Dilate the label map to make sure all traces are connected
        binary_map = ndimage.binary_dilation(binary_map)
        labeled_image, n_labels = ndimage.label(binary_map)
        X, Y = np.meshgrid(np.arange(image.shape[1], dtype=int), np.arange(image.shape[0], dtype=int))
        labeled_indices = np.arange(1, n_labels + 1)

        # Find the widths of the traces by finding the min and max x position
        x_maxes = ndimage.labeled_comprehension(X, labeled_image, labeled_indices, np.max, float, None).astype(int)
        x_mins = ndimage.labeled_comprehension(X, labeled_image, labeled_indices, np.min, float, None).astype(int)
        # Pick out only features that are wide like traces and span the center
        # Note labeled_indices is one indexed
        trace_xextents_ok = np.logical_and(x_maxes > (image.shape[1] // 2 + MIN_TRACE_HALF_WIDTH),
                                           x_mins < (image.shape[1] // 2 - MIN_TRACE_HALF_WIDTH))
        # and remove any traces whose centers are close (by a half width) to the top or bottom of the detector
        y_maxes = ndimage.labeled_comprehension(Y, labeled_image, labeled_indices, np.max, float, None).astype(int)
        y_mins = ndimage.labeled_comprehension(Y, labeled_image, labeled_indices, np.min, float, None).astype(int)
        trace_centers_not_near_edge = np.logical_and(y_maxes < (image.shape[0] - trace_half_height),
                                                     y_mins > trace_half_height)
        true_labels = labeled_indices[np.logical_and(trace_xextents_ok, trace_centers_not_near_edge)]
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
        refine_traces(image, weights=image.data,
                      trace_half_height=TRACE_HALF_HEIGHT)
        return image
