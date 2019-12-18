"""
traces.py: Stages for finding traces in an NRES frame

Authors
    Curtis McCully
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""
from banzai.stages import Stage
from banzai.images import ArrayData
import logging
from scipy import ndimage
import numpy as np

logger = logging.getLogger('banzai')


def find_y_center(y, indices, weights=None):
    if weights is None:
        weights = np.ones(y.shape, dtype=y.dtype)
    return (y * indices * weights).sum(axis=0) / (indices * weights).sum(axis=0)


def refine_traces(runtime_context, traces, weights=None):
    traces = np.zeros(traces.shape, dtype=np.uint8)
    X, Y = np.meshgrid(np.arange(traces.shape[1]), np.arange(traces.shape[0]))
    # For each label
    for i in range(1, np.max(traces) + 1):
        y_center = find_y_center(Y, traces.data == i, weights=weights)
        pixels_to_label = np.logical_and(X[traces == i],
                                         np.abs(Y - y_center) <= runtime_context.trace_half_width)
        traces[pixels_to_label] = i
    return traces


class TraceInitializer(Stage):
    def do_stage(self, image):
        if image.trace is None:
            image.trace = self.blind_solve(image)
        return image

    def blind_solve(self, image):
        # Find the peaks of each of the traces using a max filter
        peaks = ndimage.maximum_filter1d(image.data.data,
                                         size=self.runtime_context.trace_separation, axis=0)
        signal_to_noise = image.data.data / image.data.uncertainty > self.runtime_context.signal_to_noise_tracing_cutoff
        binary_map = np.logical_and(peaks == image.data.data, signal_to_noise)

        # Dilate the label map to make sure all traces are connected
        binary_map = ndimage.morphology.binary_dilation(binary_map)
        labeled_image, n_labels = ndimage.label(binary_map)
        X, Y = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
        labeled_indices = np.arange(1, n_labels + 1)

        # Find the widths of the traces by finding the min and max x position
        x_maxes = ndimage.labeled_comprehension(X, labeled_image, labeled_indices, np.max, float, None)
        x_mins = ndimage.labeled_comprehension(X, labeled_image, labeled_indices, np.min, float, None)
        true_labels = []
        for i in labeled_indices:
            # Pick out only features that are wide like traces and span the center
            if x_maxes[i] > (image.shape[1] // 2 + self.runtime_context.min_trace_half_width)\
                    and x_mins[i] < (image.shape[1] // 2 - self.runtime_context.min_trace_half_width):
                true_labels.append(i)

        for i in labeled_indices:
            if i not in true_labels:
                labeled_image[labeled_image == i] = 0
        # Reindex the labels to start at 1
        for i, label_index in enumerate(true_labels):
            labeled_image[labeled_image == label_index] = i + 1

        traces = refine_traces(self.runtime_context, labeled_image, weights=labeled_image > 0)
        return ArrayData(data=traces, meta={}, name='trace')


class TraceRefiner(Stage):
    def do_stage(self, image):
        image.trace = refine_traces(self.runtime_context, image.trace, weights=image.data)
        return image
