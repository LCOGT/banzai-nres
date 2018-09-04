"""
Scripts for assigning to each pixel:
    1. The index number of the closest trace
    2. The distance in y (perpendicular to dispersion) from said closest-trace.
Author
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai.stages import Stage
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute
from banzai_nres.utils.coordinate_utils import generate_delta_y_and_closest_trace_coordinates

import numpy as np


class MakeTraceCentricCoordinates(Stage):
    """
    This stage identifies the closest trace to each pixel in the image, and stores how far away (in y pixels)
    that pixel is from said nearest-trace.
    """
    def __init__(self, pipeline_context):
        super(MakeTraceCentricCoordinates, self).__init__(pipeline_context)

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'coordinate_transformation'

    def do_stage(self, images):
        """
        :param images: list of banzai-nres images with trace.coefficients
        Appends the array (size image.data) which gives each pixel (x,y) the index of the trace it is closest to,
        e.g. trace number 67 would be the 0th order of the second fiber (for two fibers and 67 orders each).
        Also appends the difference in y coordinate between the pixel (x,y) and the location of the trace (x, y_trace)
        y - y_trace, such that delta_y > 0 implies the point lies above the corresponding pixel.
        """
        add_class_as_attribute(images, 'coordinates', Coordinates)
        prototypical_image = images[0]
        generate_delta_y_and_closest_trace_coordinates(prototypical_image)
        for image in images:
            if (image.trace.coefficients == prototypical_image.trace.coefficients).all():
                image.coordinates.delta_y_from_trace = prototypical_image.coordinates.delta_y_from_trace
                image.coordinates.closest_trace_in_y = prototypical_image.coordinates.closest_trace_in_y
            else:
                generate_delta_y_and_closest_trace_coordinates(image)
            X, Y = np.meshgrid(np.arange(image.data.shape[1]), np.arange(image.data.shape[0]))
            image.coordinates.x, image.coordinates.y = X, Y
        return images


class Coordinates(object):
    """
    Object for storing all the coordinate related information, such as the x,y coordinates of each pixel in image.data
    (trivial example). This gets appended to each Image instance.
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.delta_y_from_trace = None
        self.closest_trace_in_y = None
