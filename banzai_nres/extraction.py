"""
Scripts for extracting spectrum as a function of x-pixel (horizontal coordinate across the CCD)
Author
    G. Mirek Brandt (gmbrandt@ucsb.edu)
"""

from banzai.stages import Stage
from banzai_nres.utils.extraction_utils import extract_spectrum_full_image, OptimalFiberProfileExtraction, BoxExtraction, VerticalExtraction
from banzai_nres.fiber_profile import FiberProfile
from banzai_nres.utils.NRES_class_utils import add_class_as_attribute
import numpy as np


class Spectra(object):
    def __init__(self):
        self.intensity_versus_x_per_order = None


class ExtractSpectrumVersusPixel(Stage):
    """
    This stage extracts the spectrum of the image versus x pixel, e.g. no wavelengths are yet assigned.
    Each method returns an list such that list[i] = flux_vs_x, x_coords     of the ith trace.
    self.extraction_window_fwhm =  number of full width half maxes you want to use to define your extraction window
    E.g. 5, for a fwhm of 1 gives a 5 pixel extraction window.
    """

    def __init__(self, pipeline_context):
        super(ExtractSpectrumVersusPixel, self).__init__(pipeline_context)
        self.extraction_method = VerticalExtraction
        self.extraction_window_fwhm = 5

    @property
    def group_by_keywords(self):
        return ['ccdsum']

    @property
    def calibration_type(self):
        return 'extracting_spectrum'

    def do_stage(self, images):
        add_class_as_attribute(images, 'fiber_profile', FiberProfile)
        add_class_as_attribute(images, 'spectra', Spectra)
        for image in images:
            fwhm = image.fiber_profile.median_full_width_half_max
            fluxes_and_x_coords = extract_spectrum_full_image(image, image.coordinates.delta_y_from_trace,
                                                                  image.coordinates.closest_trace_in_y,
                                                                  window=self.extraction_window_fwhm * fwhm,
                                                                  ExtractMethod=self.extraction_method)

            image.spectra.intensity_versus_x_per_order = \
                self.split_and_fill_spectra_with_zeros_where_missing(fluxes_and_x_coords, image)
        return images

    def split_and_fill_spectra_with_zeros_where_missing(self, fluxes_and_x_coords, image):
        fluxes = [fluxes_and_x_coords[i][0] for i in range(len(fluxes_and_x_coords))]
        x_coords = [fluxes_and_x_coords[i][1] for i in range(len(fluxes_and_x_coords))]
        full_spectrum_per_order = np.zeros((len(fluxes), image.data.shape[1]))
        for order in range(len(fluxes)):
            assert (x_coords[order] == np.arange(x_coords[order].min(), x_coords[order].max() + 1)).all()
            full_spectrum_per_order[order, x_coords[order].min():x_coords[order].max() + 1] = fluxes[order]
        return full_spectrum_per_order


class BoxExtractSpectrum(ExtractSpectrumVersusPixel):
    def __init__(self, pipeline_context):
        super(BoxExtractSpectrum, self).__init__(pipeline_context)
        self.extraction_method = BoxExtraction


class OptimallyExtractSpectrum(ExtractSpectrumVersusPixel):
    def __init__(self, pipeline_context):
        super(OptimallyExtractSpectrum, self).__init__(pipeline_context)
        self.extraction_method = OptimalFiberProfileExtraction


if __name__ == "__main__":
    None
