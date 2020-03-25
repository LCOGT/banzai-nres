import numpy as np

from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace
from xwavecal.wavelength import wavelength_calibrate, WavelengthSolution
from xwavecal.utils.wavelength_utils import find_nearest
import logging

logger = logging.getLogger('banzai')

NRES_PRINCIPLE_ORDER_NUMBER = 52
MODELS = {'initial_wavelength_model': {1: [0, 1, 2], 2: [0, 1, 2]},
          'intermediate_wavelength_model': {0: [0, 1, 2], 1: [0, 1, 2], 2: [0, 1, 2]},
          'final_wavelength_model': {0: [0, 1, 2, 3, 4, 5],
                                     1: [0, 1, 2, 3, 4, 5],
                                     2: [0, 1, 2, 3, 4, 5],
                                     3: [0, 1, 2, 3, 4, 5],
                                     4: [0]}}


class WavelengthLoader(Stage):
    """
    Loads the wavelengths and reference_ids onto from the nearest Arc-emission lamp (wavelength calibration) in the db.
    If the traces have shifted (e.g. the instrument was serviced), then we do not load wavelengths nor reference ids
    (new wavelengths will be created from scratch in WavelengthCalibrate).
    """
    def do_stage(self, image: EchelleSpectralCCDData):
        if not image.traces.blind_solve:
            # load in ref_id column.
            # load in wavelengths


class WavelengthCalibrate(Stage):
    """
    Stage to recalibrate wavelength-calibration (e.g. arc lamp) frames.
    We re-wavelength calibrate from scratch if the Arc does not have reference id's assigned to each trace.
    We lightly recalibrate the wavelength calibration if the arc has reference id's assigned to each trace.
    """
    def do_stage(self, image: EchelleSpectralCCDData):
        # image.spectrum = Table({'id': trace_ids, 'flux': flux, 'uncertainty': np.sqrt(variance)})
        # identify emission lines and get (pixel, order) positions.
        features = identify_features(image.data, image.uncertainty, image.mask, nsigma=5.0, fwhm=6.0)
        # helps to have flat fielded emission line fluxes.
        fiber, ref_id, measured_lines, line_list = None, None, None, None
        pixel, order = np.arange(image.data.shape[1]), np.arange(len(image.spectrum[image.spectrum['fiber'] == fiber]))
        if ref_id is None:
            # reidentify reference ids.
            # blind solve for the wavelengths of the emission lines
            measured_lines['wavelength'] = wavelength_calibrate(measured_lines, line_list, pixel, order,
                                                                principle_order_number=NRES_PRINCIPLE_ORDER_NUMBER)
        else:
            # adopt wavelengths for the emission lines from the existing solution, good to 1 pixel.
            # Note: we could improve this considerably with ndimage.map_coordinates(..., order=1, prefilter=False).
            measured_lines['wavelength'] = image.wavelengths[measured_lines['order'].astype(int),
                                                             measured_lines['pixel'].astype(int)]
        wavelength_solution = self.recalibrate(measured_lines, line_list, pixel, order)
        # do some qc checks. Throw warnings necessary.
        return image

    def recalibrate(self, measured_lines, line_list, pixel, order):
        """
        Recalibrate an existing wavelength solution using new measured spectral features.
        :param measured_lines:
        :param line_list:
        :param pixel:
        :param order:
        :return: wavelength_solution. xwavecal.wavelength.WavelengthSolution.
        wavelength_solution(x, order) will give the wavelength of the len(x) points with pixel and order coordinates
        x and order, respectively.
        where x and order are ndarray's of the same shape.

        """
        wavelength_solution = WavelengthSolution(model=MODELS.get('final_wavelength_model'),
                                                 min_order=np.min(order), max_order=np.max(order),
                                                 min_pixel=np.min(pixel), max_pixel=np.max(pixel),
                                                 measured_lines=measured_lines, reference_lines=line_list,
                                                 m0=NRES_PRINCIPLE_ORDER_NUMBER)
        wavelengths_to_fit = find_nearest(measured_lines['wavelength'], np.sort(line_list))
        weights = np.zeros_like(wavelengths_to_fit, dtype=float)
        # fit lines who have less than 0.1 angstrom error
        weights[np.isclose(wavelengths_to_fit, measured_lines['wavelength'], atol=0.1)] = 1
        wavelength_solution.model_coefficients = wavelength_solution.solve(measured_lines, wavelengths_to_fit, weights=weights)
        return wavelength_solution
