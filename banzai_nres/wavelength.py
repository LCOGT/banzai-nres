import numpy as np

from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace, aperture_extract
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
EMISSION_FEATURE_SIZE = 6.0 # pixels


class WavelengthLoader(Stage):
    """
    Loads the wavelengths and reference_ids onto from the nearest Arc-emission lamp (wavelength calibration) in the db.
    If the traces have shifted (e.g. the instrument was serviced), then we do not load wavelengths nor reference ids
    (new wavelengths will be created from scratch in WavelengthCalibrate).
    """
    def do_stage(self, image: EchelleSpectralCCDData):
        if not image.traces.blind_solve:
            pass
            # load in ref_id column.
            # load in wavelengths


class WavelengthCalibrate(Stage):
    """
    Stage to recalibrate wavelength-calibration (e.g. arc lamp) frames.
    We re-wavelength calibrate from scratch if the Arc does not have reference id's assigned to each trace.
    We lightly recalibrate the wavelength calibration if the arc has reference id's assigned to each trace.
    """
    def do_stage(self, image: EchelleSpectralCCDData):
        # load ref_id, id, fibers from Arc lamp
        ref_id = None
        id = None
        fibers = None
        for fiber in list(set(fibers)):
            pixel, order = np.arange(image.data.shape[1]), np.sort(ref_id[fibers == fiber])
            if ref_id is None:
                # reidentify reference ids.
                # blind solve for the wavelengths of the emission lines
                measured_lines['wavelength'] = wavelength_calibrate(measured_lines, line_list, pixel, order,
                                                                    principle_order_number=NRES_PRINCIPLE_ORDER_NUMBER)
            else:
                # assign each feature the correct reference id (order) and fiber:
                image.features['order'], image.features['fiber'] = ref_id[image.features['id'] - 1], fibers[
                    image.features['id'] - 1]
                # adopt wavelengths for the emission lines from the existing solution, good to 1 pixel.
                # Note: we could improve this considerably with ndimage.map_coordinates(..., order=1, prefilter=False).
                measured_lines['wavelength'] = image.wavelengths[measured_lines['order'].astype(int),
                                                                 measured_lines['pixel'].astype(int)]
            wavelength_solution = self.recalibrate(measured_lines, line_list, pixel, order)
        # evaluate wavelenlgth_solution at every point in the trace
        # do some qc checks. Throw warnings when necessary.
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


class IdentifyFeatures(Stage):
    """
    Stage to identify all sharp emission-like features on an Arc lamp frame.
    """
    def do_stage(self, image: EchelleSpectralCCDData):
        # identify emission feature (pixel, order) positions.
        features = identify_features(image.data, image.uncertainty, image.mask, nsigma=5.0, fwhm=EMISSION_FEATURE_SIZE)
        features = group_features_by_trace(features, image.traces)
        features = features[features['id'] != 0]  # throw out features that are outside of any trace.
        # get total flux in each emission feature
        features['flux'] = aperture_extract(features['xcentroid'], features['ycentroid'], image.data,
                                            aperture_width=EMISSION_FEATURE_SIZE, mask=image.mask)
        # blaze correct the emission features fluxes. This speeds up and improves overlap fitting in xwavecal.
        features['corrected_flux'] = features['flux'] / image.blaze['blaze'][features['id'],
                                                                             np.array(features['xcentroid'], dtype=int)]
        image.features = features