import numpy as np

from banzai.stages import Stage
from banzai.calibrations import CalibrationStacker, CalibrationUser

from banzai_nres.frames import NRESObservationFrame
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace
from xwavecal.wavelength import wavelength_calibrate, WavelengthSolution
from xwavecal.utils.wavelength_utils import find_nearest
from xwavecal.utils.fiber_utils import lit_fibers
from xwavecal.fibers import IdentifyFibers

from astropy.table import Table
import sep
import logging
import os

logger = logging.getLogger('banzai')

NRES_PRINCIPLE_ORDER_NUMBER = 52
MODELS = {'initial_wavelength_model': {1: [0, 1, 2], 2: [0, 1, 2]},
          'intermediate_wavelength_model': {0: [0, 1, 2], 1: [0, 1, 2], 2: [0, 1, 2]},
          'final_wavelength_model': {0: [0, 1, 2, 3, 4, 5],
                                     1: [0, 1, 2, 3, 4, 5],
                                     2: [0, 1, 2, 3, 4, 5],
                                     3: [0, 1, 2, 3, 4, 5],
                                     4: [0]}}


class ArcStacker(CalibrationStacker):
    def __init__(self, runtime_context):
        super(ArcStacker, self).__init__(runtime_context)

    @property
    def calibration_type(self):
        return 'DOUBLE'


class ArcLoader(CalibrationUser):
    """
    Loads the wavelengths from the nearest Arc-emission lamp (wavelength calibration) in the db.
    If the traces have shifted (e.g. the instrument was serviced), then we do not load any of this information
    (new wavelengths will be created from scratch in WavelengthCalibrate).
    """
    @property
    def calibration_type(self):
        return 'DOUBLE'

    def on_missing_master_calibration(self, image):
        if image.obstype.upper() == 'DOUBLE':
            return image
        else:
            super().on_missing_master_calibration(image)

    def apply_master_calibration(self, image: NRESObservationFrame, master_calibration_image):
        image.wavelength = master_calibration_image.wavelength
        return image


class LineListLoader(CalibrationUser):
    """
    Loads the reference line list for wavelength calibration
    """
    @property
    def calibration_type(self):
        return 'LINELIST'

    def do_stage(self, image):
        master_calibration_file_info = self.get_calibration_file_info(image)
        if master_calibration_file_info is None:
            return self.on_missing_master_calibration(image)
        line_list = np.genfromtxt(master_calibration_file_info['path'])
        logger.info('Applying master calibration', image=image,
                    extra_tags={'master_calibration':  os.path.basename(master_calibration_file_info['path'])})
        return self.apply_master_calibration(image, line_list)

    def apply_master_calibration(self, image, line_list):
        image.line_list = line_list
        return image


class WavelengthCalibrate(Stage):
    """
    Stage to recalibrate wavelength-calibration (e.g. arc lamp) frames.
    We re-wavelength calibrate from scratch if the Arc does not have reference id's assigned to each trace.
    We lightly recalibrate the wavelength calibration if the arc has reference id's assigned to each trace.
    """
    def do_stage(self, image):
        ref_id, fibers = get_ref_ids_and_fibers(image.num_traces)

        features = image.features
        features['order'], features['fiber'] = ref_id[features['id'] - 1], fibers[features['id'] - 1]
        x2d, y2d = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
        for fiber in list(set(fibers)):
            pixel, order = np.arange(image.data.shape[1]), np.sort(ref_id[fibers == fiber])
            this_fiber = features['fiber'] == fiber
            if image.wavelengths is None:
                # blind solve for the wavelengths of the emission lines.
                features['wavelength'][this_fiber], m0 = wavelength_calibrate(features[this_fiber], image.line_list,
                                                                              pixel, order, m0_range=(40, 60))
            else:
                # adopt wavelengths for the emission lines from the existing solution, good to 1 pixel.
                # Note: we could improve this considerably with ndimage.map_coordinates(..., order=1, prefilter=False).
                features['wavelength'][this_fiber] = image.wavelengths[features['order'][this_fiber].astype(int),
                                                                       features['pixel'][this_fiber].astype(int)]
                m0 = None  # TODO get m0 from wavelengths of lines?
            # update the wavelength solution
            wavelength_solution = self.recalibrate(features[this_fiber], image.line_list, pixel, order, m0)

            # overwrite old wavelengths with the new wavelengths
            for trace_id, ref_id in zip([image.spectrum['id'], ref_id]):
                this_trace = np.isclose(image.traces, trace_id)
                image.wavelengths[this_trace] = wavelength_solution(x2d[this_trace],
                                                                    ref_id * np.ones_like(x2d[this_trace]))
        # save the features with their wavelengths
        image.features = features
        return image

    def recalibrate(self, measured_lines, line_list, pixel, order, principle_order_number):
        """
        Recalibrate an existing wavelength solution using new measured spectral features.
        :param measured_lines:
        :param line_list:
        :param pixel:
        :param order:
        :param principle_order_number:
        :return: wavelength_solution. xwavecal.wavelength.WavelengthSolution.
        wavelength_solution(x, order) will give the wavelength of the len(x) points with pixel and order coordinates
        x and order, respectively.
        where x and order are ndarray's of the same shape.
        """
        wavelength_solution = WavelengthSolution(model=MODELS.get('final_wavelength_model'),
                                                 min_order=np.min(order), max_order=np.max(order),
                                                 min_pixel=np.min(pixel), max_pixel=np.max(pixel),
                                                 measured_lines=measured_lines, reference_lines=line_list,
                                                 m0=principle_order_number)
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
    nsigma = 5.0  # minimum signal to noise @ peak flux for a feature to be counted.
    fwhm = 6.0  # minimum feature size for the feature to be counted.

    def do_stage(self, image):
        # identify emission feature (pixel, order) positions.
        features = identify_features(image.data, image.uncertainty, image.mask, nsigma=self.nsigma, fwhm=self.fwhm)
        features = group_features_by_trace(features, image.traces)
        features = features[features['id'] != 0]  # throw out features that are outside of any trace.
        # mask data
        masked_data, masked_err = np.copy(image.data), np.copy(image.uncertainty)
        masked_data[np.isclose(image.mask, 1)], masked_err[np.isclose(image.mask, 1)] = 0, 0
        # get total flux in each emission feature. For now just sum_circle, although we should use sum_ellipse.
        features['flux'], features['fluxerr'], _ = sep.sum_circle(masked_data, features['xcentroid'], features['ycentroid'],
                                                                  self.fwhm, gain=1.0, err=masked_err)
        # blaze correct the emission features fluxes. This speeds up and improves overlap fitting in xwavecal.
        features['corrected_flux'] = features['flux'] / image.blaze['blaze'][features['id'],
                                                                             np.array(features['xcentroid'], dtype=int)]
        image.features = features
        return image


def get_ref_ids_and_fibers(num_traces):
    fibers, ref_id = np.zeros(num_traces), np.zeros(num_traces)
    fibers[1::2] = 1  # group alternating traces as part of the same fiber
    for fiber in [0, 1]:
        ref_id[fibers == fiber] = np.arange(np.count_nonzero(fibers == fiber))
    # note that the fiber designation does not matter, all that matters is that we separate the two AGU's wavelength
    # solutions.
    return ref_id, fibers
