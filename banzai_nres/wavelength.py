import numpy as np

from banzai.stages import Stage
from banzai.calibrations import CalibrationStacker, CalibrationUser

from banzai_nres.frames import NRESObservationFrame
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace
from xwavecal.wavelength import wavelength_calibrate, WavelengthSolution
from xwavecal.utils.wavelength_utils import find_nearest

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
    M0_RANGE = (40, 60)  # range of possible values for the integer principle order number.

    def do_stage(self, image):
        ref_id, fibers = get_ref_ids_and_fibers(image.num_traces)

        image.features['order'] = ref_id[image.features['id'] - 1]
        image.features['fiber'] = fibers[image.features['id'] - 1]
        image.features['wavelength'] = np.zeros_like(image.features['pixel'], dtype=float)
        blind_solve = True if image.wavelengths is None else False
        for fiber in list(set(fibers)):
            image = self.calibrate_fiber(fiber, fibers, ref_id, image, blind_solve)
        return image

    def calibrate_fiber(self, fiber, fibers, ref_ids, image, blind_solve=True):
        logger.info('Wavelength calibrating the {0} fiber'.format({0: 'first', 1: 'second'}[fiber]), image=image)
        features = image.features

        pixel, order = np.arange(image.data.shape[1]), np.sort(ref_ids[fibers == fiber]).astype(int)
        trace_ids = np.arange(1, image.num_traces + 1)
        this_fiber = np.equal(features['fiber'], fiber)
        x2d, y2d = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))

        if blind_solve:
            logger.info('Blind solving for the wavelengths of this fiber.')
            feature_wavelengths, m0 = wavelength_calibrate(dict(features[this_fiber]), image.line_list,
                                                           pixel, order, m0_range=self.M0_RANGE)
        else:
            # adopt wavelengths for the emission lines from the existing solution, good to 1 pixel.
            # Note: we could improve this considerably with ndimage.map_coordinates(..., order=1, prefilter=False).
            feature_wavelengths = image.wavelengths[features['ycentroid'][this_fiber].astype(int),
                                                    features['xcentroid'][this_fiber].astype(int)]
            # Get m0. This m0 value needs to be exactly the correct integer value.
            center_wavelengths = self.get_center_wavelengths(image.wavelengths, image.traces, trace_ids[fibers == fiber])
            m0 = self.get_principle_order_number(np.arange(*self.M0_RANGE), center_wavelengths, ref_ids[fibers == fiber])
        if feature_wavelengths is not None:
            features['wavelength'][this_fiber] = feature_wavelengths
            image.wavelengths = np.zeros_like(image.data, dtype=float)

            # update the wavelength solution
            wavelength_solution = self.recalibrate(features[this_fiber], image.line_list, pixel, order, m0)
            # overwrite old wavelengths with the new wavelengths
            for trace_id, ref_id in zip(trace_ids[fibers == fiber], ref_ids[fibers == fiber]):
                this_trace = np.isclose(image.traces, trace_id)
                image.wavelengths[this_trace] = wavelength_solution(x2d[this_trace],
                                                                    ref_id * np.ones_like(x2d[this_trace]))
            features['wavelength'][this_fiber] = wavelength_solution(features['pixel'][this_fiber],
                                                                     features['order'][this_fiber])
            image.features = features
        else:
            logger.error('The xwavecal wavelength solution has failed on this fiber. Check xwavecal logging.')
        return image

    @staticmethod
    def recalibrate(features, line_list, pixel, order, principle_order_number):
        """
        Recalibrate an existing wavelength solution using new measured spectral features.
        :param features:
        :param line_list:
        :param pixel:
        :param order:
        :param principle_order_number:
        :return: wavelength_solution. xwavecal.wavelength.WavelengthSolution.
        wavelength_solution(x, order) will give the wavelength of the len(x) points with pixel and order coordinates
        x and order, respectively.
        where x and order are ndarray's of the same shape.
        """
        wcs = WavelengthSolution(model=MODELS.get('final_wavelength_model'),
                                 min_order=np.min(order), max_order=np.max(order),
                                 min_pixel=np.min(pixel), max_pixel=np.max(pixel),
                                 measured_lines=dict(features), reference_lines=line_list,
                                 m0=principle_order_number)
        wavelengths_to_fit = find_nearest(features['wavelength'], np.sort(line_list))
        weights = np.ones_like(wavelengths_to_fit, dtype=float)
        # consider weights = features['flux']/features['flux_err'] or 1/features['flux_err']**2
        # reject lines who have residuals with the line list in excess of 0.1 angstroms (e.g. reject outliers)
        weights[~np.isclose(wavelengths_to_fit, wcs.measured_lines['wavelength'], atol=0.1)] = 0
        wcs.model_coefficients = wcs.solve(wcs.measured_lines, wavelengths_to_fit, weights)
        return wcs

    @staticmethod
    def get_principle_order_number(m0_values, center_wavelengths, ref_ids):
        """
        Finds the principle order number m0. Selects the m0 such that the function y(i) = (m0 + i) * central_wavelengths
        has the smallest slope. I.e. this selects the m0 that allows constant/(m0+i) to best fit central_wavelengths.
        This is exactly what CERES does. See equation 3 of Brahm et al. 2016.

        :param m0_values: ndarray of integers. 1d.
        :param center_wavelengths: ndarray of wavelengths. The wavelengths of the center pixels down the detector, for
        one fiber
        :param ref_ids: ndarray of integers. The reference id's of the traces from which central_wavelengths came.
        Same shape as central_wavelengths. i.e. central_wavelengths[0] is the center wavelength of the trace with reference
        id ref_ids[0]
        :return: m0: int.
        The principle order number for the fiber from which central_wavelengths were taken.
        This is the true order index of the the trace that corresponds to ref_id[0].
        """
        slopes = []
        for m0 in m0_values:
            # note: replacing np.ptp with some outlier resistant measure of the scatter would be more robust.
            slopes.append(np.ptp(center_wavelengths * (m0 + ref_ids)))

        if np.count_nonzero(np.isclose(slopes, np.min(slopes), rtol=0.01)) > 1:
            logger.warning('Two or more viable principle order numbers for this fiber! The m0 recovered from the '
                           'wavelength solution could be wrong!')
        return m0_values[np.argmin(slopes)]

    @staticmethod
    def get_center_wavelengths(wavelengths, traces, trace_ids):
        # get the center wavelength from each trace in trace_ids.
        cc = wavelengths.shape[1] // 2
        center_wavelengths = [wavelengths[:, cc][traces[:, cc] == i].flatten()[0] for i in trace_ids]
        return center_wavelengths


class IdentifyFeatures(Stage):
    """
    Stage to identify all sharp emission-like features on an Arc lamp frame.
    """
    nsigma = 3.0  # minimum signal to noise @ peak flux for a feature to be counted.
    fwhm = 6.0  # minimum feature size in pixels for the feature to be counted.

    def do_stage(self, image):
        # identify emission feature (pixel, order) positions.
        features = identify_features(image.data, image.uncertainty, image.mask, nsigma=self.nsigma, fwhm=self.fwhm)
        features = group_features_by_trace(features, image.traces)
        features = features[features['id'] != 0]  # throw out features that are outside of any trace.
        if len(features) == 0:
            logger.error('No emission features found on this image!', image=image)
        # mask data
        masked_data, masked_err = np.copy(image.data), np.copy(image.uncertainty)
        masked_data[np.isclose(image.mask, 1)], masked_err[np.isclose(image.mask, 1)] = 0, 0
        # get total flux in each emission feature. For now just sum_circle, although we should use sum_ellipse.
        features['flux'], features['fluxerr'], _ = sep.sum_circle(masked_data, features['xcentroid'], features['ycentroid'],
                                                                  self.fwhm, gain=1.0, err=masked_err)
        # blaze correct the emission features fluxes. This speeds up and improves overlap fitting in xwavecal.
        features['corrected_flux'] = features['flux'] / image.blaze['blaze'][features['id'] - 1,
                                                                             np.array(features['xcentroid'], dtype=int)]
        image.features = features
        return image


def get_ref_ids_and_fibers(num_traces):
    # this function always assumes two fibers are lit.
    fibers, ref_id = np.zeros(num_traces), np.zeros(num_traces)
    fibers[1::2] = 1  # group alternating traces as part of the same fiber
    for fiber in [0, 1]:
        ref_id[fibers == fiber] = np.arange(np.count_nonzero(fibers == fiber))
    # note that the fiber designation does not matter, all that matters is that we separate the two AGU's wavelength
    # solutions.
    return ref_id, fibers
