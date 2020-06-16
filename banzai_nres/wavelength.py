import numpy as np

from banzai.stages import Stage
from banzai.calibrations import CalibrationStacker, CalibrationUser

from banzai_nres.frames import NRESObservationFrame
from banzai_nres.utils.wavelength_utils import identify_features, group_features_by_trace, get_principle_order_number
from xwavecal.wavelength import find_feature_wavelengths, WavelengthSolution
from xwavecal.utils.wavelength_utils import find_nearest

import sep
import logging
import os
import pkg_resources
from astropy.table import Table


logger = logging.getLogger('banzai')

WAVELENGTH_SOLUTION_MODEL = {0: [0, 1, 2, 3, 4, 5],
                             1: [0, 1, 2, 3, 4, 5],
                             2: [0, 1, 2, 3, 4, 5],
                             3: [0, 1, 2, 3, 4, 5],
                             4: [0]}

# TODO refactor xwavecal so that we dont need this. We only need to set flux_tol to 0.5
OVERLAP_SETTINGS = {'min_num_overlaps': 5, 'overlap_linear_scale_range': (0.5, 2), 'flux_tol': 0.4, 'max_red_overlap': 1000, 'max_blue_overlap': 2000}


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
            super(ArcLoader, self).on_missing_master_calibration(image)

    def apply_master_calibration(self, image: NRESObservationFrame, master_calibration_image):
        image.wavelengths = master_calibration_image.wavelengths
        image.fibers = master_calibration_image.fibers
        image.meta['L1IDARC'] = master_calibration_image.filename, 'ID of ARC/DOUBLE frame'
        return image


class LineListLoader(CalibrationUser):
    """
    Loads the reference line list for wavelength calibration
    """
    LINE_LIST_FILENAME = pkg_resources.resource_filename('banzai_nres', 'data/ThAr_atlas_ESO.txt')
    @property
    def calibration_type(self):
        return 'LINELIST'

    def do_stage(self, image):
        if not os.path.exists(self.LINE_LIST_FILENAME):
            return self.on_missing_master_calibration(image)
        line_list = np.genfromtxt(self.LINE_LIST_FILENAME)[:, 1].flatten()
        logger.info('Applying master calibration', image=image,
                    extra_tags={'master_calibration':  os.path.basename(self.LINE_LIST_FILENAME)})
        return self.apply_master_calibration(image, line_list)

    def apply_master_calibration(self, image, line_list):
        image.line_list = line_list
        return image


class WavelengthCalibrate(Stage):
    """
    Stage to recalibrate wavelength-calibration (e.g. arc lamp) frames.
    We re-wavelength calibrate from scratch if the image does not have a pre-existing wavelength solution.
    We lightly recalibrate the wavelength calibration if the image has a pre-existing wavelength solution.

    NOTE: Pixels that are uncalibrated (e.g. regions of the CCD chip where there is not actually any spectrum present,
    or if the calibration failed on one or both fibers) are assigned a ZERO wavelength. This is because zeros compress
    well in fits fpack. So a pixel with a Zero means that pixel has an UNDEFINED (uncalibrated) wavelength.
    """
    M0_RANGE = (48, 55)  # range of possible values for the integer principle order number.

    def do_stage(self, image):
        image.features = self.init_feature_labels(image.num_traces, image.features)
        if image.wavelengths is None:
            image.features['wavelength'] = np.zeros_like(image.features['pixel'], dtype=float)  # init wavelengths
            for fiber in list(set(image.features['fiber'])):
                logger.info('Blind solving for the wavelengths for fiber {0} (arbitrary label).'.format(fiber), image=image)
                this_fiber = image.features['fiber'] == fiber
                image.features['wavelength'][this_fiber] = find_feature_wavelengths(dict(image.features[this_fiber]), image.line_list,
                                                                                    max_pixel=image.data.shape[1] - 1, min_pixel=0,
                                                                                    m0_range=self.M0_RANGE,
                                                                                    overlap_settings=OVERLAP_SETTINGS)
        else:
            logger.info('Getting feature wavelengths from past solution.')
            # Note: we could improve this estimate by using ndimage.map_coordinates(..., order=1, prefilter=False).
            image.features['wavelength'] = image.wavelengths[image.features['ycentroid'].astype(int),
                                                             image.features['xcentroid'].astype(int)]
        self.refine_wavelengths(image)
        return image

    def refine_wavelengths(self, image):
        ref_ids, fiber_ids, trace_ids = get_ref_ids_and_fibers(image.num_traces)
        image.wavelengths = np.zeros_like(image.data, dtype=float)
        for fiber in list(set(fiber_ids)):
            this_fiber = image.features['fiber'] == fiber
            if np.all(np.isnan(image.features['wavelength'][this_fiber])) or np.all(np.isclose(image.features['wavelength'][this_fiber], 0)):
                logger.error('All zeros for image.wavelengths for fiber {0}. Will not refine '
                             'wavelengths and marking image as bad.'.format(fiber), image=image)
                image.is_bad = True
                continue

            m0 = get_principle_order_number(np.arange(*self.M0_RANGE), image.features[this_fiber])
            logger.info('Principle order number is {0} for fiber {1}'.format(m0, fiber), image=image)
            wavelength_model = self.fit_wavelength_model(image.features[this_fiber], image.line_list, m0)

            # save the wavelengths of each spectral feature
            image.features['wavelength'][this_fiber] = wavelength_model(image.features['pixel'][this_fiber],
                                                                        image.features['order'][this_fiber])
            # populate a wavelengths image using the 2d polynomial function wavelength_solution
            x2d, y2d = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
            for trace_id, ref_id in zip(trace_ids[fiber_ids == fiber], ref_ids[fiber_ids == fiber]):
                this_trace = np.isclose(image.traces, trace_id)
                image.wavelengths[this_trace] = wavelength_model(x2d[this_trace], ref_id * np.ones_like(x2d[this_trace]))
            # Set the true physical order number
            ref_ids[fiber_ids == fiber] += m0
        # Calculate which fiber is which.
        # the first ref_id = 80 (arbitrary, but somewhere about the middle of the detector)
        # is largest fiber number that is lit.
        largest_fiber_index = fiber_ids[ref_ids.tolist().index(80)]
        if image.fiber2_lit:
            largest_fiber_number = 2
            smaller_fiber_number = 1
        else:
            largest_fiber_number = 1
            smaller_fiber_number = 0

        true_fibers = [largest_fiber_number if fiber_id == largest_fiber_index else smaller_fiber_number
                       for fiber_id in fiber_ids]
        # set fibers attribute on image
        image.fibers = Table({'trace': trace_ids, 'fiber': true_fibers})

    @staticmethod
    def init_feature_labels(num_traces, features):
        ref_ids, fiber_ids, trace_ids = get_ref_ids_and_fibers(num_traces)
        features['order'] = ref_ids[features['id'] - 1].astype(int)
        features['fiber'] = fiber_ids[features['id'] - 1].astype(int)
        return features

    @staticmethod
    def fit_wavelength_model(features, line_list, m0):
        """
        Fit a polynomial model to the feature wavelengths with light outlier rejection. This recalibrates
        the wavelength solution.
        :param features:
        :param line_list:
        :param m0: Real diffraction order number of the i=0 diffraction order. Also known as m0.
        :return: wavelength_solution. xwavecal.wavelength.WavelengthSolution.
        wavelength_solution(x, order) will give the wavelength of the len(x) points with pixel and order coordinates
        x and order, respectively.
        where x and order are ndarray's of the same shape.
        """
        # TODO update xwavecal so that WavelengthSolution does not need min_order, max_order, min_pixel, max_pixel.
        wcs = WavelengthSolution(model=WAVELENGTH_SOLUTION_MODEL, measured_lines=dict(features),
                                 min_order=0, max_order=np.max(features['order']),
                                 min_pixel=0, max_pixel=np.max(features['pixel']),
                                 reference_lines=line_list, m0=m0)
        wavelengths_to_fit = find_nearest(features['wavelength'], np.sort(line_list))
        weights = np.ones_like(wavelengths_to_fit, dtype=float)
        # consider weights = features['flux']/features['flux_err'] or 1/features['flux_err']**2
        # reject lines who have residuals with the line list in excess of 0.1 angstroms (e.g. reject outliers)
        weights[~np.isclose(wavelengths_to_fit, wcs.measured_lines['wavelength'], atol=0.1)] = 0
        wcs.model_coefficients = wcs.solve(wcs.measured_lines, wavelengths_to_fit, weights)
        return wcs


class IdentifyFeatures(Stage):
    """
    Stage to identify all sharp emission-like features on an Arc lamp frame.
    """
    nsigma = 15.0  # minimum signal to noise for the feature (sum of flux / sqrt(sum of error^2) within the aperture)
    # where the aperture is a circular aperture of radius fwhm.
    fwhm = 6.0  # fwhm estimate of the elliptical gaussian PSF for each feature
    num_features_max = 4000  # maximum number of features to keep. Will keep the num_features_max highest S/N features.

    def do_stage(self, image):
        # identify emission feature (pixel, order) positions.
        features = identify_features(image.data, image.uncertainty, image.mask, nsigma=self.nsigma/2, fwhm=self.fwhm)
        features = group_features_by_trace(features, image.traces)
        features = features[features['id'] != 0]  # throw out features that are outside of any trace.
        logger.info('{0} emission features found on this image'.format(len(features)), image=image)
        # get total flux in each emission feature. For now just sum_circle, although we should use sum_ellipse.
        features['flux'], features['fluxerr'], _ = sep.sum_circle(image.data, features['xcentroid'], features['ycentroid'],
                                                                  self.fwhm, gain=1.0, err=image.uncertainty, mask=image.mask)
        if image.blaze is not None:
            logger.info('Blaze correcting emission feature fluxes', image=image)
            # blaze correct the emission features fluxes. This speeds up and improves overlap fitting in xwavecal.
            features['corrected_flux'] = features['flux'] / image.blaze['blaze'][features['id'] - 1,
                                                                                 np.array(features['xcentroid'], dtype=int)]

        # cutting which lines to keep:
        # calculate the error in the centroids provided by identify_features()
        features['centroid_err'] = self.fwhm / np.sqrt(features['flux'])
        # Keep features that pass the signal to noise check.
        valid_features = features['flux']/features['fluxerr'] > self.nsigma
        features = features[valid_features]
        # Keep the highest S/N features
        features = features[np.argsort(features['flux']/features['fluxerr'])[::-1]]
        features = features[:int(self.num_features_max)]

        # report statistics
        logger.info('{0} emission features on this image will be used'.format(len(features)), image=image)
        if len(features) == 0:
            logger.error('No emission features found on this image!', image=image)

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
    return ref_id, fibers, np.arange(1, num_traces + 1)
