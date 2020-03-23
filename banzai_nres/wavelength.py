import numpy as np

from banzai.stages import Stage
from banzai_nres.frames import EchelleSpectralCCDData
from xwavecal.wavelength import wavelength_calibrate, SolutionRefineOnce, WavelengthSolution
import logging

logger = logging.getLogger('banzai')


class ArcLoader(Stage):
    def do_stage(self, image: EchelleSpectralCCDData):
        if not traces.blind_solve:
            # load in ref_id column.
            # load in wavelengths


class WavelengthCalibrate(Stage):
    def do_stage(self, image: EchelleSpectralCCDData):
        # image.spectrum = Table({'id': trace_ids, 'flux': flux, 'uncertainty': np.sqrt(variance)})
        # identify emission lines and get (pixel, order) positions.
        # helps to have flat fielded emission line fluxes.
        if ref_id is None:
            # reidentify reference ids.
            # blind solve wavelength
            measured_lines = wavelength_calibrate(measured_lines, line_list, np.arange(4096), np.arange(num_orders),
                                                  principle_order_number=52, wavelength_models=wavelength_models)
            # have wavelength_calibrate add 'wavelength' to measured_lines.
        else:
            measured_lines['wavelength'] = image.wavelengths[(measured_lines['pixel'].astype(int), measured_lines['order'])]

        wavelengths_to_fit = get_nearest(line_list, measured_lines['wavelength'])
        wavelength_solution = WavelengthSolution().solve(measured_lines, wavelengths_to_fit, weights=np.array([1]))
        # in xwavecal: make solve a cls method so that it returns WavelengthSolution
        image.wavelengths = wavelength_solution.wavelength(x, order)
        # make __call__ method.
        # do some qc checks. Throw warnings necessary.
    return image