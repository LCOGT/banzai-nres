from banzai.stages import CalibrationMaker, MasterCalibrationDoesNotExist
from astropy.io import fits
import os


def get_trace_coefficients(image, pipeline_context):
    coefficients_and_indices, fiber_order = None, None

    master_trace_filename = CalibrationMaker(pipeline_context).get_calibration_filename(image)
    master_trace_file_path = os.path.join(image.pipeline_context.processed_path, master_trace_filename)
    if image.header['OBSTYPE'] != 'TRACE' and os.path.isfile(master_trace_file_path):
        fiber_order = fits.getheader(master_trace_file_path).get('FIBRORDR')
        coefficients_and_indices = fits.getdata(master_trace_file_path)

        assert coefficients_and_indices is not None
        assert fiber_order is not None

    if image.header['OBSTYPE'] != 'LAMPFLAT' and not os.path.isfile(master_trace_file_path):
        raise MasterCalibrationDoesNotExist

    return coefficients_and_indices, fiber_order
