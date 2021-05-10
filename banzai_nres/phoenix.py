import os
import io
import numpy as np

from astropy.io import fits
import boto3

from banzai_nres.dbs import get_resource_file, PhoenixModel
from banzai.context import Context
from scipy import interpolate
from banzai_nres.continuum import ContinuumNormalizer
from banzai_nres.utils import phoenix_utils
import logging


logger = logging.getLogger('banzai')

class PhoenixModelLoader:
    def __init__(self, runtime_context: Context):
        self.runtime_context = runtime_context
        phoenix_wavelength_record = get_resource_file(self.runtime_context.db_address, 'phoenix_wavelengths')
        self._wavelengths = self.open_model_fits_file(phoenix_wavelength_record.location,
                                                      phoenix_wavelength_record.filename,
                                                      self.runtime_context)

    def load(self, model_record: PhoenixModel):
        flux = self.open_model_fits_file(model_record.location, 
                                         model_record.filename, 
                                         self.runtime_context)
        return {'wavelength': self._wavelengths, 'flux': flux}

    @staticmethod
    def open_model_fits_file(file_location: str, filename: str, runtime_context: Context):
        if 's3' in file_location:
            s3 = boto3.client('s3', 
                              aws_access_key_id=runtime_context.PHOENIX_MODEL_AWS_ACCESS_KEY_ID, 
                              aws_secret_access_key=runtime_context.PHOENIX_MODEL_AWS_SECRET_ACCESS_KEY)
            buffer = io.BytesIO()
            s3.download_fileobj(file_location.replace('s3://', ''), filename, buffer)
            buffer.seek(0)
            hdu = fits.open(buffer, memmap=False)
        else:
            hdu = fits.open(os.path.join(file_location, filename))

        data = hdu[0].data.copy()
        return data


def normalize_phoenix_model(args):
    model_file, wavelength_filename, output_dir = args
    logger.info(f'Normalizing {model_file}')
    wavelength_hdu = fits.open(wavelength_filename)
    optical = np.logical_and(wavelength_hdu[0].data >= 3000.0, wavelength_hdu[0].data <= 10000.0)
    hdu = fits.open(model_file)
    # take the data to only be the optical region
    hdu[0].data = hdu[0].data[optical]
    wavelength, flux = wavelength_hdu[0].data[optical], hdu[0].data
    # thin the spectrum by a factor of 25 so that this finishes in a few seconds instead of a few minutes.
    continuum = interpolate.interp1d(wavelength[::25],
                                     ContinuumNormalizer.get_continuum_model(flux[::25], wavelength[::25], 881),
                                     bounds_error=False)(wavelength)
    hdu[0].data /= continuum
    # Save the model file to the output directory
    Teff = hdu[0].header['PHXTEFF']
    logg = hdu[0].header['PHXLOGG']
    metallicity = hdu[0].header['PHXM_H']
    alpha = hdu[0].header['PHXALPHA']
    output_filename = phoenix_utils.parameters_to_filename(Teff, logg, metallicity, alpha)
    hdu.writeto(os.path.join(args.output_dir, output_filename), overwrite=True)

