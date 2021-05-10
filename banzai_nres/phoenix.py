import os
import io

from astropy.io import fits
import boto3

from banzai_nres.dbs import get_resource_file, PhoenixModel
from banzai.context import Context


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
