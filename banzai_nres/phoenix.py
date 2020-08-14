import os
from astropy.io import fits
import boto3
import io
from banzai_nres.dbs import get_phoenix_model_record, get_resource_file


class PhoenixModelLoader:
    def __init__(self, db_address):
        phoenix_wavelength_record = get_resource_file(db_address, 'phoenix_wavelengths')
        self._wavelengths = self.open_model_fits_file(phoenix_wavelength_record.location,
                                                      phoenix_wavelength_record.filename)

    def load(self, db_address, T_effective, log_g, metallicity, alpha):
        # Load in the template
        model_record = get_phoenix_model_record(db_address, T_effective, log_g, metallicity, alpha)
        flux = self.open_model_fits_file(model_record.location, model_record.filename)
        return {'wavelength': self._wavelengths, 'flux': flux}

    @staticmethod
    def open_model_fits_file(file_location, filename):
        if 's3' in file_location:
            s3 = boto3.client('s3')
            buffer = io.BytesIO()
            s3.download_fileobj(file_location.replace('s3://', ''), filename, buffer)
            buffer.seek(0)
            hdu = fits.open(buffer, memmap=False)
        else:
            hdu = fits.open(os.path.join(file_location, filename))

        data = hdu[0].data.copy()
        return data
