from astropy.io import fits
import numpy as np
import os
import argparse


def format_plus_minus(value, fmt):
    output_str = f'{value:{fmt}}'
    return output_str.replace('+', 'p').replace('-', 'm')

def main():
    parser = argparse.ArgumentParser('This copies out the optical region 300nm - 1000nm into fresh files for BANZAI')
    parser.add_argument('--input-dir', dest='input_dir', help='Top level directory with the PHOENIX models')
    parser.add_argument("--output-dir", dest='output_dir', help='Directory to save the reformatted models')
    args = parser.parse_args()
    model_files = []
    # Traverse the directory to find all of the phoenix files
    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            # Get the wavelength file
            if 'wave' in file.lower():
                wavelength_filename = os.path.join(root, file)
            elif file.endswith(".fits"):
                model_files.append(os.path.join(root, file))

    # Find the indices that correspond to the optical region
    wavelength_hdu = fits.open(wavelength_filename)
    optical = np.logical_and(wavelength_hdu[0].data >= 3000.0, wavelength_hdu[0].data <= 10000.0)

    for model_file in model_files:
        hdu = fits.open(model_file)
        # take the data to only be the optical region
        hdu[0].data = hdu[0].data[optical]
        # Save the model file to the output directory
        Teff = hdu[0].header['PHXTEFF']
        logg = hdu[0].header['PHXLOGG']
        metallicity = hdu[0].header['PHXM_H']
        alpha = hdu[0].header['PHXALPHA']
        output_filename = f'phoenix-{format_plus_minus(int(Teff), "05d")}-{format_plus_minus(logg, "+0.1f")}-{format_plus_minus(metallicity, "+0.1f")}-{format_plus_minus(alpha, "+0.1f")}.fits'
        hdu.writeto(os.path.join(args.output_dir, output_filename))

    # Truncate and save the wavelength file
    wavelength_hdu[0].data = wavelength_hdu[0].data[optical]
    wavelength_hdu.writeto(os.path.join(args.output_dir, 'phoenix_wavelength.fits'))
