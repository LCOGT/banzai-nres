Data Product Differences from the Commissioning Pipeline
--------------------------------------------------------
The main difference between the data products produced by the commissioning pipeline and BANZAI-NRES is
that we have formatted the data into a FITS Binary table instead of using multiple FITS extensions.
This makes the association between the different pieces of data more explicit and alleviates some issues
when compressing the files with fpack.

Below are the extension names from the FITS extensions from the previous data products and how they map
to the columns in the extracted spectra produced by BANZAI-NRES. Rather than splitting the ThAr fiber
and the science fiber into separate extensions, BANZAI-NRES interleaves the orders in the FITS binary table.


- 'SPECRAW' Extension: 'flux' column

- 'SPECFLAT' Extension: 'normflux' column is the closest analog. The 'normflux' column includes extra
  continuum normalization compared to the commissioning pipeline.

- 'SPECBLAZE' Extension: This is the blaze subtracted spectrum. We do not provide this explicitly, but
  BANZAI-NRES does provide the 'blaze' column and the 'flux' columns to calculate this.

- 'THARRAW' Extension: This is the flux of the ThAr fiber orders are interleaved with the science fiber
  in the 'flux' column .

- 'THARFLAT' Extension: Flat-field divided ThAr fiber spectrum. These are interleaved with the science
  fiber in the 'normflux' column.

- 'WAVESPEC' Extension: Wavelengths solutions per pixel for the science fiber.
  These are stored in the 'wavelength' column in BANZAI-NRES data products. Note that the commissioning
  pipeline stored these in values in nanometers. BANZAI-NRES stores wavelengths in Angstroms.

- 'WAVETHAR' Extension: Wavelengths solutions per pixel for the calibration ThAr fiber.  These are stored in
  the 'wavelength' column in BANZAI-NRES data products and are interleaved with the science spectrum orders.
  Note that the commissioning pipeline stored these in values in nanometers. BANZAI-NRES stores
  wavelengths in Angstroms.

- 'SPECXCOR' Extension: Summed cross correlation function (CCF) for a range of velocities. The velocity for a
  given pixel is saved as the WCS in the header. BANZAI-NRES stores its cross-correlation results in the
  'CCF' extension as a FITS binary table. Rather than saving the summed CCF, we store the CCF per order and
  corresponding velocity coordinates for clarity for the user.

- 'RVBLOCKFIT' Extension: Per-order/per-block RVs in a FITS Binary table. As we do not calculate radial
  velocities using blocks, BANZAI-NRES does not include this extension.
