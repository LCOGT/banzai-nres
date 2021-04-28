Differences from the Commissioning Pipeline
*******************************************

Much of the processing done by BANZAI-NRES is similar to the original pipeline used for commissioning NRES.
We note a few of the differences here.

The main difference is that the commissioning pipeline (available `here <https://github.com/LCOGT/nres-pipe>`_) is
written in IDL, while BANZAI-NRES is primarily written in Python.

Both pipelines perform optimal extraction as per Horne 1986 DOI: 10.1086/131801. The profile is estimated by
fitting polynomials to the quartz lamp flats in the commissioning pipeline, iterating to reject cosmic rays,
while BANZAI-NRES uses the lamp flat pixels directly only normalizing in the spatial direction. Cosmic ray
rejection is an early planned improvement for BANZAI-NRES.

BANZAI-NRES currently only uses the ThAr arc frames taken in during the day to obtain its wavelength solution.
The commissioning pipeline applies an additional shift to the wavelength calibration derived from the simultaneous
ThAr exposure from the calibration fiber.

BANZAI-NRES and the commissioning pipeline estimate the continuum in different ways. The commissioning pipeline
subtracted the an estimate of the blaze function (using a polynomial fit from the stacked quartz lamp flat) and
then applied a high pass filter to the science spectrum before calculating the radial velocity. BANZAI-NRES divides
out an estimate of the blaze function from the stacked quartz lamp flat. Any residual continuum is estimated by
smoothing using a median filter after masking likely absorption lines and is divided out. The continuum normalized
spectrum including both steps in included in the 'normflux' column of the extracted-spectrum data product.

The commissioning pipeline measured the radial velocity using ~12 blocks across each order. We tested this
approach but did not find it improved our radial velocity precision in the current version of BANZAI-NRES.
We plan to reexamine this as we continue to improve the radial velocity precision of BANZAI-NRES.

Currently, the commissioning pipeline does not estimate the stellar parameters of a given target. BANZAI-NRES
estimates the stellar classification by cross-correlating with the PHOENIX stellar atmosphere models from
`Husser et al. 2013, DOI: 10.1051/0004-6361/201219058 <https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract>`_.

BANZAI-NRES propagates an estimate of the formal uncertainties from all of the data processing stages and
include these in the output data products. These are used as weights in the cross correlation function to
measure the radial velocity. We adopt an optimal weighting scheme based on `Zackay & Ofek 2017
DOI: 10.3847/1538-4357/836/2/187 <https://ui.adsabs.harvard.edu/abs/2017ApJ...836..187Z/abstract>`_.

.. include:: idl-differences-data-products.rst
