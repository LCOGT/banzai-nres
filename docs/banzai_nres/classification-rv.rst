******************************************************
Stellar Classification and Radial Velocity Measurement
******************************************************

BANZAI-NRES uses cross-correlation methodology in order to measure stellar parameters and radial velocities.
In all cases we cross-correlate the normalized spectra with a Phoenix model spectrum
(`Husser et al. 2013, DOI: 10.1051/0004-6361/201219058 <https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract>`_).

Stellar Classification
~~~~~~~~~~~~~~~~~~~~~~

BANZAI-NRES obtains an initial guess for the stellar effective temperature by querying the Gaia catalog
(`Gaia Collaboration et al. 2016, DOI: 10.1051/0004-6361/201629272 <https://ui.adsabs.harvard.edu/abs/2016A%26A...595A...1G/abstract>`_;
`Gaia Collaboration et al. 2021, arXiv: <https://ui.adsabs.harvard.edu/abs/2020arXiv201201533G/abstract>`_)
at the coordinates of the observation. It then cross-correlates the spectra with models from a Phoenix model
grid in steps of :math:`T_{eff}`, :math:`\log g`, [Fe/H], and [:math:`\alpha`/Fe] around the initial guess.
Accounting for the effects of stellar rotation is a potential future development.
The set of parameters resulting in the highest cross-correlation peak is taken as the stellar parameter
estimate and the corresponding Phoenix model is used as the cross-correlation template to generate the RV measurement.

The stellar classifications are stored in the BANZAI-NRES database. If the star is observed again, it will not be re-classified and the same
template will be used for all RV measurements.

The following header keywords store the stellar classification:

- 'TEFF': Stellar effective temperature (K).

- 'LOGG': Stellar surface gravity (cgs units).

- 'FEH': Stellar metallicity [Fe/H] (dex).

- 'ALPHA': Stellar alpha abundance [alpha/Fe] (dex).

- 'CLASSIFY': Equals 1 if this spectrum was classified, or 0 if the classification was taken from an previous spectrum of this target.

Radial Velocity Measurement
~~~~~~~~~~~~~~~~~~~~~~~~~~~

BANZAI-NRES measures the radial velocity of the target by cross-correlating the best Phoenix template found by the stellar classification step
with the normalized spectra. This cross-correlation is performed on an order-by-order basis over 5140-6276 Angstroms, a range that 
has high signal-to-noise but is not affected by telluric absorption. It first computes the CCFs on a coarse grid with a wide velocity range,
selects the highest peak, and then uses a fine velocity grid around that peak. The peak of each fine CCF is taken as the per-order velocity.
The RV measurement is the sigma-clipped mean of the per-order velocities.
We compute the barycentric correction using Astropy, and correct the RVs to the barycentric frame per
`Wright & Eastman (2014), DOI: 10.1086/678541 <https://ui.adsabs.harvard.edu/abs/2014PASP..126..838W/abstract>`_.

The following header keywords contain the RV value and associated parameters:

- 'RV': Radial velocity measurement in barycentric frame (m / s)

- 'RVERR': Radial velocity formal uncertainty (m / s)

- 'BARYCORR': Barycentric velocity correction (m / s)

- 'TCORR': Mid-exposure time barycentric Julian date (BJD_TDB)


Telluric Lines
~~~~~~~~~~~~~~
Our telluric spectrum was generated with https://github.com/kgullikson88/Telluric-Fitter .
We exclude regions of the spectrum with strong telluric absorption.
