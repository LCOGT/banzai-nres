*******************
Wavelength Solution
*******************

During the day, we take a series of ThAr arc lamp exposures through both sets of fibers to use for
wavelength calibration. These are stored in the `archive <https://archive.lco.global>`_ under OBSTYPE=DOUBLE.

BANZAI-NRES stacks these frames and calibrates the wavelength solution using
`xwavecal <https://github.com/gmbrandt/xwavecal>`_ as described in
`Brandt et al. 2020 DOI: 10.3847/1538-3881/ab929c <https://ui.adsabs.harvard.edu/abs/2020AJ....160...25B/abstract>`_.

Line List
^^^^^^^^^^^^^
NOTE: NRES wavelength calibrations are in *vacuum wavelengths*.

Our line list is the ThAr atlas from ESO which was fetched from
http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html on 27 August 2020. We converted
the air wavelengths to vacuum wavelengths using the original 1966
Edlen Equation (`Murphy et al. 2007, DOI 10.1111/j.1365-2966.2007.11768.x <https://ui.adsabs.harvard.edu/abs/2007MNRAS.378..221M/abstract>`_).

Wavelength calibration metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
BANZAI-NRES calculates metrics characterizing the quality of the wavelength solution. The metrics are only calculated on features that are
matched to lines in the catalog. Features are matched to the closest feature in the line list and then outliers/non-matches
are rejected using sigma clipping. The header keywords that contain the metrics and rule-of-thumb values for
good wavelength calibrations are included in the list below. Note that these header keywords are only included in the stacked
ThAr super-calibration frames.

- 'SIGLAM': Standard deviation of measured wavelengths - line list wavelengths in Angstroms.
  Typically this should be between 0.003 and 0.01 Angstroms for NRES.

- 'RVPRECSN': Standard deviation of measured wavelengths - line list wavelengths converted to velocity units
  via :math:`\frac{\Delta \lambda}{\lambda} = \frac{\Delta v}{c}`. This metric has the advantage that it is
  normalized accross the full wavelength range and is in the same units as the radial velocity. For NRES,
  this is typically between 4 and 10 m/s.


- 'WAVRCHI2': Reduced :math:`\chi^2` of the line list wavelengths - wavelength model of the detected features.
  This should be close to one. Much larger than one would suggest the wavelength solution is not describing the
  true dispersion well. Much less than one would suggest that error bars on the centroids of the detected features
  are over-estimated.

- 'NLINEDET': Number of features detected on the frame. For NRES, this is typically ~4000.

- 'NLINES': Number of matched lines to the line list. This should be ~75% of the features detected on the frame.
  If this is less than 1000, blind solves of the wavelength solution will likely fail.

To check for a failed wavelength solution, a good test is to compare adjacent orders and ensure the
features line up in regions that the wavelengths overlap.
