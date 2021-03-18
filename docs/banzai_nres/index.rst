*************************
banzai-nres Documentation
*************************

This is the documentation for banzai-nres.

Reference/API
=============

.. automodapi:: banzai_nres


Telluric Lines
==============
Our telluric spectrum was generated with https://github.com/kgullikson88/Telluric-Fitter .
We exclude regions of the spectrum with strong telluric absorption.


Output Data Products
====================

Here we cover the various quantities contained in the banzai-nres output data products.


Master wavelength calibration files
-----------------------------------

Master wavelength calibrations are stacked arc (double) exposures, wavelength calibrated primarily
with xwavecal routines.


Wavelength calibration metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In each master wavelength calibration file, the primary hdu header holds five wavelength metrics. The keys for
these are:

* 'SIGLAM'
* 'RVPRECSN'
* 'WAVRCHI2'
* 'NLINEDET'
* 'NLINES'

Before we launch into detail in each, we have to define what a 'matched line' is. Matched lines are any detected
features that have a nearby line in the atlas (list of reference lab wavelengths) within 0.1 Angstroms.

'NLINES' is the number of matched lines.

'NLINEDET' is the total number of lines observed on the detector.

The last three quantities are calculated from the distribution of the wavelength residuals. This distribution is
calculated as follows: for each matched line, calculate the difference between its wavelength and that of
its corresponding line in the atlas.

'SIGLAM' is the standard deviation of the distribution of wavelength residuals.

'RVPRECSN' is slightly more complicated. It is the error on the estimate of the mean of the distribution of residuals,
converted
into velocity via delta lambda / lambda = v/c . I.e. it is how well you know zero-point of the wavelength residuals
(in velocity space). This naively sets the maximum precision you can attain on the instrument.

'WAVRCHI2' this is the formal chisquared statistic of the distribution of wavelength residuals, i.e. each residual divided by the
standard error (in wavelength) of the line centroid position (in wavelength).


Verifying that the wavelength solution is good.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a variety of checks to make sure your wavelength solution is good. Most of them involve checking
the numbers from the above section.

A good wavelength solution, for any NRES unit, will have the five metrics in the following ranges:


* | 'SIGLAM' This should be between 0.003 and 0.01 Angstroms. It should be roughly equal to
  | 'PRECISN' / speed of light * sqrt('MLINENUM') * (mean wavelength).
  | Adopting a mean wavelength = 5000 Angstroms is fine.
* | 'PRECISN' This should be between 4 and 10 m/s . High m/s values above 20 might indicate a failure of the wavelength
  | solution, or an extremely low S/N on the frame.
* | 'CHISQ' between 0.2 and 5. An extremely high number here might mean errors are massively underestimated and
  | vice-versa.
* | 'LINENUM' between 3500 and 5000 (i.e. 1750 and 2500 lines per fiber).
* | 'MLINENUM' should be roughly 50-80% of LINENUM, usually 70%. E.g. 3500 for 5000 total lines, or 2450 for 3500 total
  | lines. If 'MLINENUM' is below 1000, your wavelength solution is probably overconstrained in the center of the
  | detector and the edges are missed. This frame is probably too low of signal to noise.


If all the above checks pass, you most likely have a very good wavelength solution. If any fail, you will want to
look at an extracted science spectrum that uses this calibration. A good test is to look
at the wavelength overlap regions (regions where two orders overlap in wavelength). Absorption or emission lines
in any overlap region should line up. The wavelength calibration failed if that is not true.
Compare the frame with past wavelength calibrations that succeeded.


The line list
^^^^^^^^^^^^^
NOTE: NRES wavelength calibrations are in *vacuum wavelengths*.

Our line list is the ThAr atlas from ESO which was fetched from
http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html on August 27 2020. We converted
the air wavelengths to vacuum wavelengths using the original 1966
Edlen Equation (Murphy et al. 2007, DOI 10.1111/j.1365-2966.2007.11768.x).
