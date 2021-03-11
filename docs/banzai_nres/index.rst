*************************
banzai-nres Documentation
*************************

This is the documentation for banzai-nres.

Reference/API
=============

.. automodapi:: banzai_nres



Running Code Tests Locally
==========================

To test code style: "tox -e codestyle"
To test documentation: "tox -e build_docs"
To test the functionality of the code itself (on python 3.8): "tox -e py38-test-alldeps"

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
* 'PRECISN'
* 'WAVECHI2'
* 'NLINES'
* 'NLINESMC'

Before we launch into detail in each, we have to define what a 'matched line' is. Matched lines are any detected
features that have a corresponding line in the atlas within 0.1 Angstroms.

Now some more detail on each.

'NLINESMC' is the number of matched lines.

'NLINES' is the total number of lines observed on the detector.

The last three quantities are calculated from a single distribution, after applying various transformations to it.
This distribution is calculated as follows: for the wavelength
of each feature observed on the detector,
find the nearest line in the laboratory line list. Calculate that difference. Record those differences. Now restrict
the resulting distribution only to matched lines. This is the
distribution, which we will refer to as the distribution of residuals. Again, note that this distribution is restricted
to matched lines only.


'SIGLAM' is the standard deviation of the distribution of wavelength residuals.

'PRECISN' is slightly more complicated. It is the error on the estimate of the mean of the distribution of residuals,
converted
into velocity via delta lambda / lambda = v/c . I.e. it is how well you know zero-point of the wavelength residuals
(in velocity space). This naively sets the maximum precision you can attain on the instrument.

'WAVECHI2' this is the formal chisquared statistic of the distribution of wavelength residuals, i.e. each residual divided by the
standard error (in wavelength) of the line centroid position (in wavelength).


What is a good wavelength solution?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a variety of checks to make sure your wavelength solution is good. Most of them involve simply checking
the numbers from the above section.

A good wavelength solution for any NRES unit will have those quantities in the following ranges:


* | 'SIGLAM' This should be between 0.003 and 0.01 Angstroms. It should be roughly equal to
  | 'PRECISN' / speed of light * sqrt('MLINENUM') * (mean wavelength). mean wavelength = 5000 Angstroms is fine.
* | 'PRECISN' This should be between 4 and 10 m/s . High m/s values above 20 might indicate a failure of the wavelength
  | solution, or extremely low S/N on the frame.
* | 'CHISQ' between 0.2 and 5. An extremely high number here might mean errors are massively underestimated and
  | vice-versa.
* | 'LINENUM' between 3500 and 5000 (i.e. 1750 and 2500 lines per fiber). The code as is limits this value to 2500
  | features per fiber.
* | 'MLINENUM' should be roughly 50-80% of LINENUM, usually 70%. E.g. 3500 for 5000 total lines, or 2450 for 3500 total
  | lines. If 'MLINENUM' is below 1000, your wavelength solution is probably overconstrained in the center of the
  | detector and the edges are missed. This frame is probably too low of signal to noise.


If all the above checks pass, you most likely have a very good wavelength solution. If any fail, you may want to dig
into the frame deeper and look at an extracted science spectrum that uses this calibration. The end all test is to look
at the wavelength overlap regions. You want the extracted spectral lines between two orders to agree almost exactly in
the common wavelength overlap regions between those two orders. If that is not true, the wavelength calibration failed
without a doubt.


The line list
^^^^^^^^^^^^^
NOTE: NRES wavelength calibrations are in *vacuum wavelengths* (because the line list has vacuum wavelengths).

Included in banzai_nres/data is the ThAr_atlas_ESO_original_air.txt which was fetched from
http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html on August 27 2020. The line list
that we use is banzai_nres/data/ThAr_atlas_ESO_vacuum.txt, which is the same list but converted to vacuum wavelengths.
We generated ThAr_atlas_ESO_vacuum.txt as follows:

We converted every line from ThAr_atlas_ESO_original_air.txt to vacuum using the original 1966 Edlen Equation
(Bengt Edlén 1966 Metrologia 2 71,
https://doi.org/10.1088/0026-1394/2/2/002) for the index of refraction of air.
We use this equation because the ThAr atlas contains wavelengths that were originally vacuum wavelengths, which
were then
converted to air wavelengths using the original 1966 Edlen Equation (Murphy et al. 2007,
DOI 10.1111/j.1365-2966.2007.11768.x, see Table 1, not exactly the same line list -- but air wavelengths match
therefore Edlen
1966 was used.). See the note below for more details.

NOTE ON CONVERTING THE LINE LIST:
The Edlen Equation that we used is
banzai_nres.utils.wavelength_utils.index_of_refraction_Edlen . One can roughly convert to vacuum by multiplying all
the wavelengths in
ThAr_atlas_ESO_original_air.txt by the index of refraction at each line's AIR wavelength. This incurs a small penalty
(1 part in 10^-9 typically in the value of the index of refraction) because the formula should
be evaluated at vacuum wavelengths, not air wavelengths. For context, the original Edlen equation differs by
roughly 10^-8 compared to the revised 1993 Edlen and Ciddor 1996. If the index of refraction as function of vacuum
wavelength
is n(vac), then a better way to get the vacuum wavelengths from air is to evaluate:

        wavelength_vacuum = wavelength_air * n(wavelength_air * n(wavelength_air))   (Eq. 1)

Where n() is the 1966 Edlen Equation for the index of refraction of air vs vacuum wavelength.
We use Eq. 1 to convert from the ThAr_atlas_ESO_original_air.txt original air
wavelengths to vacuum wavelengths. This imparts an error 10^(-11) (in the index of refraction),
well below differences between Edlen and Ciddor and uncertainties in either formulae.
