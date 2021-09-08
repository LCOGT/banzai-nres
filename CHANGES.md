1.0.6 (2021-09-08)
------------------
- Added catches for NaNs in the telescope coordinates
- Added a catch for SIMBAD failures
- Automated deployment

1.0.5 (2021-08-30)
------------------
- More robust determination of which fiber was active

1.0.4 (2021-06-22)
------------------
- Small bugfixes for bright stars
- Minor updates to PDF formatting

1.0.3 (2021-06-10)
------------------
- Bugfix to SNR calculation in PDF generation code

1.0.2 (2021-06-01)
------------------
- Changed SNR calculation to be per resolution element rather than per pixel
- Fixed the extraction so that slightly mismatched trace regions and wavelength 
files do not yield incorrect wavelengths. This resulted from ~1 pix noise
in the center of the trace (in the super lampflats).

1.0.1 (2021-05-24)
------------------
- Add ability to configure task queue name (BANZAI 1.3.2)

1.0.0 (2021-05-13)
------------------
- Initial release

0.5.1 (2020-05-07)
------------------
- Extracted spectra now have wavelengths attached.
- The line list is stored as package data instead of in an sqlite database.

0.4.0 (2019-02-25)
------------------
- Blind master tracing on lampflats
- Fixes to account for upstream changes to BANZAI

0.3.0 (2018-11-13)
------------------
- Added flat stacking to stages.

0.2.0 (2018-10-30)
------------------
- Updated how the console entry points are called based on upstream updates.

0.1.0 (2018-07-13)
------------------
- Initial setup of the package infrastructure
