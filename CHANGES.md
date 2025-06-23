1.2.2 (2025-06-23)
------------------
- Typo fix in numpy requirement for building

1.2.1 (2025-06-19)
------------------
- Minor fix to have the docker image repo point to ghcr instead of our internal docker registry

1.2.0 (2025-06-08)
------------------
- Migrate to use poetry for dependency management
- Use updated banzai from upstream

1.1.5 (2025-03-11)
------------------
- Fixes to documentation for updated BANZAI installs

1.1.4 (2024-11-21)
-----------------
- Fix to helm chart to define large task queue as the same as normal task queue

1.1.3 (2023-11-17)
------------------
- Fixes to use BANZAI LoggingAdapter

1.1.2 (2023-07-25)
------------------
- Update NRES Frame Factory to handle UNKNOWN as a valid empty coordinate

1.1.0 (2023-04-04)
------------------
- Added calibration frame comparison to reject them if they deviate too much from previous stacks

1.0.17 (2022-08-10)
-------------------
- Bugfixes to not try to classify non-sidereal targets

1.0.17 (2022-07-09)
-------------------
- Update to use OpenSearch now that Elasticsearch is preferred
- Update to make NRES standards reduced data public

1.0.16 (2022-02-22)
-------------------
- Update to use upstream BANZAI version with new version of OCS Ingester library (3.0.4)

1.0.15 (2022-02-22)
-------------------
- Update to use upstream BANZAI version with new version of OCS Ingester library (3.0.3)

1.0.14 (2022-02-15)
-------------------
- Update upstream BANZAI install to use new version of OCS Ingester library (3.0.1)
- Helm chart updates for ingester library

1.0.13 (2022-02-14)
-------------------
- Update upstream BANZAI install to use new version of OCS Ingester library

1.0.12 (2022-01-31)
-------------------
- Bugfix to catch if simbad returns none instead of an empty list.

1.0.11 (2022-01-27)
-------------------
- Bugfix to catch if simbad returns none instead of an empty list.

1.0.10 (2021-11-09)
------------------
- Fix auto-deploy in Jenkinsfile

1.0.9 (2021-09-30)
------------------
- Bugfix for an incorrect unit conversion in RA from simbad responses during classification.

1.0.8 (2021-09-22)
------------------
- Add catch for masked values from Gaia

1.0.7 (2021-09-09)
------------------
- Bugfix to automatic deployment

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
