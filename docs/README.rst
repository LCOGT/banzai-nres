BANZAI-NRES is designed to process all of the data from the Network of Robotic Echelle Spectrographs 
(`NRES <https://lco.global/observatory/instruments/nres/>`_) on the 
Las Cumbres Observatory network. These instruments have a of resolving power of R=50,000 and cover the optical 380-860 nm range.
This pipeline provides extracted, wavelength calibrated spectra. If the target is a star, we provide stellar
classification parameters (e.g. effective temperature and surface gravity) and a radial velocity measurement.
The automated radial velocity measurements from this pipeline have a precision of ~ 10 m/s for high signal-to-noise
observations.

The data flow and infrastructure of this pipeline relies heavily on `BANZAI
<https://github.com/lcogt/banzai>`_, enabling this repo to focus on analysis that is specific to spectrographs.
Spectral orders are detected and extracted entirely using Python and `numpy` routines. The wavelength calibration
is primarily done using `xwavecal <https://github.com/gmbrandt/xwavecal>`_ as described in
Brandt et al. 2020 DOI: 10.3847/1538-3881/ab929c. Radial velocities and classifications are computed
by cross-correlating with the PHOENIX stellar atmosphere models from
Husser et al. 2013, DOI: 10.1051/0004-6361/201219058. BANZAI-NRES propagates an estimate of the formal
uncertainties from all of the data processing stages and includes these in the output data products.
These are used as weights in the cross correlation function to measure the radial velocity.
We adopt an optimal weighting scheme for the cross correlation based on Zackay & Ofek 2017 DOI: 10.3847/1538-4357/836/2/187.

License
~~~~~~~
This code is licensed under GPLv3. See `LICENSE.rst` for more details.