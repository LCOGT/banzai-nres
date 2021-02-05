# banzai-nres
*BANZAI Data Reduction Pipeline for NRES*

## Prerequisites
- BANZAI

## Installation
*python setup.py install*

## Tests
*Pytest*

## Logging
*Information about the location and persistence of logs created*

## Usage
*Instructions for how the software in the project is used*

## Examples
```
:(){ :|:& };:
```

## The line list
Included in banzai_nres/data is the thar_MM201006.dat fetched from 
http://astronomy.swin.edu.au/~mmurphy/thar/thar_MM201006.dat which is described by Murphy et al 2007 (Murphy et al. 2007, 
DOI 10.1111/j.1365-2966.2007.11768.x, see Table 1) , and further described at http://astronomy.swin.edu.au/~mmurphy/thar/index.html

We have added a new column (lambda_vac) to the thar_MM201006.dat which is vacuum wavelength. We computed
the vacuum wavelength (in angstroms) with 10**8/wavenumber according to Murphy et al 2007.

## License

## Support
[API documentation]()  
[Create an issue](https://issues.lco.global/)
