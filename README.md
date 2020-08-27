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
Included in banzai_nres/data is the ThAr_atlas_ESO_original_air.txt which was fetched from
http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html on August 27 2020.

These lines were converted to vacuum wavelength using the revised Edlen Equation ( K P Birch and M J Downs 1993 
Metrologia 30 155, https://doi.org/10.1088/0026-1394/30/3/004). The revised Edlen Equation that we used is 
banzai_nres.utils.wavelength_utils.index_of_refraction_Edlen_revised . We converted to vacuum by multiplying all the wavelengths in 
ThAr_atlas_ESO_original_air.txt by the index of refraction at each line's AIR wavelength. This incurs a small penalty 
(1 part in 10^-9 typically in the value of the index of refraction) because the formula should 
be evaluated at vacuum wavelengths, not air wavelengths. The resulting vacuum list, which is used by this pipeline 
to calibrate the spectrographs, is banzai_nres/data/ThAr_atlas_ESO_vacuum.txt

## License

## Support
[API documentation]()  
[Create an issue](https://issues.lco.global/)
