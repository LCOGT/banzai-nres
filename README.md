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
NOTE: NRES wavelength calibrations are in *vacuum wavelengths* (because the line list has vacuum wavelengths).

Included in banzai_nres/data is the ThAr_atlas_ESO_original_air.txt which was fetched from
http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/tharatlas.html on August 27 2020. The line list
that we use is banzai_nres/data/ThAr_atlas_ESO_vacuum.txt, which is the same list but converted to vacuum wavelengths. 
We generated ThAr_atlas_ESO_vacuum.txt as follows:

We converted every line from ThAr_atlas_ESO_original_air.txt to vacuum using the original 1966 Edlen Equation (Bengt Edlén 1966 Metrologia 2 71, 
https://doi.org/10.1088/0026-1394/2/2/002) for the index of refraction of air. 
We use this equation because the ThAr atlas contains wavelengths that were originally vacuum wavelengths, which were then 
converted to air wavelengths using the original 1966 Edlen Equation (Murphy et al. 2007, 
DOI 10.1111/j.1365-2966.2007.11768.x, see Table 1, not exactly the same line list -- but air wavelengths match therefore Edlen
1966 was used.). See the note below for more details.
 
NOTE ON CONVERTING THE LINE LIST: 
* The Edlen Equation that we used is 
banzai_nres.utils.wavelength_utils.index_of_refraction_Edlen . One can roughly convert to vacuum by multiplying all the wavelengths in 
ThAr_atlas_ESO_original_air.txt by the index of refraction at each line's AIR wavelength. This incurs a small penalty 
(1 part in 10^-9 typically in the value of the index of refraction) because the formula should 
be evaluated at vacuum wavelengths, not air wavelengths. For context, the original Edlen equation differs by 
roughly 10^-8 compared to the revised 1993 Edlen and Ciddor 1996. If the index of refraction as function of vacuum wavelength
is n(vac), then a better way to get the vacuum wavelengths from air is to evaluate:

        wavelength_vacuum = wavelength_air * n(wavelength_air * n(wavelength_air))   (Eq. 1)

    Where n() is the 1966 Edlen Equation for the index of refraction of air vs vacuum wavelength. 
    We use Eq. 1 to convert from the ThAr_atlas_ESO_original_air.txt original air
    wavelengths to vacuum wavelengths. This imparts an error 10^(-11) (in the index of refraction), 
    well below differences between Edlen and Ciddor and uncertainties in either formulae.

## The flag mask
1: bad pixel
2: ..
4: ..
8 : rejected because of potential spectrum edge effects
16: rejected because this wavelength coincides with a telluric line.


## License

## Support
[API documentation]()  
[Create an issue](https://issues.lco.global/)
