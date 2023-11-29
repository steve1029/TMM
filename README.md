# TMM
<<<<<<< HEAD
A short python script for the Transfer Matrix Method. I made it during my early days in graduate school.

## Installation
Just download whole package and run `$ python example.py`

## How to use
1. First, make an TMM object.
    ```Python
    import TMM
    example = TMM.TMM()
    ```

1. Set the range of wavelength.
    ```python
    wl = np.arange(550,750,0.2) * nm
    example.set_wavelength(wl)
    ```

1. Then, set the incident angle. 
    Only the incident angle is considered in this short script.
    Azimuthal angle is not tunable.

    To set an incident angle,
    ```Python
    example.set_incidentangle(angle=np.pi/6, unit='radian')
    ```

1. Set the material properties of the medium.
    If the medium is consist of 4 layers, put six arguments.
    The first and last arguments are the medium indice for the input(incident) field and output(transmitted) field.
    If the 4 layers are placed in vacuum or air, the first and last argument should be 1.
    The thickness of the mediums can be set by `set_mediumthick` method.
    It also takes 4 arguments because the thickness of the input and output medium is not considered in TMM.

    ```python
    example.set_mediumindex(1, 1.2, 1.1, 2.2, 1.4, 1)
    example.set_mediumtype('nonmagnetic')
    example.set_mediumthick(800*nm, 502.4*nm, 753.6*nm, 402.4*nm)
    ```
=======
A Julia code for TMM, an acronym for the Transfer Matrix Method.
<p align="center">
  <img src="/singleblock_LtoR_tot.png" />
  <em>Field visualization of TE wave, incidence from left-to-right on a single dielectric slab.</em>
</p>

# Features

# Installation

# Examples
>>>>>>> 4db6fff51b810e86051e926959b96276ac03b406
