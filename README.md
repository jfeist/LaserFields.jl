# LaserFields

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeist.github.io/LaserFields.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeist.github.io/LaserFields.jl/dev/)
[![Build Status](https://github.com/jfeist/LaserFields.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeist/LaserFields.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jfeist/LaserFields.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jfeist/LaserFields.jl)

`LaserFields.jl` is a library to describe the time-dependent electric fields of
a laser pulse. It implements the same pulse shapes and most of the features of
the [laserfields library](https://github.com/jfeist/laserfields) written in
Fortran. Please see the documentation of that library for the parameter
meanings, conventions used, etc. In particular, the "main" function
`make_laser_field` accepts the same parameters as the Fortran library parameter
files as keyword arguments. E.g., to create a Gaussian pulse with a duration
(defined as the FWHM of the intensity) of 6 fs, a wavelength of 800 nm, a
peak intensity of 1e14 W/cm^2, and with the peak at time t=7fs, one should call
```julia
lf = make_laser_field(form="gaussianI", is_vecpot=true, lambda_nm=800,
                      intensity_Wcm2=1e16, duration_as=6000, peak_time_as=7000)
```

The "main" interface is provided by the functions `E_field(lf,t)`,
`E_fourier(lf,ω)`, `A_field(lf,t)`, and `A_fourier(lf,ω)`, which give,
respectively, the electric field as a function of time, its Fourier transform
(implemented for most pulse shapes), the vector potential as a function of time,
and its Fourier transform. The notebook in the `examples` folder shows how to
use the library.
