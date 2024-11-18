# LaserFields

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeist.github.io/LaserFields.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeist.github.io/LaserFields.jl/dev/)
[![Build Status](https://github.com/jfeist/LaserFields.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfeist/LaserFields.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jfeist/LaserFields.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jfeist/LaserFields.jl)

`LaserFields.jl` is a library to describe the time-dependent electric fields of
a laser pulse and can be installed with `Pkg.add("LaserFields")`. It implements
the same pulse shapes and most of the features of the [laserfields
library](https://github.com/jfeist/laserfields) written in Fortran (and as its
[Python variant](https://github.com/jfeist/pylaserfields)). The "main" interface
is given by the constructor `LaserField(; kwargs...)`, which accepts (as keyword
arguments) the same parameters as the Fortran library reads from parameter
files. Please see the documentation of the Fortran library for the parameter
meanings, conventions used, etc.. `LaserField(; kwargs...)` returns an instance
of a subtype of the abstract base type `LaserField` depending on the parameters.
E.g., to create a Gaussian pulse with a duration (defined as the FWHM of the
intensity) of 6 fs, a wavelength of 800 nm, a peak intensity of 1e14 W/cm^2, and
with the peak at time t=7fs, one should call
```julia
lf = LaserField(form="gaussianI", is_vecpot=true, lambda_nm=800,
                intensity_Wcm2=1e16, duration_as=6000, peak_time_as=7000)
```

Given a `LaserField` instance `lf`, the functions `E_field(lf,t)`,
`E_fourier(lf,ω)`, `A_field(lf,t)`, and `A_fourier(lf,ω)` can be used to obtain,
respectively, the electric field as a function of time, its Fourier transform
(implemented for most pulse shapes), the vector potential as a function of time,
and its Fourier transform. Calling the instance as a function, `lf(t)` returns
the electric field, i.e., is equivalent to `E_field(lf,t)`. The notebooks in the
`examples` folder show some ways to use the library, including how to define a
set of fields through a YAML configuration file.

In addition to the pulses described by each `LaserField` instance, the library
also implements a `LaserFieldCollection` type that represents the sum over
several individual fields. It is also a `LaserField` instance and supports much
of the same interface. Note that some of the parameters it returns are just
"best-effort" values and may not be fully meaningful for the combined field -
e.g., for the carrier frequency `lf.ω0`, it returns the highest value in the
collection, to support use cases where this is used to define maximum time step
in a numerical propagation, or the maximum frequency evaluated in a Fourier
transform.

The "effective" duration of the pulse for n-photon processes can be obtained as
`Teff(lf,n_photon)`, which is the integral over the pulse intensity envelope to
the n-th power (i.e., electric field envelope to the (2n)th power) over the
pulse, see, e.g., https://doi.org/10.1103/PhysRevA.77.043420 (Eq. 14).
