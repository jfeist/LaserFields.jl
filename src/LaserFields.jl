module LaserFields

using SpecialFunctions
using DelimitedFiles
using DataInterpolations

export LaserField, LaserFieldCollection, make_laserfield
export E_field, A_field, E_fourier, A_fourier
export start_time, end_time, envelope, Teff

include("constants.jl")
include("typedef.jl")
include("fielddefs.jl")
include("make_field.jl")
include("precompile.jl")

end