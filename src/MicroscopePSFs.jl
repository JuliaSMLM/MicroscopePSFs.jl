module MicroscopePSFs
using FFTW
using LinearAlgebra
using SpecialFunctions
using Interpolations
using Zygote
using SMLMData
using HDF5
using JSON3

# Core types and interfaces
export AbstractPSF
export integrate_pixels, amplitude, integrate_pixels_amplitude
export save_psf, load_psf
export cache_field, load_cached_field

# PSF implementations
export Gaussian2D, Airy2D, Scalar3D, Vector3D

include("psfs/types.jl")
include("psfs/functions.jl")
include("interfaces.jl")
include("utils.jl")
include("integration.jl")
# include("interpolation.jl")
# include("io.jl")

# Individual PSF implementations
include("psfs/gaussian2d.jl")
include("psfs/airy2d.jl")
# include("psfs/scalar3d.jl")
# include("psfs/vector3d.jl")

end
