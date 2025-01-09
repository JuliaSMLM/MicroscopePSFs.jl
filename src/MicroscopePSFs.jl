# src/MicroscopePSFs.jl

module MicroscopePSFs

using FFTW
using LinearAlgebra
using SpecialFunctions
using Interpolations
using Zygote
using SMLMData
using HDF5
using JSON3

# First include the Zernike module
include("zernike/Zernike.jl")

# Use everything from Zernike within MicroscopePSFs
using .Zernike: ZernikeCoefficients, ZernikeIndexing, OSA, Noll,
    # Basic manipulation
    add_aberration!, reset!, propagate_field, evaluate_pupil, zernikepolynomial, radialpolynomial, max_radial_order,
    # Common aberrations
    add_defocus!, add_spherical!, add_astigmatism!, add_coma!,
    # Coefficient manipulation
    scale!, merge!, rms, trim!,
    # Analysis
    significant_terms,
    # Index conversion exports
    nl2osa, osa2nl, nl2noll, noll2nl, osa2noll, noll2osa, convert_index

# Re-export specific Zernike components
export Zernike
export ZernikeCoefficients, ZernikeIndexing, OSA, Noll

# Common operations
export add_aberration!, reset!, propagate_zernike, evaluate_pupil, zernikepolynomial, radialpolynomial

# Common aberration functions
export add_defocus!, add_spherical!, add_astigmatism!, add_coma!

# Coefficient manipulation and analysis
export scale!, merge!, rms, trim!, significant_terms

# Rest of exports unchanged...
export AbstractPSF
export integrate_pixels, amplitude, integrate_pixels_amplitude
export save_psf, load_psf
export cache_field, load_cached_field

# Pupil Function types
export PupilFunction, VectorPupilFunction

# PSF implementations
export Gaussian2D, Airy2D
export Scalar3DPSF, Scalar3DPupilPSF, Scalar3DZernikePSF
export Vector3DPSF, Vector3DPupilPSF

# Re-export from SMLMData
export AbstractCamera, IdealCamera
export AbstractEmitter, Emitter2D, Emitter3D

# Extended Emitters from SMLMData
export DipoleEmitter3D

# Dipole types and functions
export DipoleVector, DipoleEmitter3D
export calculate_pupil_field, calculate_fresnel_coefficients

# File includes...
include("emitters.jl")
include("pupil.jl")
include("psfs/types.jl")
include("interfaces.jl")
include("utils.jl")
include("integration.jl")

# Individual PSF implementations
include("psfs/gaussian2d.jl")
include("psfs/airy2d.jl")
include("psfs/scalar3d.jl")
include("psfs/vector3d.jl")

end