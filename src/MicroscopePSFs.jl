# src/MicroscopePSFs.jl

"""
    MicroscopePSFs

A Julia package for working with microscope Point Spread Functions (PSFs).
This package provides implementations of common PSF models and tools for 
integrating them with camera geometry for single-molecule localization 
microscopy applications.

# API Overview
For a comprehensive overview of the API, use the help mode on `api`:

    ?api

Or access the complete API documentation programmatically:

    docs = MicroscopePSFs.api()
"""
module MicroscopePSFs

using LinearAlgebra
using SpecialFunctions
using SMLMData
using HDF5
using Interpolations
using Dates 

# First include the Zernike module
include("zernike/Zernike.jl")

# Use everything from Zernike within MicroscopePSFs
using .Zernike: ZernikeCoefficients,
    # Basic manipulation
    evaluate_pupil, zernikepolynomial, radialpolynomial, max_radial_order,
     
    # Analysis
    rms,
    significant_terms,
    
    # Index conversion exports
    nl2osa, osa2nl, nl2noll, noll2nl, osa2noll, noll2osa

# Re-export specific Zernike components
export Zernike
export ZernikeCoefficients

# Common operations
export zernikepolynomial, radialpolynomial

# PSF exports
export AbstractPSF
export integrate_pixels, amplitude, integrate_pixels_amplitude
export save_psf, load_psf 

# Pupil Function types
export PupilFunction, VectorPupilFunction

# PSF implementations
export GaussianPSF, AiryPSF
export ScalarPSF
export VectorPSF
export SplinePSF

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
include("pupils/pupil.jl")
include("pupils/vector_tools.jl")
include("pupils/vector_pupil.jl")
include("psfs/types.jl")
include("interfaces.jl")
include("utils.jl")

include("integration/integration_core.jl")
include("integration/integration_single.jl")
include("integration/integration_multi.jl")

# Individual PSF implementations
include("psfs/gaussian2d.jl")
include("psfs/airy2d.jl")
include("psfs/scalar3d.jl")
include("psfs/vector3d.jl")
include("psfs/spline_psf.jl")

# I/O functions
include("io/io.jl")

# Include the API overview functionality
include("api.jl")

end