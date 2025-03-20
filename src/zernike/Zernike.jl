# src/zernike/Zernike.jl

"""
    Zernike

Module for working with Zernike polynomials and PSF pupil functions.
Provides tools for:
- Computing Zernike polynomials with Noll normalization (RMS = 1)
- Converting between Noll and OSA/ANSI indexing schemes
- Manipulating pupil functions via Zernike coefficients
- Representing aberrations in optical systems

Note: This module uses Noll indexing and normalization throughout.
"""
module Zernike

using LinearAlgebra
using SpecialFunctions  # For Bessel functions

# Type exports
export ZernikeCoefficients

# Function exports for polynomial computation
export zernikepolynomial, radialpolynomial, max_radial_order
export evaluate_pupil

# Index conversion exports
export nl2osa, osa2nl, nl2noll, noll2nl, osa2noll, noll2osa

# Utility exports
export rms, significant_terms

# Include all submodule files
include("types.jl")
include("indexing.jl")
include("polynomials.jl")
include("utils.jl")

end # module