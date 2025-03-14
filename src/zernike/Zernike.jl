# src/zernike/Zernike.jl

"""
    Zernike

Module for working with Zernike polynomials and PSF pupil functions.
Provides tools for:
- Converting between different Zernike indexing schemes
- Computing Zernike polynomials
- Manipulating pupil functions via Zernike coefficients
- Fitting pupil functions to Zernike bases
- Propagating Zernike-based fields to image space
"""
module Zernike

using LinearAlgebra
using SpecialFunctions  # For Bessel functions

# Type exports
export ZernikeIndexing, OSA, Noll
export ZernikeCoefficients

# Function exports for polynomial computation
export zernikepolynomial, radialpolynomial, max_radial_order
export evaluate_pupil

# Index conversion exports
export nl2osa, osa2nl, nl2noll, noll2nl, osa2noll, noll2osa
export convert_index

# Utility exports
export add_aberration!, reset!
export fit_pupil_to_zernike

# Field propagation exports
export propagate_zernike

# Include all submodule files
include("types.jl")
include("indexing.jl")
include("polynomials.jl")
include("fitting.jl")
include("utils.jl")
include("debye.jl")

end # module