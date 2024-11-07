"""
Provides conversion functions between different Zernike polynomial indexing schemes:
- OSA/ANSI: Based on radial order n and azimuthal frequency l
- Noll: Single index ordering that follows specific pattern
- (n,l): Raw radial and azimuthal indices

All functions include input validation and maintain proper mathematical constraints.
"""

#= Conversion between (n,l) and single-index schemes =#

"""
    nl2osa(n::Integer, l::Integer) -> Int

Convert from (n,l) indices to OSA/ANSI single index.

# Arguments
- `n`: Radial order (≥ 0)
- `l`: Azimuthal frequency, must satisfy:
  * |l| ≤ n
  * n - |l| must be even

# Returns
- OSA/ANSI single index j = (n(n+2) + l)/2

# Examples
```julia
julia> nl2osa(4, 2)
15
```
"""
function nl2osa(n::Integer, l::Integer)::Int
    # Input validation
    n ≥ 0 || throw(ArgumentError("Radial order n must be non-negative"))
    abs(l) ≤ n || throw(ArgumentError("Azimuthal frequency |l| must be ≤ n"))
    mod(n - abs(l), 2) == 0 || throw(ArgumentError("n - |l| must be even"))
    
    return (n * (n + 2) + l) ÷ 2
end

"""
    nl2noll(n::Integer, l::Integer) -> Int

Convert from (n,l) indices to Noll single index.

# Arguments
- `n`: Radial order (≥ 0)
- `l`: Azimuthal frequency, must satisfy:
  * |l| ≤ n
  * n - |l| must be even

# Returns
- Noll single index
"""
function nl2noll(n::Integer, l::Integer)::Int
    # Input validation
    n ≥ 0 || throw(ArgumentError("Radial order n must be non-negative"))
    abs(l) ≤ n || throw(ArgumentError("Azimuthal frequency |l| must be ≤ n"))
    mod(n - abs(l), 2) == 0 || throw(ArgumentError("n - |l| must be even"))
    
    m = abs(l)
    j = n * (n + 1) ÷ 2 + 1 + max(0, m - 1)
    
    if ((l > 0) & (mod(n, 4) >= 2)) | ((l < 0) & (mod(n, 4) <= 1))
        j += 1
    end
    
    return j
end

"""
    osa2nl(j::Integer) -> Tuple{Int,Int}

Convert from OSA/ANSI single index to (n,l) indices.

# Arguments
- `j`: OSA/ANSI index (≥ 0)

# Returns
- Tuple of (n,l) indices
"""
function osa2nl(j::Integer)::Tuple{Int,Int}
    j ≥ 0 || throw(ArgumentError("OSA index must be non-negative"))
    
    n = floor(Int, (-1 + sqrt(1 + 8j)) / 2)
    l = 2j - n * (n + 2)
    
    # Validate result
    mod(n - abs(l), 2) == 0 || throw(ArgumentError("Invalid OSA index"))
    abs(l) ≤ n || throw(ArgumentError("Invalid OSA index"))
    
    return (n, l)
end

"""
    noll2nl(j::Integer) -> Tuple{Int,Int}

Convert from Noll single index to (n,l) indices.

# Arguments
- `j`: Noll index (> 0)

# Returns
- Tuple of (n,l) indices
"""
function noll2nl(j::Integer)::Tuple{Int,Int}
    j > 0 || throw(ArgumentError("Noll index must be positive"))
    
    n = floor(Int, (-3 + sqrt(1 + 8j)) / 2)
    l = 2 * (j - (n * (n + 1) ÷ 2)) - n
    
    # Adjust sign based on Noll convention
    if mod(j, 2) == 0
        l = abs(l)
    else
        l = -abs(l)
    end
    
    # Validate result
    mod(n - abs(l), 2) == 0 || throw(ArgumentError("Invalid Noll index"))
    abs(l) ≤ n || throw(ArgumentError("Invalid Noll index"))
    
    return (n, l)
end

#= Cross-indexing conversions =#

"""
    noll2osa(j::Integer) -> Int

Convert from Noll to OSA/ANSI index.
"""
function noll2osa(j::Integer)::Int
    n, l = noll2nl(j)
    return nl2osa(n, l)
end

"""
    osa2noll(j::Integer) -> Int

Convert from OSA/ANSI to Noll index.
"""
function osa2noll(j::Integer)::Int
    n, l = osa2nl(j)
    return nl2noll(n, l)
end

"""
    convert_index(j::Integer, from::ZernikeIndexing, to::ZernikeIndexing) -> Int

Convert between indexing schemes specified by ZernikeIndexing enum.

# Examples
```julia
julia> convert_index(4, OSA, Noll)
6
julia> convert_index(6, Noll, OSA)
4
```
"""
function convert_index(j::Integer, from::ZernikeIndexing, to::ZernikeIndexing)
    from == to && return j
    if from == OSA && to == Noll
        return osa2noll(j)
    else
        return noll2osa(j)
    end
end

# Convenience method to get n,l from either indexing scheme
"""
    get_nl(j::Integer, indexing::ZernikeIndexing=OSA) -> Tuple{Int,Int}

Get (n,l) indices from either OSA or Noll index.
"""
function get_nl(j::Integer, indexing::ZernikeIndexing=OSA)::Tuple{Int,Int}
    if indexing == OSA
        return osa2nl(j)
    else
        return noll2nl(j)
    end
end
