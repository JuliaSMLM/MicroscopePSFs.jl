

"""
    ZernikeCoefficients 

A mutable struct to hold the Zernike coefficients

# Fields
- `mag`             : Zernike coefficients of pupil mangnitude
- `phase`           : Zernike coefficients of pupil phase

# Constructor
    ZernikeCoefficients(n::Int)
- `n`               : number of Zernike coefficients
"""
mutable struct ZernikeCoefficients 
    mag::Vector{Real}
    phase::Vector{Real}
end

function ZernikeCoefficients(n::Int) 
    mag=zeros(n)
    mag[1]=1.0
    return ZernikeCoefficients(mag,zeros(n)) 
end


"""
    nl2noll(n::Int,l::Int)

convert the `n` and `l` indexes into a Noll linear index

# Arguments
- `n`               : radial index
- `l`               : azimuthal index

# Returns
- `j`               : Noll linear index
"""
function nl2noll(n::Int,l::Int)::Int
    mm = abs(l)
    j = n * (n + 1) / 2 + 1 + max(0, mm - 1);
    if ((l > 0) & (mod(n, 4) >= 2)) | ((l < 0) & (mod(n, 4) <= 1))
       j = j + 1
    end
    return j
end

"""
    noll2nl(j::Int)

convert the Noll index `j` into `n` and `l` 

# Arguments
- `j`               : Noll linear index

# Returns
- `n`               : radial index
- `l`               : azimuthal index
"""
function noll2nl(j::Int)::Tuple{Int,Int}
    n = ceil((-3 + sqrt(1 + 8*j)) / 2);
    l = j - n * (n + 1) / 2 - 1;
    if mod(n, 2) != mod(l, 2)
       l = l + 1;
    end
    if mod(j, 2) == 1
       l= -l;
    end
    return n,l
end

"""
    nl2osa(n::Int,l::Int)

convert the `n` and `l` indexes into a OSE linear index

# Arguments
- `n`               : radial index
- `l`               : azimuthal index

# Returns
OSA linear index
"""
function nl2osa(n::Int,l::Int)::Int
    return (n*(n+2)+l)/2
end

"""
    osa2nl(j::Int)

convert the OSA index `j` into `n` and `l` 

# Arguments
- `j`               : OSA linear index

# Returns
- `n`               : radial index
- `l`               : azimuthal index
"""
function osa2nl(j::Int)::Tuple{Int,Int}
    if j==0
        return 0,0
    end
    n=floor( (-1+sqrt(1+8*j))/2)
    l=2*j-n*(n+2)
    return n,l
end

"""
    noll2osa(j::Int)

convert the Noll index `j` to OSA index `j` 

# Arguments
- `j`               : Noll linear index

# Returns
OSA linear index
"""
function noll2osa(j::Int)::Int
    n,l = noll2nl(j)
    return nl2osa(n,l) 
end

"""
    osa2noll(j::Int)

convert the OSA index `j` to Noll index `j` 

# Arguments
- `j`               : OSA linear index

# Returns
Noll linear index
"""
function osa2noll(j::Int)::Int
    n,l = osa2nl(j)
    return nl2noll(n,l) 
end
"""
    radialpolynomical(n::Int,m::Int,ρ)

return the value of the `n,m` radial polynomial at `ρ`

values of `ρ>1` will return zero

# Arguments
- `n`               : radial index
- `m`               : azimuthal index
- `ρ`               : radial coordinate

# Returns
- `r`               : value of the radial polynomial at `ρ`
"""
function radialpolynomial(n::Int,m::Int,ρ)
    if ρ>1
        return 0
    end

    if m==0
        g = sqrt(n+1)
    else
        g = sqrt(2*n+2)
    end

    r=0.0;
    for k=0:(n-m)/2
        p=ρ.^(n-2*k)
        coef= g*(-1)^k * prod((n - m)/2 - k + 1 : n - k) / 
                     (factorial(Int(k)) * factorial(Int((n + m)/2 - k)));
     r+=coef*p;
    end
    return r
end

"""
    zernike(n::Int,l::Int,ρ,ϕ)

return the value the `n,l` zernike polynomial at `ρ,ϕ`

Note that `l` is in `-n<:2:n`.

# Arguments
- `n`               : radial index
- `l`               : azimuthal index
- `ρ`               : radial coordinate
- `ϕ`               : azimuthal coordinate

# Returns
value of the zernike polynomial at `ρ,ϕ`
"""
function zernikepolynomial(n::Int,l::Int,ρ,ϕ)
    m=abs(l)
    r=radialpolynomial(n,m,ρ)
    if l<0
        return r*sin(m*ϕ)
    else
        return r*cos(m*ϕ)
    end
end

"""
    zernike(j::Int,ρ,ϕ; linearindex="OSA")

return the value of jth zernike polynomial at `ρ,ϕ` for a given Zernike linear index

# Arguments
- `j`               : Zernike linear index
- `ρ`               : radial coordinate
- `ϕ`               : azimuthal coordinate
- `linearindex`     : linear index type, either "OSA" or "Noll"

# Returns
value of the zernike polynomial at `ρ,ϕ`
"""
function zernikepolynomial(j::Int,ρ,ϕ; linearindex="OSA")

    if linearindex=="OSA"
        n,l=osa2nl(j)
        return zernikepolynomial(n,l,ρ,ϕ)  
    end

    if linearindex=="Noll"
        n,l=noll2nl(j)
        return zernikepolynomial(n,l,ρ,ϕ)    
    end

end






