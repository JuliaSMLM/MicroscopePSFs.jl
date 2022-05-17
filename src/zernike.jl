

"""
    ZernikeCoefficients 

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
    nm2noll(n::Int,l::Int)

    convert the `n` and `l` indexes into a Noll linear index
"""
function nl2noll(n::Int,l::Int)
    mm = abs(l)
    j = n * (n + 1) / 2 + 1 + max(0, mm - 1);
    if (l > 0 & mod(n, 4) >= 2) | (l < 0 & mod(n, 4) <= 1)
       j = j + 1
    end
    return j
end

"""
    noll2nm(j::Int)

    convert the Noll index `j` into `n` and `l` 
"""
function noll2nl(j::Int)
    n = ceil((-3 + sqrt(1 + 8*j)) / 2);
    l = j - n * (n + 1) / 2 - 1;
    if mod(n, 2) != mod(m, 2)
       l = l + 1;
    end
    if mod(j, 2) == 1
       l= -l;
    end
    return n,l
end

"""
    nm2osa(n::Int,m::Int)

    convert the `n` and `m` indexes into a OSE linear index
"""
function nl2osa(n::Int,l::Int)
    return (n*(n+2)+l)/2
end

"""
osa2nm(j::Int)

    convert the OSA index `j` into `n` and `l` 
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
    radialpolynomical(n::Int,m::Int,ρ)

    return the value of the `n,m` radial polynomial at `ρ`

    values of `ρ>1` will return zero

"""
function radialpolynomial(n::Int,m::Int,ρ)
    if ρ>1
        return 0
    end

    r=0.0;
    for k=0:(n-m)/2
        p=ρ.^(n-2*k)
        coef= (-1)^k * prod((n - m)/2 - k + 1 : n - k) / 
                     (factorial(Int(k)) * factorial(Int((n + m)/2 - k)));
     r+=coef*p;
    end
    return r
end

"""
    zernike(n::Int,l::Int,ρ,ϕ)

    return the value the `n,l` zernike polynomial at `ρ,ϕ`

    Note that l is in -n<:2:n
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






