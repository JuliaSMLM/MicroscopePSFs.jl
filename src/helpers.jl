## Helper functions 


"""
    makeobs(r)

make a set of observations points in a square window

# Arguments
- `r`   : A range  

returns a vector of Tuples e.g. (x,y)   

# Example:
obs=makeobs(-8:8)
"""
function makeobs(r)

    obspoints = Vector{Tuple}(undef, length(r)^2)
    cnt = 1
    for ii in r
        for jj in r
                obspoints[cnt] = (ii, jj)
                cnt += 1
        end
    end
    return obspoints
end

"""
    save(psffile::String,p::PSF)

save PSF type as .jld2 file

# Arguments
- `psffile` : The full saving path
- `p`       : MicroscopePSFs type
"""
function save(psffile::String,p::PSF)
    JLD2.@save psffile p
    return nothing
end

"""
    load(psffile::String)

load PSF file

# Arguments
- `psffile` : The full path to the psf file

# Return 
The MicroscopePSFs type from the file

"""
function load(psffile::String)
    JLD2.@load psffile p
    return p
end