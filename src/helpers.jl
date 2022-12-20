## Helper functions 


"""
    makeobs(r)

make a set of observations points in a square windoe  

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
    save(psffile,p)

save PSF type as .jld2 file

#Arguments
- `psffile` : The full saving path
- `p`       : MicroscopePSFs type

#Example:
save('psffile.jld2',p)
"""
function save(psffile::String,p::PSF)
    JLD2.@save psffile p
    return nothing
end

"""
    load(psffile)

load PSF file

#Arguments
- `psffile` : The full path to the psf file

return the MicroscopePSFs type from the file

#Example:
p=load('psffile.jld2')
"""
function load(psffile::String)
    JLD2.@load psffile p
    return p
end