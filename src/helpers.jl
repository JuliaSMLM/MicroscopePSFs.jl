## Helper functions 


"""
    makeobs(r)

make a set of observations points in a square windoe  

#Arguments
- `r`   : A range  

returns a vector of Tuples e.g. (x,y)   

#Example:
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

