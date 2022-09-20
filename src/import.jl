
"""
    importpsf(r)

load PSF data from PSF learning software

#Arguments
- `filename`   : file name of the PSF data, file format is .h5

#returns
- `p`           : Scalar3D PSF type
- `PSFstack`    : a 3D stack of learned PSF
- `pixelsize_x` : pixel size at sample plane in x, micron
- `pixelsize_z` : step size in z, micron

#Example:
p, PSFstack, pixelsize_x, pixelsize_z = importpsf(filename)
"""
function importpsf(filename)

    f = h5open(filename,"r")
    PSFstack = read(f["res/I_model"])
    params = JSON.parse(attrs(f)["params"])
    pixelsize_x=params["pixelsize_x"]
    pixelsize_z=params["pixelsize_z"]
    if haskey(f,"res/zernike_coeff")
        zernike_coeff = read(f["res/zernike_coeff"])    
        N = size(zernike_coeff)[1]
        j_osa = Array(0:N-1) 
        j_noll = osa2noll.(j_osa)        
        mag=zernike_coeff[j_noll,1]
        phase=zernike_coeff[j_noll,2]
     
        z=ZernikeCoefficients(mag,phase)

        # Create a scalar PSF
        na=params["option_params"]["NA"]
        n=params["option_params"]["RI_imm"]
        λ=params["option_params"]["emission_wavelength"]
        p=Scalar3D(na,λ,n,pixelsize_x;z=z)
    end
    
    return p, PSFstack, pixelsize_x,pixelsize_z
end

