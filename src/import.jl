
"""
    importpsf(filename, psftype; zstage=0.0, source="python", mvtype="bead")

import PSF data from PSF learning software

# Arguments
- `filename`   : file name of the PSF data
- `psftype`    : PSF types that are supported by julia package, options are: "scalar3D", "immPSF", "splinePSF"
- `source`     : software that generates the PSF data, default is "python"
- `zstage`     : position of the sample stage, equal to zero at the coverslip, positive when imaging inside the sample
- `mvtype`     : for immPSF only, options are: "bead", "stage" 


# returns
- `p`           : PSF Type
- `PSFstack`    : a 3D stack of learned PSF
- `z`           : ZernikeCoefficients Type
- `h`           : PupilFunction Type

# Example:
p, PSFstack, z, h = importpsf(filename,psftype)
"""
function importpsf(filename, psftype; zstage=0.0, source="python", mvtype="bead")
    if source == "python"
        f = h5open(filename, "r")
        PSFstack = read(f["res/I_model"])
        #normf = sum(PSFstack,dims=(1,2))
        #PSFstack ./= normf
        params = JSON.parse(attrs(f)["params"])
        pixelsize_x = params["pixel_size"]["x"]
        pixelsize_z = params["pixel_size"]["z"]
        na = params["option"]["imaging"]["NA"]
        RI = params["option"]["imaging"]["RI"]
        #n=collect(values(params["option"]["imaging"]["RI"]))
        n = [RI["med"], RI["cov"], RI["imm"]]
        λ = params["option"]["imaging"]["emission_wavelength"]
        p = []
        z = []
        h = []

        if haskey(f, "res/zernike_coeff")
            zernike_coeff = read(f["res/zernike_coeff"])
            N = size(zernike_coeff)[1]
            j_osa = Array(0:N-1)
            j_noll = osa2noll.(j_osa)
            mag = zernike_coeff[j_noll, 1]
            phase = zernike_coeff[j_noll, 2]
            z = ZernikeCoefficients(mag, phase)
        end

        if haskey(f, "res/pupil")
            pupilcomplex = read(f["res/pupil"])
            ksize = size(pupilcomplex, 1)
            kpixelsize = 2 * na / λ / ksize
            pupil = Float64.(cat(abs.(pupilcomplex), angle.(pupilcomplex), dims=3))
            if psftype == "scalar3D"
                p = Scalar3D(na, λ, n[3], pixelsize_x; inputpupil=pupil, ksize=ksize)
                h = PupilFunction(na, λ, n[3], pixelsize_x, kpixelsize, pupil)
            end

            if psftype == "immPSF"

                p = ImmPSF(na, λ, n, pixelsize_x; inputpupil=pupil, zstage=zstage, ksize=ksize,mvtype=mvtype)    
                h = PupilFunction(na, λ, n[1], pixelsize_x, kpixelsize, pupil)
            end
        end

        if psftype == "splinePSF"
            p = SplinePSF(PSFstack; pixelsize_z=pixelsize_z, pixelsize=pixelsize_x)
        end
    end

    return p, PSFstack, z, h
end

