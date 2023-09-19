# Saving and Importing Data

## Import

The data may be imported using the importpsf function located in src/import.jl. This function requires inputs of the raw filename and the PSF type. 
PSF types supported in this package are "scalar3D", "immPSF", and "splinePSF". You may also choose to input the zposition of the stage, the source 
software, and mvtype as "bead" or "stage."

importpsf returns the PSF type (p), the 3D PSF stack, the pupil function type (h), and the ZernikeCoefficients Type (z).

### Import Example:

    p,PSFstack,h,z = PSF.importpsf(filename,"splinePSF",zstage = 1.0) 

## Save

Saving should be done using the proper format as shown in the example below. You can create the PSF using the Scalar3D function by inputting the na, n,
lambda, and pixel size. The filename should have the extension removed before being put into the save function which will save the PSF as a jld2.

### Save Example

    # Create a scalar PSF
    na=1.2
    n=1.3
    λ=.6 
    pixelsize=.1

    p=PSF.Scalar3D(na,λ,n,pixelsize)

    # Put the filename into the proper format
    psffile = splitext(filename)[1]*".jld2"
    
    # Save it
    PSF.save(psffile,p)
