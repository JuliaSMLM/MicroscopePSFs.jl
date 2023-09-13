# Saving and Loading Data

## Load

The data may be loaded using the importpsf function located in src/import.jl. This function requires inputs of the raw filename and the PSF type. 
PSF types supported in this package are "scalar3D", "immPSF", and "splinePSF". You may also choose to input the zposition of the stage, the source 
software, and mvtype as "bead" or "stage."

### Load Example:

    filename = raw"Y:\Projects\Super Critical Angle Localization Microscopy\ZStack_TIRF_01-10-23\ZStack_psfmodel_zernike_vector_single.h5"
    p,PSFstack,z = PSF.importpsf(filename,"splinePSF",zstage = 1.0) 

## Save

Saving should be done using the proper format as shown in the example below.

### Save Example

    filename = raw"Y:\Projects\Super Critical Angle Localization Microscopy\ZStack_TIRF_01-10-23\ZStack_psfmodel_zernike_vector_single.h5"
    psffile = splitext(filename)[1]*".jld2"

    PSF.save(psffile,p)
