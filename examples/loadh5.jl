

using MicroscopePSFs
PSF=MicroscopePSFs
using Plots
using JSON

using HDF5
filename = raw"examples\psfmodel_test_insitu_single.h5"
f = h5open(filename,"r")


zernike_coeff = read(f["res/zernike_coeff"])
params = JSON.parse(attrs(f)["params"])


mag=zernike_coeff[:,1]
phase=zernike_coeff[:,2]


z=PSF.ZernikeCoefficients(mag,phase)

# Create a scalar PSF
na=params["option_params"]["NA"]
n=params["option_params"]["RI_imm"]
λ=params["option_params"]["emission_wavelength"]
pixelsize=.1

p=PSF.Scalar3D(na,λ,n,pixelsize;z=z)

# Or Create a 2D Gaussian approximation
#p=PSF.Gauss2D(p)

# calculate the PSF at a point
PSF.pdf(p,(0,0,0),(0.0,0.0,0.0))

# calculate the PSF in a region
sz=32
roi=[(x,y,0) for y=-sz/2:(sz/2-1), 
    x=-sz/2:(sz/2-1)] 

# out of focus by up to 0.5 physical unit (usually micron)
zrange=cat(dims=1,collect(LinRange(-0.5,0.5,5))
            ,collect(LinRange(0.5,-0.5,5)))
anim = @animate for z ∈ zrange
   heatmap(PSF.pdf(p,roi,(0.0,0.0,z)), aspectratio=:equal, 
   yflip = true, colorbar=:none)
end
gif(anim, fps = 1)


