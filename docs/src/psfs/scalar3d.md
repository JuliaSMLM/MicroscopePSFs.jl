# Scalar3DPSF

The `Scalar3DPSF` model implements a three-dimensional point spread function based on scalar diffraction theory. This model accounts for defocus and optical aberrations using a complex pupil function approach, providing a good balance between physical accuracy and computational efficiency. It is well-suited for 3D imaging applications where polarization effects aren't critical.

## Mathematical Model

The Scalar3DPSF uses the Fourier optics approach to calculate the complex field distribution:

```math
U(\mathbf{r}) = \int_{pupil} P(\boldsymbol{\rho}) e^{i k \boldsymbol{\rho} \cdot \mathbf{r}} d\boldsymbol{\rho}
```

where:
- ``U(\mathbf{r})`` is the complex field amplitude at position ``\mathbf{r} = (x, y, z)``
- ``P(\boldsymbol{\rho})`` is the complex pupil function at pupil coordinates ``\boldsymbol{\rho}``
- ``k = 2\pi / \lambda`` is the wave number
- The intensity is calculated as ``I(\mathbf{r}) = |U(\mathbf{r})|^2``

The pupil function can incorporate various aberrations, typically represented using Zernike polynomials.

## Constructor Options

```julia
Scalar3DPSF(na::Real, wavelength::Real, n::Real; 
            pupil::Union{Nothing, PupilFunction}=nothing,
            pupil_data::Union{Nothing, AbstractMatrix}=nothing,
            coeffs::Union{Nothing, ZernikeCoefficients}=nothing)
```

### Required Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns
- `n`: Refractive index of the medium

### Optional Parameters

- `pupil`: Pre-created `PupilFunction` instance
- `pupil_data`: Complex matrix to initialize the pupil function
- `coeffs`: `ZernikeCoefficients` instance for representing aberrations

You should provide exactly one of the optional parameters. If none are provided, an unaberrated pupil is created.

### Type Parameters
- `T`: Numeric precision type, automatically determined from input

## Basic Usage

### Creating a PSF

```julia
# Create a basic unaberrated 3D PSF
psf = Scalar3DPSF(1.4, 0.532, 1.518)  # NA=1.4, λ=532nm, n=1.518

# Create a PSF with spherical aberration
zc = ZernikeCoefficients(15)  # Up to 15th Zernike polynomial
add_spherical!(zc, 0.5)       # Add 0.5 waves of spherical aberration
psf_aberrated = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)

# Create a PSF with a pre-computed pupil function
pupil = PupilFunction(1.4, 0.532, 1.518, zernike_coeffs)
psf_from_pupil = Scalar3DPSF(1.4, 0.532, 1.518, pupil=pupil)
```

### Direct Evaluation

```julia
# Evaluate PSF at a specific 3D position
x = 0.1  # μm
y = 0.2  # μm
z = 0.5  # μm (distance from focal plane)
intensity = psf(x, y, z)

# Get complex amplitude
amp = amplitude(psf, x, y, z)
```

### Creating PSF Images

```julia
# Create a grid of positions for a 2D slice at a specific z
x = range(-1, 1, length=101)  # μm
y = range(-1, 1, length=101)  # μm
z = 0.0  # focal plane

# Compute PSF values at each position
intensity_2d = [psf(xi, yi, z) for yi in y, xi in x]

# Visualize with CairoMakie
using CairoMakie
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="Scalar3DPSF (z=0μm)",
          xlabel="x (μm)", ylabel="y (μm)")
hm = heatmap!(ax, x, y, intensity_2d, colormap=:viridis)
Colorbar(fig[1, 2], hm)

# Create a 3D PSF stack
z_stack = range(-2, 2, length=5)  # z positions in microns
intensity_3d = Array{Float64}(undef, length(y), length(x), length(z_stack))
for (k, zi) in enumerate(z_stack)
    for (j, yi) in enumerate(y)
        for (i, xi) in enumerate(x)
            intensity_3d[j, i, k] = psf(xi, yi, zi)
        end
    end
end

# Visualize z-stack slices
fig2 = Figure(size=(800, 200))
for (i, zi) in enumerate(z_stack)
    ax = Axis(fig2[1, i], aspect=DataAspect(),
              title="z = $(zi)μm",
              xlabel=i==1 ? "x (μm)" : "",
              ylabel=i==1 ? "y (μm)" : "")
    heatmap!(ax, x, y, intensity_3d[:,:,i], colormap=:viridis)
end
fig2
```

## Integration with Camera

```julia
# Create camera with 100nm pixels (20×20 pixel grid)
pixel_size = 0.1  # μm
camera = IdealCamera(1:20, 1:20, pixel_size)

# Create 3D emitter at position (1μm, 1μm, 0.5μm) with 1000 photons
emitter = Emitter3D(1.0, 1.0, 0.5, 1000.0)  # x, y, z, photons

# Integrate PSF over pixels with 2×2 subsampling
pixels = integrate_pixels(psf, camera, emitter, sampling=2)

# Visualize the camera image
using CairoMakie
fig = Figure(size=(500, 400))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="Integrated Camera Image (z=0.5μm)",
          xlabel="x (μm)", ylabel="y (μm)")
ax.yreversed = true  # Flip y-axis to match camera convention

# Get physical coordinates of pixel centers
x_centers = (1:20) * pixel_size - pixel_size/2
y_centers = (1:20) * pixel_size - pixel_size/2

hm = heatmap!(ax, x_centers, y_centers, pixels', colormap=:viridis)
scatter!(ax, [emitter.x], [emitter.y], color=:red, marker=:cross, markersize=15)
Colorbar(fig[1, 2], hm)
fig
```

## Aberration Modeling

One of the key features of the Scalar3DPSF model is its ability to incorporate optical aberrations through Zernike polynomials:

```julia
# Create Zernike coefficients object with up to 15th term
zc = ZernikeCoefficients(15)

# Add common aberrations
add_defocus!(zc, 1.0)         # 1 wave of defocus
add_astigmatism!(zc, 0.5)     # 0.5 waves of astigmatism
add_coma!(zc, 0.3)            # 0.3 waves of coma
add_spherical!(zc, 0.2)       # 0.2 waves of spherical aberration

# Create PSF with these aberrations
psf = Scalar3DPSF(1.4, 0.532, 1.518, coeffs=zc)

# Compare aberrated and unaberrated PSFs
psf_unaberrated = Scalar3DPSF(1.4, 0.532, 1.518)

# Create comparison images at different z positions
z_positions = [-1.0, 0.0, 1.0]  # μm
x = y = range(-1, 1, length=101)  # μm

fig = Figure(size=(900, 600))
for (i, z) in enumerate(z_positions)
    # Unaberrated PSF
    ax1 = Axis(fig[1, i], aspect=DataAspect(),
               title="Unaberrated (z=$(z)μm)",
               xlabel="x (μm)", ylabel=i==1 ? "y (μm)" : "")
    img1 = [psf_unaberrated(xi, yi, z) for yi in y, xi in x]
    heatmap!(ax1, x, y, img1, colormap=:viridis)
    
    # Aberrated PSF
    ax2 = Axis(fig[2, i], aspect=DataAspect(),
               title="Aberrated (z=$(z)μm)",
               xlabel="x (μm)", ylabel=i==1 ? "y (μm)" : "")
    img2 = [psf(xi, yi, z) for yi in y, xi in x]
    heatmap!(ax2, x, y, img2, colormap=:viridis)
end
fig
```

## Performance Considerations

The Scalar3DPSF model has the following performance characteristics:

- More computationally intensive than 2D models (Gaussian2D, Airy2D)
- Significantly faster than the full vectorial model (Vector3DPSF)
- Performance scales with the resolution of the pupil function
- Computing many PSF values at different positions is more efficient after pre-computing the pupil function
- For performance-critical applications, consider creating a `SplinePSF` from a Scalar3DPSF

## Limitations

The Scalar3DPSF model has several limitations:

1. **No Polarization Effects**: Doesn't account for polarization, which becomes significant at high NA (>1.2)
2. **Simplified Refractive Index Interfaces**: Simplified treatment of refractive index interfaces
3. **Scalar Approximation**: Uses scalar diffraction theory instead of full vector theory
4. **Moderate NA Assumption**: Most accurate for moderate NA values and regions near the focal plane

## Relationship to Other PSFs

The Scalar3DPSF is related to other PSF models in the following ways:

1. **Airy2D**: The in-focus (z=0) slice of an unaberrated Scalar3DPSF is equivalent to the Airy pattern
2. **Vector3DPSF**: The Scalar3DPSF is a simplified version of Vector3DPSF that doesn't account for polarization effects
3. **SplinePSF**: For performance-critical applications, a Scalar3DPSF can be converted to a SplinePSF

For applications requiring higher physical accuracy with high-NA objectives, consider using the `Vector3DPSF` model.