# Vector3DPSF

The `Vector3DPSF` model implements a comprehensive three-dimensional point spread function based on vectorial diffraction theory. It accounts for polarization effects, high-NA phenomena, dipole emission patterns, and arbitrary aberrations, providing the highest level of physical accuracy among the PSF models. This model is particularly important for high-NA objectives and applications where polarization effects cannot be ignored.

## Mathematical Model

The Vector3DPSF uses the Richards-Wolf vector diffraction theory to calculate the electromagnetic field near focus:

```math
\mathbf{E}(\mathbf{r}) = \int_{pupil} \mathbf{P}(\theta, \phi) e^{i k (\mathbf{s} \cdot \mathbf{r})} \sin\theta \, d\theta \, d\phi
```

where:
- ``\mathbf{E}(\mathbf{r})`` is the electric field vector at position ``\mathbf{r} = (x, y, z)``
- ``\mathbf{P}(\theta, \phi)`` is the complex vector pupil function
- ``\mathbf{s}`` is the unit vector in the direction of wave propagation
- ``\theta`` and ``\phi`` are the polar and azimuthal angles in the pupil
- ``k = 2\pi / \lambda`` is the wave number
- The intensity is calculated as ``I(\mathbf{r}) = |\mathbf{E}(\mathbf{r})|^2``

The model accounts for polarization, refractive index interfaces, and emitter dipole orientation.

## Constructor Options

```julia
Vector3DPSF(na::Real, wavelength::Real, dipole::DipoleVector; 
            n_medium::Real=1.33,
            n_immersion::Real=1.518, 
            n_coverslip::Real=1.518,
            focal_z::Real=0.0, 
            base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
            grid_size::Integer=128)
```

### Required Parameters

- `na`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns
- `dipole`: Dipole orientation vector specified using `DipoleVector`

### Optional Parameters

- `n_medium`: Refractive index of the sample medium (default: 1.33, water)
- `n_immersion`: Refractive index of the immersion medium (default: 1.518, oil)
- `n_coverslip`: Refractive index of the coverslip (default: 1.518, glass)
- `focal_z`: Focal plane position in microns (default: 0.0)
- `base_zernike`: Zernike coefficients for aberrations (default: none)
- `grid_size`: Size of grid for pupil function (default: 128)

### Type Parameters
- `T`: Numeric precision type, automatically determined from input

## Basic Usage

### Creating a PSF

```julia
# Create dipole vectors for different orientations
dipole_x = DipoleVector(1.0, 0.0, 0.0)  # X-oriented dipole
dipole_y = DipoleVector(0.0, 1.0, 0.0)  # Y-oriented dipole
dipole_z = DipoleVector(0.0, 0.0, 1.0)  # Z-oriented dipole
dipole_xy = DipoleVector(0.707, 0.707, 0.0)  # 45° in XY plane

# Create a basic Vector3DPSF with Z-oriented dipole
psf = Vector3DPSF(
    1.4,                # Numerical aperture
    0.532,              # Wavelength in microns
    dipole_z,           # Z-oriented dipole
    n_medium=1.33,      # Sample is water
    n_immersion=1.518   # Immersion oil
)

# Create a PSF with aberrations
zc = ZernikeCoefficients(15)
add_spherical!(zc, 0.5)  # Add 0.5 waves of spherical aberration
psf_aberrated = Vector3DPSF(
    1.4,
    0.532,
    dipole_z,
    n_medium=1.33,
    base_zernike=zc
)
```

### Direct Evaluation

```julia
# Evaluate PSF at a specific 3D position
x = 0.1  # μm
y = 0.2  # μm
z = 0.5  # μm (distance from focal plane)
intensity = psf(x, y, z)

# Get complex vector field amplitude (returns [Ex, Ey])
E_field = amplitude(psf, x, y, z)
Ex = E_field[1]  # X-component of electric field
Ey = E_field[2]  # Y-component of electric field

# Calculate intensity manually if needed
intensity_manual = abs2(Ex) + abs2(Ey)
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
          title="Vector3DPSF (z=0μm, z-oriented dipole)",
          xlabel="x (μm)", ylabel="y (μm)")
hm = heatmap!(ax, x, y, intensity_2d, colormap=:viridis)
Colorbar(fig[1, 2], hm)

# Compare different dipole orientations
psf_x = Vector3DPSF(1.4, 0.532, dipole_x, n_medium=1.33)
psf_y = Vector3DPSF(1.4, 0.532, dipole_y, n_medium=1.33)
psf_z = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33)

fig2 = Figure(size=(900, 300))
psfs = [psf_x, psf_y, psf_z]
titles = ["X-oriented dipole", "Y-oriented dipole", "Z-oriented dipole"]

for (i, (p, title)) in enumerate(zip(psfs, titles))
    ax = Axis(fig2[1, i], aspect=DataAspect(),
              title=title,
              xlabel="x (μm)", ylabel=i==1 ? "y (μm)" : "")
    img = [p(xi, yi, 0.0) for yi in y, xi in x]
    heatmap!(ax, x, y, img, colormap=:viridis)
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

# For dipole emitters, you can use the DipoleEmitter3D type
dipole_emitter = DipoleEmitter3D(1.0, 1.0, 0.5, 1000.0, 0.0, 0.0, 1.0)  # x, y, z, photons, dx, dy, dz

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

The Vector3DPSF supports the same aberration modeling capabilities as Scalar3DPSF:

```julia
# Create Zernike coefficients object with up to 21 terms
zc = ZernikeCoefficients(21)

# Add common aberrations
add_defocus!(zc, 1.0)           # 1 wave of defocus
add_astigmatism!(zc, 0.5, π/4)  # 0.5 waves of astigmatism at 45°
add_coma!(zc, 0.3)              # 0.3 waves of coma
add_spherical!(zc, 0.2)         # 0.2 waves of spherical aberration

# Create PSF with these aberrations
psf = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33, base_zernike=zc)

# Update aberrations for an existing PSF
add_spherical!(psf.zernike_coeffs, 0.3)  # Add more spherical aberration
update_pupils!(psf)  # Update the pupil function to reflect the changes
```

## Performance Considerations

The Vector3DPSF is the most computationally intensive model in MicroscopePSFs.jl:

- Significantly slower than all other PSF models
- Performance scales with the resolution of the pupil function (grid_size parameter)
- Calculating the vector pupils is computationally intensive but only needs to be done once
- For repeated evaluations at many positions, the pupil computation is amortized
- For performance-critical applications, consider creating a `SplinePSF` from a Vector3DPSF:

```julia
# Create a SplinePSF for faster repeated evaluations
x_range = range(-2, 2, length=81)  # μm
y_range = range(-2, 2, length=81)  # μm
z_range = range(-2, 2, length=41)  # μm

# This step is slow but only needs to be done once
spline_psf = SplinePSF(psf, x_range, y_range, z_range)

# Now evaluations are much faster
intensity = spline_psf(0.1, 0.2, 0.3)
```

## Dipole Orientations and Polarization

The Vector3DPSF model allows for simulating different dipole orientations, which is crucial for applications like single-molecule orientation studies:

```julia
# Create PSFs for different dipole orientations
orientations = [
    ("X-oriented", DipoleVector(1.0, 0.0, 0.0)),
    ("Y-oriented", DipoleVector(0.0, 1.0, 0.0)),
    ("Z-oriented", DipoleVector(0.0, 0.0, 1.0)),
    ("XY-plane 45°", DipoleVector(0.707, 0.707, 0.0)),
    ("3D 45°", DipoleVector(0.577, 0.577, 0.577))
]

# Compare how different orientations appear at different z positions
z_positions = [-1.0, 0.0, 1.0]  # μm
x = y = range(-1, 1, length=51)  # μm

fig = Figure(size=(1200, 800))

for (i, (name, dipole)) in enumerate(orientations)
    psf = Vector3DPSF(1.4, 0.532, dipole, n_medium=1.33)
    
    for (j, z) in enumerate(z_positions)
        ax = Axis(fig[i, j], aspect=DataAspect(),
                 title=j==1 ? name : "z = $(z)μm",
                 xlabel=i==length(orientations) ? "x (μm)" : "",
                 ylabel=j==1 ? "y (μm)" : "")
                 
        img = [psf(xi, yi, z) for yi in y, xi in x]
        heatmap!(ax, x, y, img, colormap=:viridis)
    end
end

fig
```

## Refractive Index Interfaces

The Vector3DPSF model accounts for refractive index mismatches at the interfaces:

```julia
# Compare different refractive index configurations
configs = [
    ("Matched (water)", 1.33, 1.33, 1.33),     # All water
    ("Oil immersion", 1.33, 1.518, 1.518),    # Water sample, oil immersion
    ("Glycerol immersion", 1.33, 1.47, 1.47)  # Water sample, glycerol immersion
]

fig = Figure(size=(900, 300))

for (i, (name, n_medium, n_coverslip, n_immersion)) in enumerate(configs)
    psf = Vector3DPSF(
        1.4, 0.532, dipole_z,
        n_medium=n_medium,
        n_coverslip=n_coverslip,
        n_immersion=n_immersion
    )
    
    ax = Axis(fig[1, i], aspect=DataAspect(),
              title=name,
              xlabel="x (μm)", ylabel=i==1 ? "y (μm)" : "")
              
    img = [psf(xi, yi, 0.0) for yi in y, xi in x]
    heatmap!(ax, x, y, img, colormap=:viridis)
end

fig
```

## Limitations

The Vector3DPSF model, while most physically accurate, has several limitations:

1. **Computational Intensity**: Significantly more computationally demanding than other PSF models
2. **Memory Requirements**: Requires more memory for storing vector pupil functions
3. **Simplified Interfaces**: Uses a simplified model of the objective/coverslip/sample interfaces
4. **Limited Defocus Range**: Accuracy decreases at large defocus values
5. **Fixed Dipole Orientation**: Each PSF instance has a fixed dipole orientation

## Relationship to Other PSFs

The Vector3DPSF is related to other PSF models in the following ways:

1. **Scalar3DPSF**: Vector3DPSF reduces to Scalar3DPSF when polarization effects are negligible (low NA)
2. **Airy2D**: For z=0 and low NA, the Vector3DPSF approaches the Airy pattern
3. **SplinePSF**: For performance-critical applications, create a SplinePSF from a Vector3DPSF

## When to Use Vector3DPSF

The Vector3DPSF model is recommended for:

1. High-NA objectives (NA > 1.2)
2. Polarization-sensitive applications
3. Applications requiring dipole emission pattern modeling
4. Systems with significant refractive index mismatches
5. Research requiring the highest level of physical accuracy