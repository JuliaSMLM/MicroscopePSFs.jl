# Vector3DPSF

The `Vector3DPSF` model implements a comprehensive three-dimensional point spread function based on vectorial diffraction theory. It accounts for polarization effects, high-NA phenomena, dipole emission patterns, and arbitrary aberrations, providing the highest level of physical accuracy among the PSF models. This model is particularly important for high-NA objectives and applications where polarization effects cannot be ignored.

## Mathematical Model

The Vector3DPSF implements a full vectorial diffraction model based on the Richards-Wolf vector diffraction theory, accounting for dipole emission patterns, polarization effects, refractive index interfaces, and optical aberrations.

### Complete PSF Formula

$$PSF(x - x_i, y - y_i, z_i, z_s) = \sum_{m=x,y} \sum_{n=p_x, p_y, p_z} |F[h(k_x, k_y)w_{m,n}A(k_x, k_y)e^{\iota 2\pi(k_x x_i+k_y y_i+k_{z\mathrm{med}} z_i-k_{z\mathrm{imm}}z_s)}]|^2$$

Where:

-  $(x_i, y_i, z_i)$ is the emitter position relative to the focal plane 
-  $z_s$ is the focal plane relative to the coverslip (nominal focal plane) 
-  $F$ denotes the Fourier transform operation 
-  $h(k_x, k_y)$ is the complex pupil function incorporating aberrations 
-  $w_{m,n}$ represents the electric field components at the pupil plane from a fixed dipole and incorporates Fresnel coefficients 
-  $A(k_x, k_y)$ is the apodization factor for energy conservation 
-  $k_{z\mathrm{med}}$ and $k_{z\mathrm{imm}}$ are the z-components of the wave vectors in the sample and immersion media respectively 

### Electric Field Components

The vectorial field components $w_{m,n}$ at the pupil plane are calculated as:

$$w_{x,n} = P_n \cos(\phi) - S_n \sin(\phi)$$
$$w_{y,n} = P_n \sin(\phi) + S_n \cos(\phi)$$

Where $\phi$ is the angular component in the polar coordinate of the frequency space, and $P_n$ and $S_n$ are electric field components in $p$ and $s$ polarizations relative to the incident plane at the sample space:

$$\begin{align}
P_{px} &= T_p \cos \theta_1 \cos \phi \\
P_{py} &= T_p \cos \theta_1 \sin \phi \\
P_{pz} &= -T_p \sin \theta_1 \\
S_{px} &= -T_s \sin \phi \\
S_{py} &= T_s \cos \phi \\
S_{pz} &= 0
\end{align}$$

Where $P_x$, $P_y$, $P_z$ are the cartesian components of the dipole moments.

### Fresnel Coefficients

The total transmission coefficients $T_p$ and $T_s$ for p- and s-polarized light account for the interfaces between different media (sample, coverslip, and immersion medium):

$$\begin{align}
T_P &= \tau_{P12} \tau_{P23} \\
T_S &= \tau_{S12} \tau_{S23}
\end{align}$$

With interface-specific Fresnel transmission coefficients:

$$\begin{align}
\tau_{Pij} &= \frac{2n_i\cos\theta_i}{n_i\cos\theta_j + n_j\cos\theta_i} \\
\tau_{Sij} &= \frac{2n_i\cos\theta_i}{n_i\cos\theta_i + n_j\cos\theta_j}
\end{align}$$

Where $n_i$ and $\theta_i$ are the refractive index and the light propagation angle in medium $i$, and the subscriptions 1, 2, 3 denote the sample medium, the coverslip, and the immersion medium.

### Apodization Factor

The apodization factor $A(k_x, k_y)$ maintains energy conservation across media interfaces:

$$A(k_x, k_y) = \frac{\sqrt{\cos\theta_{\mathrm{imm}}}}{\cos\theta_{\mathrm{med}}}$$

Where:

$$\cos\theta_{\mathrm{med}} = \sqrt{1 - \frac{(k_x^2 + k_y^2)\lambda^2}{4\pi^2 n_{\mathrm{med}}^2}}$$

$$\cos\theta_{\mathrm{imm}} = \sqrt{1 - \frac{(k_x^2 + k_y^2)\lambda^2}{4\pi^2 n_{\mathrm{imm}}^2}}$$

This factor accounts for the change in solid angle that occurs during refraction, critical for accurately modeling high-NA systems.

### Wave Vectors

The z-components of the wave vectors in the different media are:

$$k_{z\mathrm{med}} = \sqrt{\frac{n^2_{\mathrm{med}}}{\lambda^2} - \frac{n^2_{\mathrm{imm}}}{\lambda^2}\sin^2 \theta_{\mathrm{imm}}}$$

$$k_{z\mathrm{imm}} = \frac{n_{\mathrm{imm}}}{\lambda}\cos \theta_{\mathrm{imm}}$$

When $\theta_{\mathrm{imm}}$ exceeds the critical angle, $k_{z\mathrm{med}}$ becomes imaginary, resulting in evanescent waves and giving rise to supercritical angle fluorescence (SAF) effects.

## Constructor Options

```julia
Vector3DPSF(nₐ::Real, λ::Real, dipole::DipoleVector;
            base_pupil::Union{Nothing, PupilFunction}=nothing,
            base_zernike::Union{Nothing, ZernikeCoefficients}=nothing,
            n_medium::Real=1.33,
            n_coverslip::Real=1.52,
            n_immersion::Real=1.52,
            focal_z::Real=0.0, 
            grid_size::Integer=128)
```

### Required Parameters

- `nₐ`: Numerical aperture of the objective
- `wavelength`: Wavelength of light in microns
- `dipole`: Dipole orientation vector specified using `DipoleVector`

### Optional Parameters

- `base_pupil`: Optional base aberration pupil function (default: none)
- `base_zernike`: Zernike coefficients for aberrations (default: none)
- `n_medium`: Refractive index of the sample medium (default: 1.33, water)
- `n_immersion`: Refractive index of the immersion medium (default: 1.52, oil)
- `n_coverslip`: Refractive index of the coverslip (default: 1.52, glass)
- `focal_z`: Focal plane position in microns (default: 0.0)
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
    n_immersion=1.52   # Immersion oil
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
    ("Oil immersion", 1.33, 1.52, 1.52),    # Water sample, oil immersion
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

## Supercritical Angle Fluorescence

A notable feature of the Vector3DPSF model is its ability to accurately capture Supercritical Angle Fluorescence (SAF), which occurs when light emitted by fluorophores near the coverslip is collected at angles exceeding the critical angle:

```julia
# Compare SAF effect at different distances from the coverslip
z_positions = [0.0, 0.1, 0.3, 0.5, 1.0]  # μm from coverslip
x = y = range(-1, 1, length=51)  # μm

fig = Figure(size=(1200, 300))

# Create PSF with oil immersion (favorable for SAF)
psf = Vector3DPSF(
    1.4, 0.532, dipole_z,
    n_medium=1.33,   # Water
    n_immersion=1.52  # Oil
)

for (i, z) in enumerate(z_positions)
    ax = Axis(fig[1, i], aspect=DataAspect(),
             title="z = $(z)μm from coverslip",
             xlabel="x (μm)", 
             ylabel=i==1 ? "y (μm)" : "")
             
    # Emitter at z from coverslip, focus at z
    img = [psf(xi, yi, 0.0) for yi in y, xi in x]
    hm = heatmap!(ax, x, y, img, colormap=:viridis)
    
    # Add radial profile
    r = range(0, 1, length=50)
    intensity = [psf(ri, 0.0, 0.0) for ri in r]
    lines!(ax, r, intensity .* 1.5, color=:red, linewidth=2)
end

fig
```

As the emitter moves away from the coverslip, you'll observe the SAF contribution decreasing, which appears as a change in both the total intensity and the shape of the PSF.

## Performance Considerations

The Vector3DPSF is the most computationally intensive model in MicroscopePSFs.jl:

- Significantly slower than all other PSF models
- Performance scales with the resolution of the pupil function (grid_size parameter)
- Calculating the vector pupils is computationally intensive but only needs to be done once
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
6. Studies involving near-interface fluorescence effects like SAF