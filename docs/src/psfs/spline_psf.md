# SplinePSF

The `SplinePSF` model provides an efficient representation of point spread functions using B-spline interpolation. Unlike other PSF models that are based on physical principles, SplinePSF is a computational acceleration technique that can significantly speed up PSF evaluations for complex models like Vector3DPSF or Scalar3DPSF. It can also be used to represent experimentally measured PSFs from calibration data.

## Mathematical Model

The SplinePSF uses B-spline interpolation to rapidly evaluate a pre-computed PSF grid:

```math
I(x, y, z) = \sum_{i,j,k} c_{i,j,k} \beta^n(x - x_i) \beta^n(y - y_j) \beta^n(z - z_k)
```

where:
- ``c_{i,j,k}`` are the B-spline coefficients
- ``\beta^n`` is the B-spline basis function of order n (typically cubic, n=3)
- ``x_i, y_j, z_k`` are the knot points of the spline

This representation allows for:
1. Continuous representation with smooth derivatives
2. Fast evaluation through efficient lookup and computation
3. Accurate interpolation between sampled grid points

## Constructor Options

SplinePSF offers several constructor methods for different use cases:

### From Another PSF Model (Most Common)

```julia
# 3D PSF
SplinePSF(psf::AbstractPSF, 
          x_range::AbstractRange,
          y_range::AbstractRange,
          z_range::AbstractRange;
          order::Integer=3)

# 2D PSF
SplinePSF(psf::AbstractPSF, 
          x_range::AbstractRange,
          y_range::AbstractRange;
          order::Integer=3)
```

### From Pre-computed PSF Data

```julia
# From 3D PSF stack
SplinePSF(psf_stack::AbstractArray{<:Real,3}, 
          x_range::AbstractRange,
          y_range::AbstractRange,
          z_range::AbstractRange;
          order::Integer=3)

# From 2D PSF image
SplinePSF(psf_image::AbstractArray{<:Real,2}, 
          x_range::AbstractRange,
          y_range::AbstractRange;
          order::Integer=3)
```

### Convenience Constructor with Auto-ranges

```julia
SplinePSF(psf::AbstractPSF; 
          lateral_range::Float64=2.0,
          axial_range::Float64=1.0,
          lateral_step::Float64=0.05,
          axial_step::Float64=0.1,
          order::Integer=3)
```

### Parameters

- `psf`: Source PSF to sample or pre-computed PSF data
- `x_range`, `y_range`, `z_range`: Coordinate ranges for the grid points
- `order`: Interpolation order (default: 3 for cubic B-splines)
- `lateral_range`: Half-width of lateral (xy) sampling range in microns
- `axial_range`: Half-width of axial (z) sampling range in microns
- `lateral_step`: Step size in microns for lateral sampling
- `axial_step`: Step size in microns for axial sampling

## Basic Usage

### Creating a SplinePSF

```julia
using MicroscopePSFs

# Create a source PSF to accelerate
dipole_z = DipoleVector(0.0, 0.0, 1.0)
vector_psf = Vector3DPSF(
    1.4,                # Numerical aperture
    0.532,              # Wavelength in microns
    dipole_z,           # Dipole orientation
    n_medium=1.33       # Sample medium refractive index
)

# Define coordinate sampling grid with explicit step sizes
x_range = range(-2.0, 2.0, step=0.1)  # μm, 0.1μm step size
y_range = range(-2.0, 2.0, step=0.1)  # μm, 0.1μm step size
z_range = range(-2.0, 2.0, step=0.2)  # μm, 0.2μm step size

# Create a spline representation (this step may take time)
spline_psf = SplinePSF(vector_psf, x_range, y_range, z_range)

# Alternatively, use the convenience constructor with auto-ranges
spline_psf_auto = SplinePSF(vector_psf, 
                           lateral_range=2.0, 
                           axial_range=1.0,
                           lateral_step=0.1, 
                           axial_step=0.2)
```

### Direct Evaluation

```julia
# Evaluate at specific 3D position
x = 0.1  # μm
y = 0.2  # μm
z = 0.5  # μm
intensity = spline_psf(x, y, z)

# For 2D PSFs or z=0 evaluation
intensity_2d = spline_psf(x, y)

# Get complex amplitude (returns sqrt of intensity as complex value)
amp = amplitude(spline_psf, x, y, z)
```

### Creating PSF Images

```julia
# Create a grid of positions with explicit step size
x = range(-1, 1, step=0.02)  # μm, 0.02μm step size
y = range(-1, 1, step=0.02)  # μm, 0.02μm step size
z = 0.0  # focal plane

# Compute PSF values at each position
intensity_2d = [spline_psf(xi, yi, z) for yi in y, xi in x]

# Visualize with CairoMakie
using CairoMakie
fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1], aspect=DataAspect(),
          title="SplinePSF Interpolation (z=0μm)",
          xlabel="x (μm)", ylabel="y (μm)")
hm = heatmap!(ax, x, y, intensity_2d, colormap=:viridis)
Colorbar(fig[1, 2], hm)
fig
```

## Integration with Camera

```julia
# Create camera with 100nm pixels (20×20 pixel grid)
pixel_size = 0.1  # μm
camera = IdealCamera(1:20, 1:20, pixel_size)

# Create 3D emitter at position (1μm, 1μm, 0.5μm) with 1000 photons
emitter = Emitter3D(1.0, 1.0, 0.5, 1000.0)  # x, y, z, photons

# Integrate PSF over pixels with 2×2 subsampling
pixels = integrate_pixels(spline_psf, camera, emitter, sampling=2)

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

## Performance Comparison

One of the main purposes of SplinePSF is to accelerate computation. Here's how to compare performance:

```julia
using MicroscopePSFs
using BenchmarkTools

# Create computationally expensive PSF
dipole_z = DipoleVector(0.0, 0.0, 1.0)
vector_psf = Vector3DPSF(1.4, 0.532, dipole_z, n_medium=1.33)

# Create a spline version with same parameters
x_range = range(-2.0, 2.0, length=41)
y_range = range(-2.0, 2.0, length=41)
z_range = range(-2.0, 2.0, length=21)
spline_psf = SplinePSF(vector_psf, x_range, y_range, z_range)

# Compare evaluation speed for single point
@btime $vector_psf(0.5, 0.5, 0.1)  # Slow
@btime $spline_psf(0.5, 0.5, 0.1)  # Fast

# Compare speed for camera integration
pixel_size = 0.1  # μm
camera = IdealCamera(1:16, 1:16, pixel_size)
emitter = Emitter3D(1.0, 1.0, 0.0, 1000.0)

@btime integrate_pixels($vector_psf, $camera, $emitter)  # Slow
@btime integrate_pixels($spline_psf, $camera, $emitter)  # Fast
```

Typical speedups range from 10-100× depending on the complexity of the original PSF model and grid size.

## Working with Experimental PSFs

SplinePSF can also be used with experimentally measured PSF data:

```julia
# Load experimental PSF data
# (In practice, this would be loaded from a file)
function load_experimental_psf(filename)
    # Example function - actual implementation would depend on file format
    # Here we assume you have a 3D stack in HDF5 format
    using HDF5
    h5open(filename, "r") do file
        return read(file, "psf_stack")
    end
end

# Define physical coordinates for the experimental PSF
pixel_size = 0.1  # μm in the image
z_step = 0.2      # μm between z-planes
measured_psf = load_experimental_psf("path/to/measured_psf.h5")

# Create coordinate ranges based on PSF stack dimensions
nx, ny, nz = size(measured_psf)
x_center = nx ÷ 2
y_center = ny ÷ 2
z_center = nz ÷ 2

x_range = range(-(x_center) * pixel_size, (nx-x_center-1) * pixel_size, length=nx)
y_range = range(-(y_center) * pixel_size, (ny-y_center-1) * pixel_size, length=ny)
z_range = range(-(z_center) * z_step, (nz-z_center-1) * z_step, length=nz)

# Create the SplinePSF
experimental_psf = SplinePSF(measured_psf, x_range, y_range, z_range)
```

## B-Spline Interpolation Details

The SplinePSF uses different interpolation orders with specific trade-offs:

| Order | Name | Properties | Use Case |
|-------|------|------------|----------|
| 0 | Constant | Simple nearest neighbor, no continuity | Very fast, low accuracy |
| 1 | Linear | Continuous function, discontinuous derivatives | Good compromise |
| 3 | Cubic | Continuous second derivatives | Best accuracy (default) |

Cubic splines (order=3) provide:
1. Continuous second derivatives for smooth interpolation
2. Local support for efficient computation
3. Good approximation of smooth functions
4. Analytic derivatives for optimization
5. Minimal memory requirements compared to accuracy

## Accuracy Considerations

The accuracy of SplinePSF depends on the density of the sampling grid:

```julia
# Compare SplinePSF with original PSF
scalar_psf = Scalar3DPSF(1.4, 0.532, 1.518)

# Create a SplinePSF with recommended sampling resolution
spline_psf = SplinePSF(scalar_psf, 
                      range(-2, 2, step=0.05),  # 0.05μm steps in x
                      range(-2, 2, step=0.05),  # 0.05μm steps in y
                      range(-1, 1, step=0.1))   # 0.1μm steps in z

# Create test positions
x_eval = range(-1, 1, step=0.02)  # Fine evaluation grid for comparison
y = 0.0
z = 0.0

# Calculate profiles
original = [scalar_psf(xi, y, z) for xi in x_eval]
spline = [spline_psf(xi, y, z) for xi in x_eval]

# Visualize comparison
using CairoMakie
fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], xlabel="Position (μm)", ylabel="Intensity",
          title="SplinePSF Accuracy Comparison")

lines!(ax, x_eval, original, label="Original PSF", linewidth=2)
lines!(ax, x_eval, spline, label="SplinePSF (0.05μm xy, 0.1μm z steps)", 
       linestyle=:dash, linewidth=2)

axislegend(ax)
fig
```

## Performance Considerations

SplinePSF is designed for computational efficiency:

- **Initial Construction**: Creating the SplinePSF can be slow for complex source PSFs
- **Evaluation Speed**: Once created, evaluation is very fast (typically 10-100× faster than source PSF)
- **Memory Usage**: Memory scales with grid size (O(nx×ny×nz))
- **Interpolation Order**: Higher order provides more accuracy at slight computational cost
- **Boundary Handling**: Returns 0 for positions outside the grid boundaries

## Limitations

The SplinePSF model has several limitations:

1. **Sampling Range**: Limited to positions within the sampled grid range
2. **Memory Requirements**: Can use significant memory for large, dense grids
3. **Initial Computation**: Requires upfront computation time to create
4. **No Physics Propagation**: Cannot directly model physical effects not in the source PSF
5. **Phase Information**: Limited phase accuracy in complex field representation

## When to Use SplinePSF

Consider using SplinePSF when:

1. **Performance is Critical**: For fitting or simulation tasks requiring many PSF evaluations
2. **Working with Complex Models**: To accelerate Vector3DPSF or aberrated Scalar3DPSF evaluations
3. **Using Experimental PSFs**: When working with measured PSF data
4. **Running on Limited Resources**: To make complex models feasible on modest hardware
5. **Building Interactive Applications**: Where real-time performance is necessary