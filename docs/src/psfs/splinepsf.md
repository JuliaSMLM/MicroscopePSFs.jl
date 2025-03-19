# SplinePSF

The `SplinePSF` model provides an efficient representation of point spread functions using B-spline interpolation. Unlike other PSF models that are based on physical principles, SplinePSF is a computational acceleration technique that can significantly speed up PSF evaluations for complex models, or represent experimentally measured PSFs.

## Mathematical Model

The SplinePSF uses B-spline interpolation to rapidly evaluate a pre-computed PSF grid:

```math
I(x, y, z) = \sum_{i,j,k} c_{i,j,k} \beta^n(x - x_i) \beta^n(y - y_j) \beta^n(z - z_k)
```

where:
- ``c_{i,j,k}`` are the B-spline coefficients
- ``\beta^n`` is the B-spline basis function of order n (typically cubic, n=3)
- ``x_i, y_j, z_k`` are the knot points of the spline

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

## Key Features

- **Performance Acceleration**: Significantly faster evaluation of complex PSFs
- **Experimental PSF Support**: Can represent measured PSF data from calibration beads
- **Smooth Interpolation**: Cubic splines provide continuous second derivatives
- **Flexible Precision**: Trade-off between accuracy and speed via grid density

## B-Spline Interpolation Orders

The SplinePSF uses different interpolation orders with specific trade-offs:

| Order | Name | Properties | Use Case |
|-------|------|------------|----------|
| 0 | Constant | Simple nearest neighbor, no continuity | Very fast, low accuracy |
| 1 | Linear | Continuous function, discontinuous derivatives | Good compromise |
| 3 | Cubic | Continuous second derivatives | Best accuracy (default) |

## Examples

Creating a SplinePSF:

```julia
# Create a source PSF to accelerate
scalar_psf = ScalarPSF(1.4, 0.532, 1.518)

# Define coordinate sampling grid
x_range = y_range = range(-2.0, 2.0, step=0.1)  # μm, 0.1μm step size
z_range = range(-2.0, 2.0, step=0.2)  # μm, 0.2μm step size

# Create a spline representation (this step may take time)
spline_psf = SplinePSF(scalar_psf, x_range, y_range, z_range)

# Alternative: use the convenience constructor with auto-ranges
spline_psf_auto = SplinePSF(scalar_psf, 
                           lateral_range=2.0, 
                           axial_range=1.0,
                           lateral_step=0.1, 
                           axial_step=0.2)
```

## Performance Considerations

- **Initial Construction**: Creating the SplinePSF can be slow for complex source PSFs
- **Evaluation Speed**: Once created, evaluation is very fast (typically 10-100× faster than source PSF)
- **Memory Usage**: Memory scales with grid size (O(nx×ny×nz))
- **Interpolation Order**: Higher order provides more accuracy at slight computational cost
- **Boundary Handling**: Returns 0 for positions outside the grid boundaries

## Working with Experimental PSFs

SplinePSF can be used with experimentally measured PSF data:

```julia
# Define physical coordinates for the experimental PSF
pixel_size = 0.1  # μm in the image
z_step = 0.2      # μm between z-planes
measured_psf = load_experimental_psf("path/to/measured_psf.h5")

# Create coordinate ranges
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

## Limitations

1. **Sampling Range**: Limited to positions within the sampled grid range
2. **Memory Requirements**: Can use significant memory for large, dense grids
3. **Initial Computation**: Requires upfront computation time to create
4. **No Physics Propagation**: Cannot directly model physical effects not in the source PSF
5. **Phase Information**: Limited phase accuracy in complex field representation

For standard usage patterns, camera integration, and comparison with other PSF types, see the [PSF Overview](overview.md).