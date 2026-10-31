# Multi-channel PSFs

`PlanePSF` shifts an emitter's axial coordinate, `AnalyzerPSF` selects a vector
field component, and `FourPiPSF` interferes two opposed objectives. Each is a
complete 3D `AbstractPSF` that can be evaluated, integrated, tabulated as a
`SplinePSF`, and saved. `MultiChannelPSF` holds complete channel PSFs and their
registration and photon fractions. It deliberately does not subtype
`AbstractPSF`, so a one-channel fitter must handle channels explicitly.

## Coordinates and the forward law

The one definition of the container is

```text
expected image i = N * mc.transmission[i] *
    integrate_pixels(mc.channels[i], camera_i, unit-photon emitter at channel_position(mc, i, x, y, z))
```

`N` is the emitter's total photons. Registration matrices have size 2×3 and map
**absolute camera coordinates in microns**, as used by `SMLMData` pixel edges:
`[xᵢ, yᵢ] = A * [x, y] + t`. They do not map pixels or crop-local coordinates;
keeping the original pixel edges makes the same registration valid after crops.
The axial law is `zᵢ = z + mc.z_offset[i]`. The reference channel must have
identity registration and zero axial offset. Missing transforms, fractions,
and axial offsets default to identity, equal splitting, and zeros.

Each channel must integrate to one alone. For 4Pi phase channels this means
unit energy **on average over the fringe**, rather than at every emitter z.
Transmission is a nonnegative fraction of one total photon budget, with sum
at most one. The constructor cannot check the channel normalization. In
particular, raw H/V `AnalyzerPSF`s carry the polarization split themselves and
are refused, including through `PlanePSF` and either arm of `FourPiPSF`.
Tabulate these channels into a `SplinePSF` scaled to unit integral and put
the split in `transmission`.

## Dual focal plane

```julia
using MicroscopePSFs

p = ScalarPSF(1.2, 0.68, 1.33)
mc = MultiChannelPSF(; channels=[p, PlanePSF(p; dz=0.2)],
                     transmission=[0.5, 0.5])
camera = IdealCamera(collect(0.0:0.1:2.0), collect(0.0:0.1:2.0))
N = 1000.0
images = [mc.transmission[i] * integrate_pixels(mc.channels[i], camera,
    Emitter3D(channel_position(mc, i, 1.0, 1.0, 0.1)..., N)) for i in 1:2]
```

`PlanePSF(p; dz)` evaluates every method at `z + dz`. This is emitter-side
high-NA defocus. It is not the same phase as translating a camera in low-NA
image space; a learned pupil absorbs that difference. Use either a plane
wrapper or a container `z_offset` for a given shift to avoid adding it twice.
A plane can wrap a 3D intensity spline, but splines cannot supply coherent
fields to an analyzer or an interferometer.

## H/V polarization

An analyzer passes `Ex` for `component=1` (H) or `Ey` for `component=2` (V),
zeroing the rejected component for every dipole. It does not renormalize, so
`h(x,y,z) + v(x,y,z) ≈ p(x,y,z)`. Rotating dipoles remain mutually incoherent.
`dipole_fields` preserves that distinction, returning one component vector
per dipole. `amplitude` is meaningful only for a single dipole.

This example tabulates a focal slice and normalizes each channel alone.
Increase the field of view for more accurate collected photon fractions.
For a 3D calibration, normalize each channel using its calibrated energy
convention, preserving the relative intensities of its axial slices.

```julia
using MicroscopePSFs

p = VectorPSF(1.2, 0.68; n_medium=1.33, n_coverslip=1.33,
              n_immersion=1.33, grid_size=64)
h, v = AnalyzerPSF(p; component=1), AnalyzerPSF(p; component=2)
r = -1.0:0.05:1.0
area = step(r)^2
grids = [[a(x, y, 0.0) for y in r, x in r] for a in (h, v)]
energies = sum.(grids) .* area
channels = [SplinePSF(g ./ e, r, r) for (g, e) in zip(grids, energies)]
mc = MultiChannelPSF(; channels, transmission=energies ./ sum(energies))
```

## 4Pi interference

```julia
using MicroscopePSFs

arm = ScalarPSF(1.2, 0.68, 1.33)
p₀ = FourPiPSF(arm, arm; phase=0.0)
pπ = FourPiPSF(arm, arm; phase=π)
mc = MultiChannelPSF(; channels=[p₀, pπ], transmission=[0.5, 0.5])
intensity = mc.channels[1](0.1, 0.0, 0.2)
```

The objectives share a focal plane; the top sees `z` and the bottom sees `-z`.
The two images overlap with the same orientation. One arm therefore has one
more mirror than the other, and both arms' camera x and y components lie in a
common frame. For dipole `d`,

```math
E_d = \frac{E_{t,d}(x,y,z) + e^{i\phi}s_d E_{b,d}(x,y,-z)}{\sqrt{2}},
\qquad p = \sum_d \lVert E_d\rVert^2.
```

`phase` is the cavity phase. Polarization-dependent cavity phases use one
`AnalyzerPSF(FourPiPSF(top, bottom; phase=φ); component=c)` per channel,
followed by the same channel normalization convention used for H/V.
Vector arms require `n_medium ≥ NA`: otherwise evanescent fields grow
without bound at negative emitter z, which one of the arms always sees.

The factor `1/√2` makes unit-energy arms give a unit-energy channel averaged
over the fringe. Complementary phase channels obey
`p(φ) + p(φ + π) ≈ top(x,y,z) + bottom(x,y,-z)` when the mirror signs are all
ones or for a rotating dipole.

### Why the axial dipole changes sign

A dipole radiates the transverse far field `p - dot(p, k̂) * k̂`. A lens rotates
this field onto its optical axis. Reflection of the opposed geometry changes
`pz` to `-pz` while leaving `px` and `py` unchanged. After the arm mirrors align
the camera components, the x and y dipoles are even and the axial dipole is
odd. A rotating dipole thus needs signs `(1, 1, -1)` in x, y, z dipole order.
Its three intensities add incoherently, but each dipole's two arms interfere.

For a fixed dipole, build the bottom arm with `DipoleVector(px, py, -pz)`;
the mirror is already in its orientation, so its remaining field sign is +1.
The constructor compares normalized dipoles, including dipoles constructed
with same-typed arguments that did not normalize automatically.

Emitter z and dipole z share one frame, with +z pointing into the sample,
away from the objective. A code dipole `DipoleVector(1.0, 0, 1.0)` matches
lab orientation `[1, 0, 1]/√2`. The opposed-objective mirror rule still applies:

```julia
using MicroscopePSFs

optics = (; n_medium=1.33, n_coverslip=1.33, n_immersion=1.33, grid_size=64)
top = VectorPSF(1.2, 0.68, DipoleVector(1.0, 0, 1.0); optics...)
bottom = VectorPSF(1.2, 0.68, DipoleVector(1.0, 0, -1.0); optics...)
p = FourPiPSF(top, bottom; phase=0.3)
fields = dipole_fields(p, 0.1, 0.0, 0.2)
```
