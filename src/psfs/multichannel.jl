# Complex-field wrappers and the multi-channel forward model.

_psf_float(::Abstract3DPSF{T}) where {T} = T
function _psf_float(psf::SplinePSF{T}) where {T}
    isnothing(psf.z_range) && throw(ArgumentError(
        "$(typeof(psf)) has no z axis; use a 3D PSF"))
    return T
end
_psf_float(psf::AbstractPSF) = throw(ArgumentError(
    "$(typeof(psf)) has no z axis; use a 3D PSF"))

"""
    PlanePSF(psf::AbstractPSF; dz::Real)

Create a detection plane focused `dz` microns away.

# Arguments
- `psf`: A 3D PSF, including a 3D `SplinePSF`.
- `dz`: Required emitter-side focus shift in microns.

# Returns
A 3D PSF whose intensity, amplitude, and dipole fields evaluate `psf` at
`z + dz`. This is emitter-side high-NA defocus, which is not the same phase
as moving a camera in low-NA image space. A learned pupil absorbs the difference.

# Examples
```julia
plane = PlanePSF(ScalarPSF(1.2, 0.68, 1.33); dz=0.2)
intensity = plane(0.0, 0.0, 0.1)
```
"""
struct PlanePSF{T, P <: AbstractPSF} <: Abstract3DPSF{T}
    psf::P
    dz::T
    function PlanePSF(psf::P; dz::Real) where {P <: AbstractPSF}
        T = _psf_float(psf)
        isfinite(dz) && isfinite(T(dz)) || throw(ArgumentError(
            "dz must be finite; supply a finite focus shift"))
        return new{T,P}(psf, T(dz))
    end
end
(p::PlanePSF)(x::Real, y::Real, z::Real) = p.psf(x, y, z + p.dz)
amplitude(p::PlanePSF, x::Real, y::Real, z::Real) = amplitude(p.psf, x, y, z + p.dz)
dipole_fields(p::PlanePSF, x::Real, y::Real, z::Real) = dipole_fields(p.psf, x, y, z + p.dz)

_n_components(::AbstractPSF) = 0
_n_components(::ScalarPSF) = 1
_n_components(::VectorPSF) = 2
_n_components(p::PlanePSF) = _n_components(p.psf)

"""
    AnalyzerPSF(psf::AbstractPSF; component::Integer)

Pass one camera-frame polarization component of a vector field.

# Arguments
- `psf`: `VectorPSF`, or a `PlanePSF`, `FourPiPSF`, or `AnalyzerPSF` over one.
- `component`: Required component: 1 passes `Ex` (H), 2 passes `Ey` (V).

# Returns
A 3D PSF with the rejected component set to zero in every dipole field.
It is not renormalized: H and V intensities sum to the input intensity.
`amplitude` returns `[Ex, Ey]` for one dipole and errors for several dipoles.

# Examples
```julia
psf = VectorPSF(1.2, 0.68)
h = AnalyzerPSF(psf; component=1)
v = AnalyzerPSF(psf; component=2)
intensity = h(0.1, 0.0, 0.0) + v(0.1, 0.0, 0.0)
```
"""
struct AnalyzerPSF{T, P <: AbstractPSF} <: Abstract3DPSF{T}
    psf::P
    component::Int
    function AnalyzerPSF(psf::P; component::Integer) where {P <: AbstractPSF}
        _n_components(psf) == 2 || throw(ArgumentError(
            "An analyzer needs a vector field; use a VectorPSF or a wrapper over one"))
        component in (1, 2) || throw(ArgumentError(
            "component must be 1 (H) or 2 (V); select one of these components"))
        return new{_psf_float(psf),P}(psf, Int(component))
    end
end
_n_components(p::AnalyzerPSF) = _n_components(p.psf)
(p::AnalyzerPSF)(x::Real, y::Real, z::Real) =
    sum(abs2(E[p.component]) for E in dipole_fields(p.psf, x, y, z))
function dipole_fields(p::AnalyzerPSF, x::Real, y::Real, z::Real)
    fields = dipole_fields(p.psf, x, y, z)
    for E in fields
        E[3 - p.component] = zero(eltype(E))
    end
    return fields
end

function _single_dipole_field(fields)
    length(fields) == 1 || error(
        "amplitude() not meaningful for multiple dipole orientations; " *
        "use intensity evaluation directly")
    return only(fields)
end
function _single_dipole_amplitude(fields, ::Val{1})
    return only(_single_dipole_field(fields))
end
function _single_dipole_amplitude(fields, ::Val{2})
    return _single_dipole_field(fields)
end
amplitude(p::AnalyzerPSF, x::Real, y::Real, z::Real) =
    _single_dipole_amplitude(dipole_fields(p, x, y, z), Val(_n_components(p)))

_mirror_signs(::ScalarPSF) = (1,)
function _mirror_signs(p::VectorPSF)
    count = length(p.vector_pupils)
    count == 1 && return (1,)
    count == 3 && return (1, 1, -1)
    throw(ArgumentError(
        "A 4Pi VectorPSF arm has $count pupils; " *
        "use 1 fixed-dipole pupil or 3 rotating-dipole pupils"))
end
_mirror_signs(p::Union{PlanePSF,AnalyzerPSF}) = _mirror_signs(p.psf)
_mirror_signs(p::AbstractPSF) = throw(ArgumentError(
    "A 4Pi arm needs a complex field; " *
    "replace $(typeof(p)) with a ScalarPSF or VectorPSF arm"))
_arm_wavelength(p::Union{ScalarPSF,VectorPSF}) = p.λ
_arm_wavelength(p::Union{PlanePSF,AnalyzerPSF}) = _arm_wavelength(p.psf)
_arm_tolerance(top, bottom) = sqrt(max(eps(_psf_float(top)), eps(_psf_float(bottom))))
_vector_psf(::AbstractPSF) = nothing
_vector_psf(p::VectorPSF) = p
_vector_psf(p::Union{PlanePSF,AnalyzerPSF}) = _vector_psf(p.psf)

function _check_fourpi_vectors(top, bottom)
    a, b = _vector_psf(top), _vector_psf(bottom)
    for v in (a, b)
        isnothing(v) && continue
        v.nₐ <= v.n_medium || throw(ArgumentError(
            "A 4Pi VectorPSF arm has evanescent components that grow " *
            "at negative emitter z; use n_medium ≥ NA"))
    end
    if !isnothing(a) && !isnothing(b) && length(a.vector_pupils) == 1
        pa = [a.dipole.px, a.dipole.py, -a.dipole.pz]
        pb = [b.dipole.px, b.dipole.py, b.dipole.pz]
        # Same-typed DipoleVector arguments can bypass its normalizing factory.
        na, nb = norm(pa), norm(pb)
        τ = _arm_tolerance(top, bottom)
        if !(isfinite(na) && isfinite(nb) && na > 0 && nb > 0 &&
             isapprox(pa / na, pb / nb; atol=τ, rtol=τ))
            throw(ArgumentError(
                "Fixed 4Pi dipoles must be mirrored; build the bottom arm with " *
                "DipoleVector(px, py, -pz) from the top dipole"))
        end
    end
end

"""
    FourPiPSF(top::AbstractPSF, bottom::AbstractPSF; phase::Real=0)

Interfere two opposed objectives about a shared focal plane.

# Arguments
- `top, bottom`: Complex-field arms with equal dipole and component counts
  and the same wavelength within relative tolerance `sqrt(eps(T))` of the coarser arm type.
  A fixed bottom dipole must mirror the top: `DipoleVector(px, py, -pz)`.
  Vector arms require `n_medium ≥ NA` to avoid evanescent growth at negative z.
- `phase`: Cavity phase in radians.

# Returns
A 3D PSF evaluating the top at `z` and bottom at `-z`, measured from the shared
focal plane. The images overlap with the same orientation: one arm has one
more mirror than the other, so the camera sees both arms' x and y components
in a common frame. For each dipole, fields combine as
`(Et + cis(phase) * s * Eb) / sqrt(2)`. The mirror signs are `(1, 1, -1)` for
rotating x, y, z dipoles and `(1,)` for fixed mirrored dipoles or scalar arms.
Dipole intensities add incoherently. `amplitude` returns a bare scalar for
scalar arms, `[Ex, Ey]` for one vector dipole, and errors for rotating dipoles.

Arms of unit energy give a channel of unit energy averaged over the fringe.
`p(φ) + p(φ + π) == top(x,y,z) + bottom(x,y,-z)` always holds
(up to rounding). Polarization-dependent cavity
phases use one `AnalyzerPSF(FourPiPSF(top, bottom; phase=φ); component=c)` per channel.

# Examples
```julia
arm = ScalarPSF(1.2, 0.68, 1.33)
p = FourPiPSF(arm, arm; phase=π/2)
intensity = p(0.0, 0.0, 0.1)
```
"""
struct FourPiPSF{T, A <: AbstractPSF, B <: AbstractPSF} <: Abstract3DPSF{T}
    top::A
    bottom::B
    phase::T
    function FourPiPSF(top::A, bottom::B; phase::Real=0) where {
        A <: AbstractPSF,B <: AbstractPSF}
        st, sb = _mirror_signs(top), _mirror_signs(bottom)
        length(st) == length(sb) && _n_components(top) == _n_components(bottom) ||
            throw(ArgumentError(
                "4Pi arms differ in dipole count or component count; " *
                "use two scalar arms, two fixed vector arms, or two rotating vector arms"))
        _check_fourpi_vectors(top, bottom)
        T = promote_type(_psf_float(top), _psf_float(bottom))
        λt, λb = _arm_wavelength(top), _arm_wavelength(bottom)
        isapprox(λt, λb; rtol=_arm_tolerance(top, bottom), atol=zero(T)) ||
            throw(ArgumentError(
                "4Pi arms have different wavelengths λ=$λt and λ=$λb; " *
                "use matching wavelengths for one optical frequency"))
        isfinite(phase) && isfinite(T(phase)) || throw(ArgumentError(
            "phase must be finite; supply a finite cavity phase"))
        return new{T,A,B}(top, bottom, T(phase))
    end
end
_n_components(p::FourPiPSF) = _n_components(p.top)
function dipole_fields(p::FourPiPSF, x::Real, y::Real, z::Real)
    Et = dipole_fields(p.top, x, y, z)
    Eb = dipole_fields(p.bottom, x, y, -z)
    s = _mirror_signs(p.bottom)
    return [(a .+ cis(p.phase) .* sign .* b) ./ sqrt(2) for (a, b, sign) in zip(Et, Eb, s)]
end
(p::FourPiPSF)(x::Real, y::Real, z::Real) =
    sum(sum(abs2, E) for E in dipole_fields(p, x, y, z))
amplitude(p::FourPiPSF, x::Real, y::Real, z::Real) =
    _single_dipole_amplitude(dipole_fields(p, x, y, z), Val(_n_components(p)))

_n_dipoles(p::VectorPSF) = length(p.vector_pupils)
_n_dipoles(p::Union{PlanePSF,AnalyzerPSF}) = _n_dipoles(p.psf)
_n_dipoles(p::FourPiPSF) = _n_dipoles(p.top)

function _check_vector_amplitude_integration(psf, support)
    support == Inf || throw(ArgumentError(
        "Vector amplitude integration does not support finite support; " *
        "omit support or use support=Inf"))
    _n_dipoles(psf) > 1 && error(
        "integrate_pixels_amplitude() not meaningful for multiple dipole orientations")
end

"""
    integrate_pixels_amplitude(psf::Union{PlanePSF,AnalyzerPSF,FourPiPSF}, camera, emitter;
                               support=Inf, sampling=2, threaded=true)
    integrate_pixels_amplitude(psf::Union{PlanePSF,AnalyzerPSF,FourPiPSF}, camera, emitters;
                               support=Inf, sampling=2, threaded=true)

Integrate wrapper field amplitudes while preserving vector polarization components.

# Arguments
- `psf`: A plane, analyzer, or 4Pi PSF with a single coherent dipole.
- `camera`: Camera geometry defining pixel edges in microns.
- `emitter` or `emitters`: One emitter or a vector of emitters, summed coherently.
- `support`: Scalar wrappers retain the generic support options. Vector wrappers
  require `Inf`, as for `VectorPSF` integration without a support keyword.
- `sampling`: Subpixel samples per dimension (default 2).
- `threaded`: Enable threaded integration (default true).

# Returns
A complex `[y, x, pol]` array for vector fields, with components `Ex` and `Ey`.
Scalar wrappers retain the generic complex matrix. Several incoherent dipoles
throw the same error as `VectorPSF` amplitude integration.

# Examples
```julia
p = PlanePSF(VectorPSF(1.2, 0.68, DipoleVector(1.0, 0.0, 0.0)); dz=0.1)
camera = IdealCamera(4, 4, 0.1)
fields = integrate_pixels_amplitude(p, camera, Emitter3D(0.2, 0.2, 0.0, 1.0))
```
"""
function integrate_pixels_amplitude(
    psf::Union{PlanePSF,AnalyzerPSF,FourPiPSF}, camera::AbstractCamera,
    emitter::AbstractEmitter;
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf,
    sampling::Integer=2, threaded::Bool=true
)
    if _n_components(psf) != 2
        return invoke(integrate_pixels_amplitude,
            Tuple{AbstractPSF,AbstractCamera,AbstractEmitter}, psf, camera, emitter;
            support=support, sampling=sampling, threaded=threaded)
    end
    _check_vector_amplitude_integration(psf, support)
    return _integrate_pixels_generic(psf, camera.pixel_edges_x, camera.pixel_edges_y,
        emitter, amplitude, Complex{typeof(emitter.photons)}, (2,);
        sampling=sampling, threaded=threaded)
end

function integrate_pixels_amplitude(
    psf::Union{PlanePSF,AnalyzerPSF,FourPiPSF}, camera::AbstractCamera,
    emitters::Vector{<:AbstractEmitter};
    support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf,
    sampling::Integer=2, threaded::Bool=true
)
    if _n_components(psf) != 2
        return invoke(integrate_pixels_amplitude,
            Tuple{AbstractPSF,AbstractCamera,Vector{<:AbstractEmitter}},
            psf, camera, emitters;
            support=support, sampling=sampling, threaded=threaded)
    end
    _check_vector_amplitude_integration(psf, support)
    T = isempty(emitters) ? Complex{Float64} : Complex{typeof(first(emitters).photons)}
    result = zeros(T, length(camera.pixel_edges_y)-1, length(camera.pixel_edges_x)-1, 2)
    buffer = similar(result)
    for emitter in emitters
        fill!(buffer, zero(T))
        _integrate_pixels_generic!(buffer, psf,
            camera.pixel_edges_x, camera.pixel_edges_y, emitter, amplitude;
            sampling=sampling, threaded=threaded)
        result .+= buffer
    end
    return result
end

_carries_split(::AbstractPSF) = false
_carries_split(::AnalyzerPSF) = true
_carries_split(p::PlanePSF) = _carries_split(p.psf)
_carries_split(p::FourPiPSF) = _carries_split(p.top) || _carries_split(p.bottom)

"""
    MultiChannelPSF(; channels, registration=nothing, transmission=nothing,
                    z_offset=nothing, reference=1)

Construct a collection of complete PSFs with a shared emitter coordinate frame.
This container is deliberately not an `AbstractPSF`; one-channel consumers must
handle its channels explicitly.

# Arguments
- `channels`: Nonempty collection of `AbstractPSF`s, each normalized alone.
- `registration`: One 2×3 matrix `[A t]` per channel, mapping absolute reference
  camera coordinates in microns to channel coordinates: `[xᵢ,yᵢ] = A*[x,y]+t`.
  Defaults to identity; absolute coordinates keep registration valid after crops.
- `transmission`: Nonnegative fractions of the emitter's total photons; their
  sum must not exceed one. Defaults to equal splitting.
- `z_offset`: Per-channel emitter z shifts in microns, with the `PlanePSF.dz`
  sign convention. Defaults to zero.
- `reference`: Reference channel index (default 1), whose registration and z
  shift must be identity and zero within `sqrt(eps(T))`.

# Returns
A `MultiChannelPSF{T,P}` with channel element type narrowed when homogeneous.
`T` is the promoted float type of supplied numeric arrays, restricted to
`Float32` or `Float64`; integer and rational inputs promote to `Float64`,
as do all-default inputs. All numeric entries must be finite. Each channel integrates
to one alone (4Pi phase channels on average over the fringe); this cannot be
checked by the constructor. Raw `AnalyzerPSF` channels and wrappers carrying
their photon split are refused: tabulate them into a `SplinePSF` scaled to unit
integral and put the split in `transmission`. See [`channel_position`](@ref)
for the forward law.

# Examples
```julia
p = ScalarPSF(1.2, 0.68, 1.33)
mc = MultiChannelPSF(; channels=[p, PlanePSF(p; dz=0.2)])
position = channel_position(mc, 2, 1.0, 1.0, 0.1)
```
"""
struct MultiChannelPSF{T <: AbstractFloat, P <: AbstractPSF}
    channels::Vector{P}
    registration::Vector{Matrix{T}}
    transmission::Vector{T}
    z_offset::Vector{T}
    reference::Int
    function MultiChannelPSF(; channels, registration=nothing, transmission=nothing,
                             z_offset=nothing, reference=1)
        ch = collect(channels)
        N = length(ch)
        N > 0 || throw(ArgumentError("No channels supplied; provide at least one PSF"))
        all(p -> p isa AbstractPSF, ch) || throw(ArgumentError(
            "Every channel must be an AbstractPSF; replace invalid channels with PSFs"))
        any(_carries_split, ch) && throw(ArgumentError(
            "An AnalyzerPSF carries the photon split itself; " *
            "tabulate it into a SplinePSF scaled to unit integral " *
            "and put the split in transmission"))
        for (name, values) in (
            ("registration", registration), ("transmission", transmission),
            ("z_offset", z_offset))
            isnothing(values) || length(values) == N || throw(ArgumentError(
                "$name length must match channels; provide $N entries"))
        end
        if !isnothing(registration)
            all(r -> r isa AbstractMatrix && size(r) == (2, 3), registration) ||
                throw(ArgumentError(
                    "Each registration must be 2×3; supply matrices [A t] in microns"))
        end
        # Defaults do not widen supplied Float32 (or other float) data.
        numeric_type = nothing
        for (values, is_registration) in (
            (registration, true), (transmission, false), (z_offset, false))
            isnothing(values) && continue
            entries = is_registration ? Iterators.flatten(values) : values
            for v in entries
                v isa Real && isfinite(v) || throw(ArgumentError(
                    "Channel numeric values must be finite real numbers; " *
                    "replace non-finite or non-real values"))
                V = v isa Union{Integer,Rational} ? Float64 : typeof(float(v))
                numeric_type = isnothing(numeric_type) ? V : promote_type(numeric_type, V)
            end
        end
        T = isnothing(numeric_type) ? Float64 : numeric_type
        T in (Float32, Float64) || throw(ArgumentError(
            "Channel numeric type $T cannot be saved; pass Float32 or Float64 values"))
        reg = isnothing(registration) ?
              [T[1 0 0; 0 1 0] for _ in 1:N] : [Matrix{T}(r) for r in registration]
        tr = isnothing(transmission) ? fill(one(T) / N, N) : Vector{T}(transmission)
        dz = isnothing(z_offset) ? zeros(T, N) : Vector{T}(z_offset)
        all(r -> all(isfinite, r), reg) && all(isfinite, tr) && all(isfinite, dz) ||
            throw(ArgumentError(
                "Channel numeric values must remain finite after conversion; " *
                "use representable finite values"))
        all(t -> t >= 0, tr) || throw(ArgumentError(
            "transmission must be nonnegative; supply photon fractions ≥ 0"))
        τ = sqrt(eps(T))
        sum(tr) <= one(T) + τ || throw(ArgumentError(
            "transmission sums to more than one; supply fractions of one total " *
            "photon budget, and normalize analyzer channels alone"))
        reference isa Integer && 1 <= reference <= N || throw(ArgumentError(
            "reference must be a channel index in 1:$N; select an existing channel"))
        all(abs.(reg[reference] .- T[1 0 0; 0 1 0]) .<= τ) && abs(dz[reference]) <= τ ||
            throw(ArgumentError(
                "The reference channel must have identity registration and zero " *
                "z_offset; express all transforms relative to it"))
        P = reduce(typejoin, typeof.(ch))
        return new{T,P}(Vector{P}(ch), reg, tr, dz, Int(reference))
    end
end

"""
    channel_position(mc::MultiChannelPSF, i, x, y, z) -> (xᵢ, yᵢ, zᵢ)

Map an emitter from the reference camera frame into channel `i`.

# Arguments
- `mc`: Multi-channel PSF container.
- `i`: Channel index.
- `x, y, z`: Reference-frame emitter position in microns, using absolute camera
  coordinates as in `SMLMData` pixel edges (not pixel indices).

# Returns
`(xᵢ, yᵢ, zᵢ)` with `[xᵢ,yᵢ] = A*[x,y]+t` from `registration[i]` and
`zᵢ = z + z_offset[i]`. The defining forward law is

    expected image i = N * mc.transmission[i] *
        integrate_pixels(mc.channels[i], camera_i, emitter at channel_position)

where `N` is the emitter's total photons and the integrated emitter has unit
photons. Each channel must integrate to one alone, with 4Pi phase channels
normalized on average over the fringe. The constructor cannot check arbitrary
channel normalization; it refuses raw `AnalyzerPSF` channels and wrappers carrying
their photon split. Tabulate these into a `SplinePSF` scaled to unit integral
and put the split in `transmission`.

# Examples
```julia
p = ScalarPSF(1.2, 0.68, 1.33)
mc = MultiChannelPSF(; channels=[p, p], z_offset=[0.0, 0.2])
channel_position(mc, 2, 1.0, 2.0, 0.1)  # (1.0, 2.0, 0.3), up to rounding
```
"""
function channel_position(mc::MultiChannelPSF, i::Integer, x::Real, y::Real, z::Real)
    1 <= i <= length(mc.channels) || throw(ArgumentError(
        "Channel index must be in 1:$(length(mc.channels)); select an existing channel"))
    r = mc.registration[i]
    return (r[1,1] * x + r[1,2] * y + r[1,3],
            r[2,1] * x + r[2,2] * y + r[2,3], z + mc.z_offset[i])
end

Base.show(io::IO, p::PlanePSF) = print(io, "PlanePSF(", p.psf, "; dz=", p.dz, ")")
Base.show(io::IO, p::AnalyzerPSF) =
    print(io, "AnalyzerPSF(", p.psf, "; component=", p.component, ")")
Base.show(io::IO, p::FourPiPSF) =
    print(io, "FourPiPSF(", p.top, ", ", p.bottom, "; phase=", p.phase, ")")
Base.show(io::IO, mc::MultiChannelPSF) =
    print(io, "MultiChannelPSF(", length(mc.channels),
          " channels; reference=", mc.reference, ")")
