# Recursive wrapper and channel storage; loaders enforce constructor invariants.
function _save_psf_impl(file::_H5Parent, p::PlanePSF)
    params = create_group(file, "parameters")
    params["dz"] = p.dz
    _io_save_child(file, "psf", p.psf)
end
function _load_psf_impl(file::_H5Parent, ::Type{PlanePSF})
    return PlanePSF(_io_load_child(file, "psf"); dz=read(file["parameters/dz"]))
end

function _save_psf_impl(file::_H5Parent, p::AnalyzerPSF)
    params = create_group(file, "parameters")
    params["component"] = p.component
    _io_save_child(file, "psf", p.psf)
end
function _load_psf_impl(file::_H5Parent, ::Type{AnalyzerPSF})
    return AnalyzerPSF(_io_load_child(file, "psf");
                       component=read(file["parameters/component"]))
end

function _save_psf_impl(file::_H5Parent, p::FourPiPSF)
    params = create_group(file, "parameters")
    params["phase"] = p.phase
    _io_save_child(file, "top", p.top)
    _io_save_child(file, "bottom", p.bottom)
end
function _load_psf_impl(file::_H5Parent, ::Type{FourPiPSF})
    return FourPiPSF(_io_load_child(file, "top"), _io_load_child(file, "bottom"),
        read(file["parameters/phase"]), Val(false))
end

function _save_psf_impl(file::_H5Parent, mc::MultiChannelPSF)
    params = create_group(file, "parameters")
    params["n_channels"] = length(mc.channels)
    params["reference"] = mc.reference
    params["transmission"] = mc.transmission
    params["z_offset"] = mc.z_offset
    # Explicit allocation preserves Julia size (2, 3, N), including N == 1.
    registration = Array{eltype(mc.transmission)}(undef, 2, 3, length(mc.channels))
    for i in eachindex(mc.channels)
        registration[:, :, i] = mc.registration[i]
    end
    params["registration"] = registration
    channels = create_group(file, "channels")
    for (i, p) in enumerate(mc.channels)
        _io_save_child(channels, "channel_$i", p)
    end
end
function _load_psf_impl(file::_H5Parent, ::Type{MultiChannelPSF})
    params = file["parameters"]
    N = read(params["n_channels"])
    registration = read(params["registration"])
    size(registration) == (2, 3, N) || throw(ArgumentError(
        "Stored registration must have size (2, 3, n_channels); " *
        "repair the registration dataset"))
    return MultiChannelPSF(;
        channels=[_io_load_child(file["channels"], "channel_$i") for i in 1:N],
        registration=[registration[:, :, i] for i in 1:N],
        transmission=read(params["transmission"]),
        z_offset=read(params["z_offset"]), reference=read(params["reference"]))
end
