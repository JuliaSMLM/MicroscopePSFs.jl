"""
IO functionality for PSFs including serialization and caching
"""

"""
    save_psf(filename::String, psf::AbstractPSF; format=:hdf5)

Save a PSF to file. Supports multiple formats:
- HDF5 (.h5, .hdf5)
- JSON (.json)
"""
function save_psf(filename::String, psf::AbstractPSF; format=:hdf5)
    error("save_psf not implemented for $(typeof(psf))")
end

"""
    load_psf(filename::String) -> AbstractPSF

Load a PSF from file. Format is determined by file extension.
"""
function load_psf(filename::String)
    ext = lowercase(splitext(filename)[2])
    if ext in [".h5", ".hdf5"]
        return load_psf_hdf5(filename)
    elseif ext == ".json"
        return load_psf_json(filename)
    else
        error("Unsupported file format: $ext")
    end
end

# Cache utilities for precomputed fields

"""
    cache_field(
        filename::String,
        psf::AbstractPSF,
        roi::ROI;  # Region of interest
        sampling_nm::Real=1.0
    )

Cache precomputed complex field for faster evaluation.
"""
function cache_field(filename::String, psf::AbstractPSF, roi::ROI; sampling_nm::Real=1.0)
    error("cache_field not implemented for $(typeof(psf))")
end

"""
    load_cached_field(filename::String) -> Tuple{AbstractPSF, Array{Complex{Float64}}}

Load a cached field and associated PSF parameters.
"""
function load_cached_field(filename::String)
    error("Loading cached field not implemented")
end

# Implement format-specific IO
function save_psf_hdf5(filename::String, psf::AbstractPSF)
    error("HDF5 saving not implemented for $(typeof(psf))")
end

function load_psf_hdf5(filename::String)
    error("HDF5 loading not implemented")
end

function save_psf_json(filename::String, psf::AbstractPSF)
    error("JSON saving not implemented for $(typeof(psf))")
end

function load_psf_json(filename::String)
    error("JSON loading not implemented")
end
