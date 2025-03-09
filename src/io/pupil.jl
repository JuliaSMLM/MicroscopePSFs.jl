# src/io/pupil.jl

"""
I/O implementations for pupil function types (PupilFunction and VectorPupilFunction).
These represent the complex field in the pupil plane of the microscope.
"""

# ---- PupilFunction ----

"""
    _save_psf_impl(file::HDF5.File, pupil::PupilFunction)

Save a PupilFunction to an HDF5 file, including physical parameters and complex field.
"""
function _save_psf_impl(file::HDF5.File, pupil::PupilFunction)
    # Save basic parameters
    params = create_group(file, "parameters")
    params["na"] = pupil.nₐ
    params["lambda"] = pupil.λ
    params["n"] = pupil.n
    
    # Save grid dimensions
    dims = size(pupil.field)
    params["grid_size_y"] = dims[1]
    params["grid_size_x"] = dims[2]
    
    # Save field data
    data = create_group(file, "data")
    _io_save_complex_array(data, "field", pupil.field)
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{PupilFunction{T}}) where T

Load a PupilFunction from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{PupilFunction{T}}) where T
    # Load parameters
    params = file["parameters"]
    _io_check_required_fields(params, ["na", "lambda", "n"])
    
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n = read(params["n"])
    
    # Load field data
    data = file["data"]
    field = _io_load_complex_array(data, "field")
    
    # Create and return the PupilFunction
    return PupilFunction(nₐ, λ, n, field)
end

# ---- VectorPupilFunction ----

"""
    _save_psf_impl(file::HDF5.File, pupil::VectorPupilFunction)

Save a VectorPupilFunction to an HDF5 file, including all physical parameters
and both Ex and Ey field components.
"""
function _save_psf_impl(file::HDF5.File, pupil::VectorPupilFunction)
    # Save parameters
    params = create_group(file, "parameters")
    params["na"] = pupil.nₐ
    params["lambda"] = pupil.λ
    params["n_medium"] = pupil.n_medium
    params["n_coverslip"] = pupil.n_coverslip
    params["n_immersion"] = pupil.n_immersion
    
    # Save grid dimensions
    dims = size(pupil.Ex.field)
    params["grid_size_y"] = dims[1]
    params["grid_size_x"] = dims[2]
    
    # Save Ex field data
    ex_data = create_group(file, "ex_field")
    _io_save_complex_array(ex_data, "field", pupil.Ex.field)
    
    # Save Ey field data
    ey_data = create_group(file, "ey_field")
    _io_save_complex_array(ey_data, "field", pupil.Ey.field)
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{VectorPupilFunction{T}}) where T

Load a VectorPupilFunction from an HDF5 file.
"""
function _load_psf_impl(file::HDF5.File, ::Type{VectorPupilFunction{T}}) where T
    # Load parameters
    params = file["parameters"]
    _io_check_required_fields(params, ["na", "lambda", "n_medium", 
                                     "n_coverslip", "n_immersion"])
    
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n_medium = read(params["n_medium"])
    n_coverslip = read(params["n_coverslip"])
    n_immersion = read(params["n_immersion"])
    
    # Load Ex field
    ex_data = file["ex_field"]
    ex_field = _io_load_complex_array(ex_data, "field")
    
    # Load Ey field
    ey_data = file["ey_field"]
    ey_field = _io_load_complex_array(ey_data, "field")
    
    # Create PupilFunction objects for Ex and Ey
    Ex = PupilFunction(nₐ, λ, n_medium, ex_field)
    Ey = PupilFunction(nₐ, λ, n_medium, ey_field)
    
    # Create and return the VectorPupilFunction
    return VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, Ex, Ey)
end