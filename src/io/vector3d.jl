# src/io/vector3d.jl

"""
    _save_psf_impl(file::HDF5.File, psf::VectorPSF)

Save a VectorPSF to an HDF5 file, including all optical parameters,
dipole orientation, and vector pupil fields.

# Arguments
- `file`: Open HDF5 file handle
- `psf`: VectorPSF to save
"""
function _save_psf_impl(file::HDF5.File, psf::VectorPSF)
    # Save optical parameters
    params = create_group(file, "parameters")
    params["na"] = psf.nₐ
    params["lambda"] = psf.λ
    params["n_medium"] = psf.n_medium
    params["n_coverslip"] = psf.n_coverslip
    params["n_immersion"] = psf.n_immersion
    params["z_stage"] = psf.z_stage
    params["num_pupils"] = length(psf.vector_pupils)
    
    # Save Zernike coefficients if available
    if !isnothing(psf.zernike_coeffs)
        _io_save_zernike_coeffs(file, psf.zernike_coeffs)
    end
    
    # Save dipole orientation
    dipole = create_group(file, "dipole")
    dipole["px"] = psf.dipole.px
    dipole["py"] = psf.dipole.py
    dipole["pz"] = psf.dipole.pz
    
    # Save pupil functions
    data = create_group(file, "data")
    
    # Save each pupil
    for (i, pupil) in enumerate(psf.vector_pupils)
        pupil_group = create_group(data, "pupil_$i")
        
        # Save Ex field
        pupil_ex = create_group(pupil_group, "pupil_ex")
        _io_save_complex_array(pupil_ex, "field", pupil.Ex.field)
        
        # Save Ey field
        pupil_ey = create_group(pupil_group, "pupil_ey")
        _io_save_complex_array(pupil_ey, "field", pupil.Ey.field)
    end
    
    # Save base pupil if available
    if !isnothing(psf.base_pupil)
        base_pupil = create_group(data, "base_pupil")
        _io_save_pupil_params(base_pupil, psf.base_pupil)
        _io_save_complex_array(base_pupil, "field", psf.base_pupil.field)
    end
end

"""
    _load_psf_impl(file::HDF5.File, ::Type{VectorPSF}) 

Load a VectorPSF from an HDF5 file, reconstructing all components.

# Arguments
- `file`: Open HDF5 file handle
- `::Type{VectorPSF}`: Type to load

# Returns
- VectorPSF reconstructed from stored data
"""
function _load_psf_impl(file::HDF5.File, ::Type{VectorPSF}) 
    # Load optical parameters
    params = file["parameters"]
    _io_check_required_fields(params, ["na", "lambda", "n_medium",
                                     "n_coverslip", "n_immersion", "z_stage", "num_pupils"])
    
    nₐ = read(params["na"])
    λ = read(params["lambda"])
    n_medium = read(params["n_medium"])
    n_coverslip = read(params["n_coverslip"])
    n_immersion = read(params["n_immersion"])
    z_stage = read(params["z_stage"])
    num_pupils = read(params["num_pupils"])
    
    # Load dipole orientation
    dipole_group = file["dipole"]
    _io_check_required_fields(dipole_group, ["px", "py", "pz"])
    
    px = read(dipole_group["px"])
    py = read(dipole_group["py"])
    pz = read(dipole_group["pz"])
    dipole = DipoleVector(px, py, pz)
    
    # Load Zernike coefficients if available
    zernike_coeffs = _io_load_zernike_coeffs(file)
    
    # Load pupil fields
    data = file["data"]
    
    # Create vector to hold all pupils
    vector_pupils = Vector{VectorPupilFunction{typeof(nₐ)}}(undef, num_pupils)
    
    # Load each pupil
    for i in 1:num_pupils
        pupil_group = data["pupil_$i"]
        
        # Load Ex field
        pupil_ex = pupil_group["pupil_ex"]
        ex_field = _io_load_complex_array(pupil_ex, "field")
        
        # Load Ey field
        pupil_ey = pupil_group["pupil_ey"]
        ey_field = _io_load_complex_array(pupil_ey, "field")
        
        # Create pupil functions
        Ex = PupilFunction(nₐ, λ, n_medium, ex_field)
        Ey = PupilFunction(nₐ, λ, n_medium, ey_field)
        vector_pupils[i] = VectorPupilFunction(nₐ, λ, n_medium, n_coverslip, n_immersion, Ex, Ey)
    end
    
    # Load base pupil if available
    base_pupil = nothing
    if haskey(data, "base_pupil")
        bp_group = data["base_pupil"]
        bp_field = _io_load_complex_array(bp_group, "field")
        base_pupil = PupilFunction(nₐ, λ, n_medium, bp_field)
    end
    
    # Create the VectorPSF
    return VectorPSF{typeof(nₐ)}(
        nₐ, λ, n_medium, n_coverslip, n_immersion,
        dipole, z_stage, vector_pupils, base_pupil, zernike_coeffs
    )
end