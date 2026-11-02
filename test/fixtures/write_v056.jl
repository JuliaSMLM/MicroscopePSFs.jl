# Run in a separate Julia process, with an unmodified git archive of tag v0.5.6:
# julia --project=. test/fixtures/write_v056.jl ARCHIVED_TAG_DIRECTORY test/fixtures
# Loading the archived module directly uses the active project's installed dependencies.
include(joinpath(abspath(ARGS[1]), "src", "MicroscopePSFs.jl"))
using .MicroscopePSFs

for (name, orientation) in (("tilted", (0.6, 0.0, 0.8)),
    ("bottom", (0.6, 0.0, -0.8)), ("transverse", (1.0, 0.0, 0.0)),
    ("rotating", nothing))
    dipole = isnothing(orientation) ? () : (DipoleVector(orientation...),)
    psf = VectorPSF(1.2, 0.68, dipole...; n_medium=1.33, n_coverslip=1.33,
        n_immersion=1.33, grid_size=16, base_zernike=ZernikeCoefficients(6))
    save_psf(joinpath(ARGS[2], "vector-v056-$name.h5"), psf)
    if name == "tilted"
        save_psf(joinpath(ARGS[2], "vector-v056-standalone.h5"), only(psf.vector_pupils))
    end
end
