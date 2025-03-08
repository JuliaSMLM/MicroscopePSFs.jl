using Documenter

# Temporarily disable using MicroscopePSFs due to package errors
# using MicroscopePSFs
# DocMeta.setdocmeta!(MicroscopePSFs, :DocTestSetup, :(using MicroscopePSFs); recursive=true)

makedocs(;
    # modules=[MicroscopePSFs], # Temporarily disabled
    authors="klidke@unm.edu",
    repo="https://github.com/JuliaSMLM/MicroscopePSFs.jl/blob/{commit}{path}#{line}",
    sitename="MicroscopePSFs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSMLM.github.io/MicroscopePSFs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Interface" => "interface.md",
        "Conventions" => "conventions.md",
        "PSF Types" => [
            "Overview" => "psfs/overview.md",
            "Gaussian2D" => "psfs/gaussian2d.md",
            "Airy2D" => "psfs/airy2d.md",
            "Scalar3D" => "psfs/scalar3d.md",
            "Vector3D" => "psfs/vector3d.md",
            "Spline PSF" => "psfs/spline_psf.md",
        ],
        "API Reference" => "api.md",
    ],
    # Disable doctests due to package errors
    doctest = false,
)

# Temporarily disable deploydocs due to package errors
# deploydocs(;
#     repo="github.com/JuliaSMLM/MicroscopePSFs.jl",
#     devbranch="main",
# )