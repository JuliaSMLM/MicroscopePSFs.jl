using Documenter
using MicroscopePSFs

DocMeta.setdocmeta!(MicroscopePSFs, :DocTestSetup, :(using MicroscopePSFs); recursive=true)

makedocs(;
    modules=[MicroscopePSFs],
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
    doctest = false,  # Disable doctests for now
    checkdocs = :none,  # Don't throw errors for missing docstrings
    warnonly = [:missing_docs, :docs_block, :autodocs_block, :cross_references],  # Only warn for these error types
)

# Uncomment this when you're ready to deploy
deploydocs(;
    repo="github.com/JuliaSMLM/MicroscopePSFs.jl",
    devbranch="main",
)
