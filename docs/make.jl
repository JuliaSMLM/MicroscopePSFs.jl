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
            "GaussianPSF" => "psfs/gaussianpsf.md",
            "AiryPSF" => "psfs/airypsf.md",
            "ScalarPSF" => "psfs/scalarpsf.md",
            "VectorPSF" => "psfs/vectorpsf.md",
            "SplinePSF" => "psfs/splinepsf.md",
        ],
        "Examples" => "examples.md",
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
