using MicroscopePSFs
using Documenter

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
        "Tutorials" => "Tutorials.md",
        "PSF Types" => "psftypes.md",
        "Aberrations" => "Aberrations.md",
        "Interpolation" => "Interpolation.md",
        "Library" => "Library.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSMLM/MicroscopePSFs.jl",
    devbranch="main",
)
