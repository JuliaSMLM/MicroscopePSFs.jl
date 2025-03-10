# Run this script to build the documentation locally.

using Pkg
Pkg.activate("docs")

#Extract package name from pwd()
pkgname = splitpath(pwd())[end]

try
    @info "trying to make docs"
    include("../docs/make.jl")
catch
    @info "trying again with setting up packages for docs"
    @info "this may take a while"
    sleep(5)
    #check if the package is already installed
    try
        Pkg.rm(pkgname)
    catch
    end
    devpath = pwd()
    Pkg.develop(path=devpath)
    Pkg.instantiate()

    @info "making docs"
    include("../docs/make.jl")

end

function open_html_file(filepath::String)
    if Sys.iswindows()
        run(`cmd /c start $filepath`)
    elseif Sys.isapple()
        run(`open $filepath`)
    else # assuming it's a Unix-like OS otherwise
        run(`xdg-open $filepath`)
    end
end

@info "opening docs"
html_file_path = normpath(joinpath(pwd(), "docs/build/index.html"))
open_html_file(html_file_path)