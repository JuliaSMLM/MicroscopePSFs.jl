# Simple API overview functionality for Julia packages

# Determine the package root directory
function _package_root()
    # Get the directory of the current file (api.jl)
    src_dir = dirname(@__FILE__)
    # Go up one level to the package root
    return abspath(joinpath(src_dir, ".."))
end

# Path to the api_overview.md file
const _API_OVERVIEW_PATH = joinpath(_package_root(), "api_overview.md")

# Load the content of the api_overview.md file if it exists
const _API_OVERVIEW = if isfile(_API_OVERVIEW_PATH)
    read(_API_OVERVIEW_PATH, String)
else
    """
    API overview documentation not found.
    
    Expected file: $(basename(_API_OVERVIEW_PATH))
    Expected location: $(dirname(_API_OVERVIEW_PATH))
    """
end

"""
$(_API_OVERVIEW)

---
`api_overview()` returns this documentation as a plain `String`.
"""
function api_overview()
    return _API_OVERVIEW
end

# Note: No export statement - this function remains internal to the package