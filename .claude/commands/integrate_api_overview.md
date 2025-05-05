# Integrating API Overview into Package

This document provides instructions for integrating the API overview into the package structure.

## Integration Steps

### Step 1: Create src/api.jl

Create a new file at `src/api.jl` with the following generic code:

```julia
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
```

### Step 2: Update Your Main Package File

Add the following to your main package file (`src/MicroscopePSFs.jl`):

1. Include the `api.jl` file near the end of your module:

```julia
# Include the API overview functionality
include("api.jl")
```

2. Update the module docstring to mention this documentation feature.

Your module docstring should look something like this:

```julia
"""
    MicroscopePSFs

A Julia package for working with microscope Point Spread Functions (PSFs).
This package provides implementations of common PSF models and tools for 
integrating them with camera geometry for single-molecule localization 
microscopy applications.

# API Overview
For a comprehensive overview of the API, use the help mode on `api_overview`:

    ?api_overview

Or access the complete API documentation programmatically:

    docs = MicroscopePSFs.api_overview()
"""
module MicroscopePSFs
# Rest of your module code...
```

### Step 3: Test The Implementation

After implementing these changes:

1. Build/load your package
2. Test the help mode: `?MicroscopePSFs.api_overview`
3. Test programmatic access: `docs = MicroscopePSFs.api_overview()`

## Why Integrate the API Overview?

### For Humans
- Provides a **concise reference** without diving into full documentation
- Offers **quick-start examples** for common use cases
- Shows **relevant patterns** more clearly than individual docstrings
- Creates an **at-a-glance understanding** of package capabilities
- Complements detailed documentation with a **higher-level view**

### For AI Assistants
- Enables **better code generation** with correct API patterns
- Provides **structured context** about type hierarchies and relationships
- Offers **consistent examples** to learn from when generating code
- Creates **clear patterns** to follow when suggesting package usage
- Helps avoid **common pitfalls** or misunderstandings about the API