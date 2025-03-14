# MicroscopePSFs Development Guide

## Build & Test Commands
- Run all tests: `julia --project -e "using Pkg; Pkg.test()"`
- Run specific test: Create a custom script with specific testset, e.g.:
```julia
using Test, MicroscopePSFs
@testset "MyTest" begin
    # Test code here
end
```
- Build docs: `julia --project=docs/ docs/make.jl`

## Style Guidelines
- **Naming**: Types use PascalCase; functions/variables use snake_case; mutating functions use bang suffix (!)
- **Typing**: Use strong typing for function parameters and struct fields
- **Documentation**: Write docstrings for all public functions with description, arguments, returns, and examples
- **Imports**: Group by standard library followed by third-party packages
- **Error handling**: Use descriptive error messages with `throw(ArgumentError("message"))`
- **Physics variables**: Use Unicode math symbols (λ, nₐ) for readability when appropriate
- **Testing**: Write tests that verify mathematical correctness using `isapprox` for floating-point comparisons

## Visualization
- Use CairoMakie for examples in dev folder and documentation, but don't make it a dependency in the main package