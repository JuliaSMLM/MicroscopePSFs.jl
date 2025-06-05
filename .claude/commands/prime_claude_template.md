# Development of Package MyPackage.jl 

We are actively developing this Julia package. 

This project-specific overview should give Claude context when starting a new code session. Be as specific as possible about your project's purpose, architecture, and requirements.

## High Level Custom Instructions 
*Examples: Always include docstrings, use PascalCase for types, prefer functions over macros, etc.*

## Conceptual Overview 
*Example: This package implements several numerical optimization algorithms for high-dimensional data, with a focus on performance and memory efficiency.*

## Dependencies
*Example: Uses LinearAlgebra for matrix operations, StaticArrays for performance-critical code, and Plots for visualization support.*

## User Interfaces and Exports
*Example: Exposes three main functions: optimize(), gradient_descent(), and plot_convergence(). All functions follow a consistent parameter order.*

## Code Conventions
*Example: Function names use snake_case, types use PascalCase. Error messages are specific and actionable. Core algorithms are documented with mathematical explanations in docstrings.*

## Testing Standards
*Example: Unit tests for all public functions, integration tests for workflows, benchmarks for performance-critical sections. Test coverage target is >90%.*

## File Organization/Architecture
*Example: src/algorithms/ contains implementations, src/interfaces/ has API definitions, src/visualization/ handles plotting functionality.*

## Priorities and Preferences 
*Example: Prioritize correctness over performance, except in core algorithm implementations. Prefer immutable data structures where appropriate.*

*Note: Delete any sections that aren't relevant to your package, and add any project-specific sections needed.*