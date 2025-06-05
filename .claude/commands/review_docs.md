Please perform a verification of the documentation for this Julia package against its implementation. Access the package files to ensure all examples, function signatures, types, and other documented elements accurately reflect the current state of the code.

## Verification Tasks

1. Identify the package structure:
   - Examine the `src/` directory for source code
   - Review the `docs/` directory for documentation files
   - Locate the `README.md` file at the project root

2. Extract all exported items:
   - Parse `src/Package.jl` to identify all explicitly exported names
   - Check submodules for additional exports
   - Create a comprehensive list of all public API elements

3. Verify README examples:
   - Extract all code blocks from README.md
   - Check that all referenced functions and types exist in exports
   - Verify function signatures match actual implementation
   - Confirm constructor patterns match actual type definitions
   - Flag any discrepancies for correction

4. Verify documentation structure:
   - Review `docs/make.jl` for correct module name and page structure
   - Confirm that all pages referenced in `make.jl` exist in `docs/src/`
   - Check that DocTest setup is properly configured

5. Verify API documentation:
   - Compare documented functions, types, and constants against exports
   - Identify items that are exported but not documented
   - Find items that are documented but not exported
   - For each documented function, verify parameter names, types, and descriptions match implementation
   - Check return value descriptions against actual return types

6. Verify examples in documentation:
   - Test all code blocks marked with `@example` or equivalent
   - Ensure all imported modules are listed in dependencies
   - Verify function calls against actual implementations
   - Check type constructors match current definitions
   - Confirm parameter names match implementation

7. Verify type hierarchies:
   - Compare documented type hierarchies against actual implementation
   - Check supertypes and subtypes for consistency
   - Verify that type field descriptions match implementation

8. Verify mathematical formulas:
   - Compare documented formulas with corresponding implementations
   - Check that formula symbols correspond to variable names
   - Verify formula logic against implementation

9. Check docstrings:
   - Extract docstrings from source code
   - Compare against external documentation
   - Verify parameter and return value descriptions
   - Test examples in docstrings against current implementation

10. Run documented examples:
    - Attempt to execute example code when possible
    - Verify outputs match what's described in documentation

## Reporting Format

Please compile a report of all findings, including:

1. A summary of the verification process
2. A list of discrepancies categorized by type:
   - Missing documentation for exported items
   - Documented items that don't exist in code
   - Function signature mismatches
   - Parameter inconsistencies
   - Type hierarchy inconsistencies
   - Broken or outdated examples
   - Formula inconsistencies

3. Recommendations for documentation improvements, including:
   - Items needing documentation
   - Outdated documentation requiring updates
   - Examples needing revision

When assessing discrepancies, please use your judgment to determine the likely cause (code updates, planned features, deprecated functionality, or simple errors).

Thank you for conducting this thorough documentation verification to ensure the package's documentation remains accurate and helpful to users.