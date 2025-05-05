# Overview Review Command

This command helps review the api_overview.md file for accuracy against the source code.

## Usage
```
/overview_review
```

## What This Does

When you run this command, it will:

1. Perform a detailed verification of the llm_api_overview.md document against the actual source code implementation
2. Check all type definitions, function signatures, and constructor patterns
3. Verify exact constructor signatures (including positional vs keyword parameters)
4. Compare code examples against source implementations line-by-line
5. Check for unit discrepancies between documentation and code
6. Verify all available constructor forms, not just simplified usage examples
7. Create a comprehensive checklist of discrepancies that need to be fixed

## Focus Areas

The review will specifically examine:

- Type hierarchies and inheritance relationships
- Constructor signatures (including positional vs keyword arguments)
- Default parameter values and their consistency
- Parameter names and naming conventions (e.g., Unicode characters vs ASCII)
- Units used in examples vs implementation (e.g., nanometers vs microns)
- Missing or extraneous parameters in examples
- Completeness of API coverage
- Code style consistency

## Output Format

The command will produce:

1. A comprehensive checklist of all examples in the overview
2. Each item marked as correct or incorrect
3. For incorrect items, specific details about the discrepancy
4. Suggested corrections based on the actual implementation
5. Priority levels for fixing each issue