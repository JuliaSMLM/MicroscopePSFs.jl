# Create a prime_claude.md Context File for This Package

Please help me create a prime_claude.md file for this Julia package. This file will provide context about the package's purpose, architecture, and requirements for future Claude sessions.

## Initial Setup

1. First, check if prime_claude.md already exists by reading it
2. If it exists, summarize the current content and ask:
   - "Would you like to continue from where you left off?"
   - "Is there a specific section you'd like to revise?"
   - "Would you prefer to start fresh?"
3. If it doesn't exist, offer to create it using the prime_claude_template.md as a starting point
4. Assume the user has already created their Julia package using PkgTemplates
5. Remember the slash commands (including this one) are installed in their package

## Interactive Approach

Throughout this process:
- Draw on your knowledge of Julia packages and best practices
- Use available tools to explore the user's package structure
- Make specific suggestions based on observed code conventions
- Ask if they want to proceed to the next section before moving on
- Provide examples from popular Julia packages when relevant

## Information Gathering Process

Guide the user through each section systematically:

### Project Overview
- Determine the package name from Project.toml or src/ directory files
- "Can you describe your package's primary functionality in 1-2 sentences?"
- "Who are the intended users of this package?"
- Ask if they want to proceed to the next section

### High Level Instructions
- "Are there any global coding practices everyone should follow?"
- "Do you have specific documentation or testing requirements?"
- "Are there particular Julia idioms you prefer or avoid?"
- Ask if they want to proceed to the next section

### Conceptual Overview
- "What are the core concepts or algorithms your package implements?"
- "Are there mathematical or computational models that underpin the code?"
- "What problem does this package solve for users?"
- Ask if they want to proceed to the next section

### Dependencies
- "What are the main Julia packages your code depends on?"
- "Are there any version constraints or compatibility issues to be aware of?"
- "Do certain dependencies handle critical functionality?"
- Ask if they want to proceed to the next section

### User Interfaces and Exports
- "What functions/types will users interact with most?"
- "Is there a consistent API design pattern across your package?"
- "Are there grouped functionalities users should know about?"
- Ask if they want to proceed to the next section

### Code Conventions
- "Do you have any specific code conventions that differ from standard Julia practices?"
- "If not, I'll assume standard Julia conventions (snake_case for functions/variables, PascalCase for types, etc.)"
- "Any specific documentation style preferences beyond standard Julia docstrings?"
- Ask if they want to proceed to the next section

### Testing Standards
- "Do you have any specific testing approaches beyond standard Julia testing practices?"
- "If not, I'll assume you're using the standard Julia Test package with typical testing patterns"
- "Any particular test coverage goals or benchmarking requirements?"
- Ask if they want to proceed to the next section

### File Organization
- "Do you have any special file organization needs beyond the standard Julia package structure?"
- "If not, I'll assume you're following standard Julia package organization conventions"
- "Any specific modules or directories you'd like to highlight?"
- Ask if they want to proceed to the next section

### Priorities and Preferences
- "What do you prioritize: performance, readability, extensibility, etc.?"
- "Are there specific Julia patterns or features you prefer to use or avoid?"
- "What trade-offs are acceptable in this codebase?"
- Ask if they want to proceed to the final step

## Synthesis and Refinement

1. After gathering information for each section, summarize your understanding
2. Present a draft of each section for confirmation or refinement
3. For sparse areas, suggest typical practices from Julia community standards
4. Remove any sections the user indicates aren't relevant
5. Ask if they'd like to make any final adjustments before finalizing

## Output Generation

1. Format the complete file using Markdown with clear section headings
2. Replace any example text with the user's actual preferences
3. Present the final document and create it in the project's .claude/commands directory:
   - Use the Write tool to create or update `.claude/commands/prime_claude.md`
   - If the `.claude/commands` directory doesn't exist, suggest creating it with `mkdir -p .claude/commands`
   - Remind users to keep this file separate from prime_claude_template.md
4. Explain that they can use this with `/claude prime_claude.md` at the start of sessions
5. Emphasize that prime_claude.md can be updated over time as the project evolves

## If User is Unsure

If the user is uncertain about certain aspects:
- Analyze their existing code to suggest appropriate conventions
- Suggest common Julia community practices
- Provide examples of what other packages typically do
- Offer to revisit sections later if they want to think about it

## Multi-Session Support

1. Acknowledge that building this file might take multiple sessions
2. At the end of each session, summarize progress and what remains to be done
3. Suggest saving the current state of prime_claude.md
4. When starting a new session, read the existing file and suggest where to continue

The goal is to create a prime_claude.md file that will effectively guide future Claude sessions with context specific to the user's project and preferences.