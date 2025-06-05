# Create a Mathematical and Algorithmic Reference for This Package

Please create a comprehensive mathematical reference document for this scientific/computational package. Analyze the codebase to identify the mathematical foundations, algorithms, and physical principles implemented in this code.

## Initial Setup

1. First, check if math_reference.md already exists by reading it
2. If it exists, summarize the current content and ask:
   - "Would you like to continue from where you left off?"
   - "Is there a specific section you'd like to revise?"
   - "Would you prefer to start fresh?"
3. If it doesn't exist, offer to create it using the math_reference_template.md as a starting point
4. Understand the scientific domain of the package (optics, quantum physics, Bayesian statistics, etc.)
5. Determine the level of mathematical detail appropriate for the package

## Interactive Approach

Throughout this process:
- Draw on your knowledge of scientific computing, physics, mathematics, and statistics
- Use LaTex notation for all mathematical expressions
- Request specific equations from the user when needed
- Ask for clarification on specialized terms or concepts
- Provide examples of well-documented mathematical foundations from similar packages

## Information Gathering Process

Guide the user through each section systematically:

### Package Identification
- "Could you describe the scientific domain of your package?"
- "What is the primary mathematical or physical framework you're implementing?"
- "Who are the intended users of this mathematical reference?"

### Mathematical Foundations
- "What are the core equations that your implementation is based on?"
- "Are there specific coordinate systems or conventions users should know about?"
- "What units are used throughout the implementation?"
- "What are the key mathematical assumptions made?"

### Algorithm Details
- "What are the primary algorithms implemented in your package?"
- "What numerical methods do you employ?" (e.g., integration, optimization)
- "Are there specific approximations or simplifications made in your implementation?"
- "What is the computational complexity of your key algorithms?"

### Physical Models
- "What physical principles does your package model?"
- "What are the key assumptions and limitations of your physical model?"
- "Are there specific parameter ranges where your model is valid?"
- "What physical constants are used in your implementation?"

### Data Analysis
- "What statistical or data analysis methods does your package employ?"
- "How does your package handle uncertainty or error propagation?"
- "Do you implement any Bayesian methods or inference techniques?"
- "What signal processing techniques are used, if any?"

### Visualization
- "What are the standard ways to visualize the outputs from your package?"
- "Do you use specific color maps or scaling approaches?"
- "Are there particular features users should look for in visualizations?"

### Validation
- "What test cases or validation approaches confirm your implementation is correct?"
- "Are there analytical solutions available for comparison?"
- "What benchmarks demonstrate the performance or accuracy of your implementation?"

### References
- "What are the key papers that your implementation is based on?"
- "Are there standard textbooks that provide background for your methods?"
- "What related software packages implement similar methods?"

## Synthesis and Refinement

1. After gathering information for each section, summarize your understanding
2. Present a draft of each section for confirmation or refinement
3. For sparse areas, suggest standard mathematical documentation from similar domains
4. Remove any sections the user indicates aren't relevant
5. Ask if they'd like to make any final adjustments before finalizing

## Output Generation

1. Format all equations using LaTeX notation (with proper delimiters: $$ for display equations, $ for inline)
2. Ensure all variables are properly defined after first use
3. Present the final document and create it in the project's .claude/commands directory:
   - Use the Write tool to create or update `.claude/commands/math_reference.md`
   - If the `.claude/commands` directory doesn't exist, suggest creating it with `mkdir -p .claude/commands`
   - Remind users to keep this file separate from math_reference_template.md
4. Explain that they can use this with `/claude math_reference.md` alongside their prime_claude.md to provide Claude with both code structure and mathematical context
5. Suggest integrating relevant portions into formal documentation for end users

## If User is Unsure

If the user is uncertain about certain mathematical aspects:
- Suggest standard approaches from the relevant scientific domain
- Provide examples of how similar packages document their mathematical foundations
- Offer simplified versions that can be expanded later
- For complex equations, suggest breaking them down into constituent parts

## Multi-Session Support

1. Acknowledge that creating detailed mathematical documentation might take multiple sessions
2. At the end of each session, summarize progress and what remains to be done
3. Suggest saving the current state of math_reference.md
4. When starting a new session, read the existing file and suggest where to continue

The goal is to create a math_reference.md file that effectively documents the mathematical and algorithmic foundations of the package, serving as both a reference for developers and a resource for understanding the implementation details.