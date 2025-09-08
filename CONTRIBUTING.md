# Contributing to TenK10K Causal Inference Analysis

Thank you for your interest in contributing to this project! This document provides guidelines for contributing to the TenK10K causal inference analysis repository.

## Types of Contributions

We welcome several types of contributions:

### 🐛 Bug Reports
- Report issues with code, documentation, or analysis results
- Include reproducible examples and system information
- Check existing issues before creating new ones

### 📖 Documentation
- Improve README files and code documentation
- Add examples or clarify instructions
- Fix typos or formatting issues

### 🔧 Code Improvements
- Fix bugs or improve existing functionality
- Optimize performance or memory usage
- Add new features that enhance the analysis pipeline

### 📊 Analysis Enhancements
- Improve statistical methods or visualizations
- Add new validation approaches
- Extend analyses to new datasets or traits

## Getting Started

### 1. Fork and Clone
```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR-USERNAME/tenk10k-causal-paper.git
cd tenk10k-causal-paper
git remote add upstream https://github.com/powellgenomicslab/tenk10k-causal-paper.git
```

### 2. Set Up Development Environment
```bash
# Create conda environment
conda env create -f environment.yml
conda activate tenk10k-causal-paper

# Or install requirements manually
pip install -r requirements.txt
```

### 3. Create a Branch
```bash
git checkout -b feature/your-feature-name
# or
git checkout -b bugfix/issue-description
```

## Development Guidelines

### Code Style

#### R Code
- Follow [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable names
- Comment complex logic
- Keep functions focused and modular

```r
# Good
calculate_mr_effect <- function(beta_gwas, beta_eqtl, se_gwas, se_eqtl) {
  # Calculate MR effect using ratio method
  mr_beta <- beta_gwas / beta_eqtl
  mr_se <- abs(mr_beta) * sqrt((se_gwas / beta_gwas)^2 + (se_eqtl / beta_eqtl)^2)
  
  return(list(beta = mr_beta, se = mr_se))
}
```

#### Python Code
- Follow [PEP 8](https://pep8.org/) style guide
- Use type hints where appropriate
- Write docstrings for functions and classes

```python
def process_single_cell_data(adata: AnnData, 
                           min_genes: int = 200,
                           min_cells: int = 3) -> AnnData:
    """
    Process single-cell data with basic quality control.
    
    Parameters
    ----------
    adata : AnnData
        Raw single-cell data
    min_genes : int
        Minimum genes per cell
    min_cells : int
        Minimum cells per gene
        
    Returns
    -------
    AnnData
        Filtered single-cell data
    """
    # Implementation here
    return adata
```

### Documentation

- Update README files when adding new functionality
- Include examples in function documentation
- Explain the rationale behind analysis choices
- Document data requirements and expected outputs

### Testing

While formal unit tests are not required, please:
- Test your changes with representative data
- Verify that existing functionality still works
- Document any new dependencies or requirements

## Submitting Changes

### 1. Commit Guidelines
- Write clear, descriptive commit messages
- Use present tense ("Add feature" not "Added feature")
- Reference issues when applicable

```bash
git add .
git commit -m "Add drug target enrichment visualization

- Implement heatmap for target categories
- Add statistical significance testing
- Update documentation with examples

Fixes #123"
```

### 2. Push and Create Pull Request
```bash
git push origin feature/your-feature-name
```

Then create a pull request on GitHub with:
- Clear description of changes
- Reference to related issues
- Screenshots for visual changes
- List of breaking changes (if any)

### 3. Pull Request Review
- Address reviewer feedback promptly
- Update documentation as needed
- Ensure all discussions are resolved

## Specific Areas for Contribution

### High Priority
- **Documentation improvements** - README files, code comments, examples
- **Reproducibility enhancements** - Better dependency management, containerization
- **Visualization improvements** - More informative plots, interactive figures
- **Performance optimization** - Memory usage, computational efficiency

### Medium Priority
- **Additional validation methods** - New statistical approaches, external datasets
- **Extended analyses** - Additional traits, cell types, or populations
- **Code refactoring** - Modularization, reusability improvements

### Ideas for New Contributors
- Fix documentation typos or unclear instructions
- Add examples to existing functions
- Improve error messages and input validation
- Create tutorial notebooks for specific analyses

## Community Guidelines

### Be Respectful
- Use welcoming and inclusive language
- Be respectful of different viewpoints and experiences
- Focus on constructive feedback

### Be Collaborative
- Help others learn and contribute
- Share knowledge and resources
- Credit others for their contributions

### Be Patient
- Understand that reviews take time
- Be patient with new contributors
- Ask for clarification when needed

## Questions?

If you have questions about contributing:
- Check existing issues and discussions
- Open a new issue with the "question" label
- Reach out to maintainers directly

## Recognition

Contributors will be acknowledged in:
- Git commit history
- Future publications (for significant contributions)
- Project documentation

Thank you for contributing to advancing causal inference methods in single-cell genomics!