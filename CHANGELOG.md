# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Modomics Integration**: Pre-computed modification annotations from MODOMICS database
  - 435 modification positions mapped across E. coli, S. cerevisiae, and H. sapiens
  - 100% mapping success to global coordinate system
  - Includes modification names, positions, Sprinzl coordinates, and structural regions
  - Output files: `outputs/modomics/*.tsv`
  - Documentation: `docs/development/MODOMICS_INTEGRATION.md`
  - Enables validation of detected modifications against known reference data

## [1.0.0] - 2025-01-15

### Added
- Initial release of tRNAs in space
- Core coordinate transformation script (`trnas_in_space.py`)
- Pre-computed global coordinates for three organisms:
  - *Escherichia coli* K12 MG1655 (87 tRNAs)
  - *Saccharomyces cerevisiae* S288C (292 tRNAs, nuclear + mitochondrial)
  - *Homo sapiens* GRCh38 (454 tRNAs, nuclear + mitochondrial)
- Comprehensive documentation:
  - README with quick start and usage examples
  - docs/OUTPUT_FORMAT.md with detailed column descriptions
  - docs/FAQ.md with practical Q&A
  - CONTRIBUTING.md for community guidelines
- Example Jupyter notebook for visualization and analysis
- Python package structure (pip installable)
- Automated testing and CI/CD:
  - Test suite with pytest
  - GitHub Actions workflows (tests, build, validation)
  - Multi-version Python support (3.9-3.12)
- Dataset metadata and provenance (METADATA.json)
- Standard project files:
  - MIT LICENSE
  - CITATION.cff for academic citations
  - pyproject.toml for modern Python packaging

### Infrastructure
- GitHub Actions workflows for continuous integration
- Automated testing across Python 3.9, 3.10, 3.11, 3.12
- Package build validation
- Repository health checks

### Documentation
- Complete API documentation in docs/OUTPUT_FORMAT.md
- Interactive examples in Jupyter notebook
- FAQ covering common use cases in docs/FAQ.md
- Citation information in standard formats

[Unreleased]: https://github.com/lkwhite/tRNAs-in-space/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/lkwhite/tRNAs-in-space/releases/tag/v1.0.0
