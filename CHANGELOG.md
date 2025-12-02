# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Mitochondrial tRNA Coordinate Support**: New `--mito` flag for separate mito coordinate files
  - Human mito: 22 tRNAs, 111 unique positions, no collisions
  - Yeast mito: 19 tRNAs, 167 unique positions, no collisions
  - Separate files: `{species}_mito_global_coords.tsv`
  - T-loop and anticodon validation skipped for mito (different biology)
  - Documentation updated: `docs/OUTPUT_FORMAT.md`
- **Biological Validation Tests**: Strict quality checks for coordinate accuracy
  - Anticodon validation: positions 34-35-36 must match tRNA name
  - T-loop validation: positions 54-55-56 must be TTC/TTT or *TC variant
  - 100% pass rate required (no threshold-based testing)
  - Documentation: `docs/BIOLOGICAL_VALIDATION.md`
- **T-loop Variant Discovery**: Identified 7 human tRNAs with valid *TC T-loop variants
  - CTC: Ile-GAU (3 copies), Gly-UCC, Lys-CUU
  - ATC: Val-AAC
  - GTC: Val-UAC
- **Excluded tRNAs Documentation**: Comprehensive list of excluded tRNAs with reasons
  - Documentation: `docs/EXCLUDED_TRNAS.md`
- **Modomics Integration**: Pre-computed modification annotations from MODOMICS database
  - 435 modification positions mapped across E. coli, S. cerevisiae, and H. sapiens
  - 100% mapping success to global coordinate system
  - Includes modification names, positions, Sprinzl coordinates, and structural regions
  - Output files: `outputs/modomics/*.tsv`
  - Documentation: `docs/development/MODOMICS_INTEGRATION.md`
  - Enables validation of detected modifications against known reference data
- **Multi-species Processing Infrastructure**: Tools for batch processing multiple organisms
  - Organism configuration system (`config/organisms.yaml`)
  - Automated processing pipeline (`scripts/process_organisms.py`)
  - Download helper for GtRNAdb FASTA files (`scripts/download_gtrnadb_fastas.py`)
  - Comprehensive download guide (`docs/DOWNLOAD_GUIDE.md`)
  - Support for 14 Tier 1 model organisms (11 new: mouse, fly, worm, zebrafish, arabidopsis, etc.)

### Changed
- **tRNA Counts Updated**: Reflect exclusions for annotation quality
  - E. coli K12: 82 tRNAs (unchanged)
  - S. cerevisiae: 267 tRNAs (was 268, excluded 1)
  - H. sapiens: 416 tRNAs (was 422, excluded 6)
- **Documentation Organization**: Reorganized documentation structure
  - Created `docs/archive/` for historical planning documents
  - Moved completed planning docs to `docs/archive/planning/`
  - Moved completed PR summaries to `docs/archive/completed-prs/`
  - Added archive README explaining historical context
  - Keeps `docs/development/` focused on active reference documentation

### Fixed
- **Removed fix_label_index_mismatch() function**: Was incorrectly "fixing" R2DT labels
  - R2DT labels were actually correct; the pattern was misdiagnosed
  - The function was shifting anticodons by +1, breaking 21 yeast tRNAs
  - Anticodon validation tests now catch this class of bug

### Removed
- Archived `docs/R2DT_LABEL_INDEX_MISMATCH_BUG.md` - documented a misdiagnosed issue
- Archived `docs/SESSION_HANDOFF_ARG_CCU_BUG.md` - related to same misdiagnosis

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
