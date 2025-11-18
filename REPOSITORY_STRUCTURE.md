# Repository Structure

This document provides an overview of the tRNAs-in-space repository organization.

## Repository Layout

```
tRNAs-in-space/
├── README.md                          # Main project documentation
├── CHANGELOG.md                       # Version history and release notes
├── CONTRIBUTING.md                    # Contribution guidelines
├── ANALYSIS_GUIDELINES.md             # Research usage guidelines
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation information
│
├── scripts/                           # Production scripts
│   ├── trnas_in_space.py             # Main coordinate generation script
│   ├── process_organisms.py          # Multi-organism processing
│   ├── download_gtrnadb_fastas.py    # FASTA download utility
│   ├── fix_e_position_global_index.py # Position fixing utility
│   ├── regenerate_coords.sh          # Batch regeneration script
│   └── modomics/                     # Modomics integration scripts
│
├── outputs/                           # Production-ready outputs
│   ├── README.md                     # Output file documentation
│   ├── METADATA.json                 # Output metadata
│   ├── ecoliK12_global_coords_fixed.tsv    # E. coli coordinates
│   ├── sacCer_global_coords_fixed.tsv      # S. cerevisiae coordinates
│   ├── hg38_global_coords_fixed.tsv        # H. sapiens coordinates
│   └── modomics/                     # Modomics annotation files
│       └── modomics_to_sprinzl_mapping.tsv
│
├── docs/                              # Documentation
│   ├── OUTPUT_FORMAT.md              # Output specification
│   ├── FAQ.md                        # Frequently asked questions
│   ├── DOWNLOAD_GUIDE.md             # Data download instructions
│   └── archive/                      # Historical documentation
│       ├── README.md                 # Archive organization guide
│       ├── coordinate-fixes/         # Coordinate system development
│       │   ├── COORDINATE_SYSTEM_SUCCESS_REPORT.md
│       │   ├── COORDINATE_FIX_VALIDATION_REPORT.md
│       │   ├── TRNAS_IN_SPACE_COORDINATE_FIX_BRIEFING.md
│       │   ├── COORDINATE_SYSTEM_SCOPE.md
│       │   ├── COORDINATE_SYSTEM_VALIDATION.md
│       │   └── GLOBAL_COORDINATE_CAPABILITIES.md
│       ├── modomics-integration/     # Modomics development docs
│       │   ├── MODOMICS_INTEGRATION.md
│       │   ├── modomicscodes.csv
│       │   ├── modified_tRNA_all_all_rna_sequences.fasta
│       │   └── unmodified_tRNA_all_all_rna_sequences.fasta
│       ├── completed-prs/            # PR documentation
│       │   ├── PHASE2_SUMMARY.md
│       │   ├── PR_DESCRIPTION.md
│       │   ├── PR_SUMMARY.md
│       │   └── QUICK_WINS_SUMMARY.md
│       └── planning/                 # Historical planning docs
│           ├── NEXT_SESSION_START.md
│           ├── RECOMMENDATIONS.md
│           └── REMAINING_PLAN.md
│
├── examples/                          # Usage examples
│   └── 01_basic_visualization.ipynb  # Interactive tutorial
│
├── fastas/                            # Input FASTA files (user-provided)
│
├── config/                            # Configuration files
│
├── .github/                           # GitHub workflows and actions
│
├── .claude/                           # Claude Code configuration
│
├── test_trnas_in_space.py            # Test suite
├── pyproject.toml                    # Python package configuration
├── requirements.txt                  # Python dependencies
├── Makefile                          # Build automation
└── build.sh                          # Build script
```

## Key Directories

### `/outputs/`
Contains **production-ready coordinate files** for three model organisms:
- All files use the `*_fixed.tsv` naming convention
- Each file contains validated, collision-free global coordinates
- Ready for use in downstream analyses and visualizations

### `/scripts/`
Python scripts for generating and processing tRNA coordinates:
- `trnas_in_space.py` - Core coordinate generation from R2DT output
- `process_organisms.py` - Batch processing for multiple organisms
- `modomics/` - Modomics database integration utilities

### `/docs/`
User-facing documentation:
- Active documentation for current users
- `/archive/` contains historical and development documentation

### `/docs/archive/`
Historical documentation organized by topic:
- **`coordinate-fixes/`** - Development of the coordinate system (2024-2025)
- **`modomics-integration/`** - Modomics database integration work
- **`completed-prs/`** - Documentation from merged pull requests
- **`planning/`** - Historical planning and recommendations

### `/examples/`
Jupyter notebooks demonstrating usage:
- `01_basic_visualization.ipynb` - Interactive tutorial for basic usage

## File Naming Conventions

### Output Files
- **Production files**: `{organism}_global_coords_fixed.tsv`
  - `ecoliK12` - *Escherichia coli* K12
  - `sacCer` - *Saccharomyces cerevisiae*
  - `hg38` - *Homo sapiens* (GRCh38)

### Documentation Files
- **User documentation**: Uppercase with underscores (e.g., `OUTPUT_FORMAT.md`)
- **Archive documentation**: Preserved with original naming

## What's NOT in the Repository

The following are excluded via `.gitignore`:
- `__pycache__/` - Python bytecode cache
- `.pytest_cache/` - Pytest cache
- `*.pyc`, `*.pyo` - Compiled Python files
- `.DS_Store` - macOS system files
- Virtual environment directories

## Archive Organization

The `docs/archive/` directory preserves historical context while keeping the main repository clean:

1. **Coordinate system development** (`coordinate-fixes/`) - Documents the evolution and validation of the global coordinate system
2. **Modomics integration** (`modomics-integration/`) - Development notes for integrating the Modomics modification database
3. **Completed work** (`completed-prs/`, `planning/`) - Historical planning and PR documentation

All archived documents remain accessible but don't clutter the main repository navigation.

## Updates and Maintenance

This structure was established in November 2025 following project completion. Key changes:
- Consolidated 6 coordinate system docs → `docs/archive/coordinate-fixes/`
- Removed duplicate/old output files (kept only `*_fixed.tsv` versions)
- Archived development documentation → `docs/archive/modomics-integration/`
- Updated all documentation links to reflect new structure

For questions about repository organization, see [CONTRIBUTING.md](CONTRIBUTING.md) or open an issue.
