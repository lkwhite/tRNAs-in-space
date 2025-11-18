# Pre-computed tRNA Global Coordinates

This directory contains pre-computed global coordinate mappings for tRNAs from common model organisms.

## Available Datasets

**Production-ready coordinate files** (validated, collision-free):

| File | Organism | tRNAs | Nucleotides | Global Positions | Size |
|------|----------|-------|-------------|------------------|------|
| `ecoliK12_global_coords_fixed.tsv` | *E. coli* K12 MG1655 | 82 | 7,075 | 130 | 552 KB |
| `sacCer_global_coords_fixed.tsv` | *S. cerevisiae* S288C | 268 | 23,119 | 150 | 1.9 MB |
| `hg38_global_coords_fixed.tsv` | *H. sapiens* GRCh38 | 422 | 36,356 | 265 | 3.0 MB |

**Note**: All files contain only **nuclear elongator tRNAs**. Selenocysteine, mitochondrial, and initiator tRNAs are excluded by design (see [ANALYSIS_GUIDELINES.md](../ANALYSIS_GUIDELINES.md)).

## Quick Start

### Python

```python
import pandas as pd

# Load E. coli data
df = pd.read_csv('outputs/ecoliK12_global_coords_fixed.tsv', sep='\t')

# Create alignment matrix
alignment = df.pivot_table(
    index='trna_id',
    columns='global_index',
    values='residue'
)

# Filter by region
anticodon_loop = df[df['region'] == 'anticodon-loop']
```

### R

```r
library(readr)
library(dplyr)

# Load yeast data
df <- read_tsv('outputs/sacCer_global_coords_fixed.tsv')

# Filter nuclear tRNAs (exclude mito)
nuclear <- df %>%
  filter(!grepl('mito', trna_id, ignore.case = TRUE))

# Extract anticodon sequences
anticodons <- df %>%
  filter(sprinzl_index %in% 34:36) %>%
  group_by(trna_id) %>%
  summarize(anticodon = paste(residue, collapse = ''))
```

## File Format

Each TSV file contains 10 columns per nucleotide position. These columns come from three sources:

### Direct R2DT Outputs

- **`sprinzl_index`**: Canonical Sprinzl position (1-76), functionally consistent across tRNAs
- **`sprinzl_label`**: Template alignment position, adjusts for insertions (e.g., "20a")
- **`residue`**: Nucleotide base

R2DT provides two position systems because:
- `sprinzl_index` represents the **functional position** (position 34 = first anticodon base, always)
- `sprinzl_label` represents the **alignment position** in R2DT's template (may vary with insertions)

### Gap-Filled Columns

- **`sprinzl_index`**: Enhanced with inference for positions R2DT couldn't annotate

### Derived Coordinate Transformations

- **`sprinzl_ordinal`**: Global ordering of positions (1, 2, ..., 20, 20a, 21, ...)
- **`sprinzl_continuous`**: Per-tRNA fractional coordinates with interpolation
- **`global_index`**: Equal-spaced integer coordinates for cross-tRNA alignment (**use this for plotting**)
- **`region`**: Structural region annotation (acceptor-stem, D-loop, anticodon-loop, etc.)

### Metadata Columns

- **`trna_id`**: Unique tRNA identifier
- **`source_file`**: Original R2DT JSON filename
- **`seq_index`**: Sequential position (1-based, 5' to 3')

**For detailed explanations**, see [OUTPUT_FORMAT.md](../docs/OUTPUT_FORMAT.md).

## Metadata

Detailed metadata for all datasets is available in `METADATA.json`, including:
- Generation dates and versions
- Source databases
- Processing pipeline details
- Quality notes and warnings
- Validation information

## Important Notes

### Scope and Filtering

**These datasets contain only nuclear elongator tRNAs** for structural consistency:

**✅ Included:**
- Type I tRNAs (standard structure, ~76 nt)
- Type II tRNAs (extended variable arm: Leu, Ser, Tyr, ~90 nt)

**🚫 Excluded by design:**
- **Selenocysteine tRNAs** - Structurally incompatible (~95 nt, unique binding requirements)
- **Mitochondrial tRNAs** - Different architecture (60-75 nt, missing canonical features)
- **Initiator methionine tRNAs** - Modified structure for specialized ribosome binding

For analysis of excluded tRNA types, see [ANALYSIS_GUIDELINES.md](../ANALYSIS_GUIDELINES.md) for alternative approaches.

### Coordinate System Success

All files have been validated with **zero collisions** across all organisms:
- ✅ E. coli: 82 nuclear elongator tRNAs, 7,075 positions
- ✅ Yeast: 268 nuclear elongator tRNAs, 23,119 positions
- ✅ Human: 422 nuclear elongator tRNAs, 36,356 positions

See [docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SUCCESS_REPORT.md](../docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SUCCESS_REPORT.md) for technical details.

## Data Provenance

All sequences were obtained from [GtRNAdb](https://gtrnadb.org) and processed using:
1. **R2DT 2.0** - For Sprinzl position annotation and structural alignment
2. **trnas_in_space.py** - For coordinate transformation and global index generation
3. **Filtering** - Nuclear elongator tRNAs only (excludes SeC, mitochondrial, initiator)

## Regenerating Data

To regenerate these files or create coordinates for new organisms:

```bash
# Step 1: Run R2DT
docker run --rm \
  -v "$(pwd):/data" \
  rnacentral/r2dt \
  r2dt.py gtrnadb draw /data/fastas/your_trnas.fa /data/outputs/your_jsons

# Step 2: Generate global coordinates
python scripts/trnas_in_space.py outputs/your_jsons outputs/your_coords.tsv
```

## Quality Assurance

All pre-computed files have been validated to ensure:
- Correct column structure
- Continuous global indices
- Valid region annotations
- No missing critical data

Run validation tests:
```bash
pytest test_trnas_in_space.py::test_output_file_structure
pytest test_trnas_in_space.py::test_global_index_continuity
```

## Examples

See the [examples](../examples/) directory for:
- `01_basic_visualization.ipynb` - Interactive data exploration
- Heatmap generation
- Cross-isodecoder comparisons
- Region-specific analysis

## Citation

If you use these pre-computed coordinates, please cite:

```bibtex
@software{trnas_in_space,
  author = {White, Laura K.},
  title = {tRNAs in space: Standardized coordinates for tRNA analysis},
  year = {2025},
  url = {https://github.com/lkwhite/tRNAs-in-space},
  license = {MIT}
}
```

## License

These datasets are released under the MIT License. See [LICENSE](../LICENSE) for details.

## Contact

Questions or issues? Please [open an issue](https://github.com/lkwhite/tRNAs-in-space/issues) on GitHub.
