# Pre-computed tRNA Global Coordinates

This directory contains pre-computed global coordinate mappings for tRNAs from common model organisms.

## Available Datasets

| File | Organism | tRNAs | Nucleotides | Global Positions | Size |
|------|----------|-------|-------------|------------------|------|
| `ecoliK12_global_coords.tsv` | *E. coli* K12 MG1655 | 87 | 6,793 | 130 | 588 KB |
| `sacCer_global_coords.tsv` | *S. cerevisiae* S288C | 292 | 21,877 | 150 | 2.1 MB |
| `hg38_global_coords.tsv` | *H. sapiens* GRCh38 | 454 | 33,819 | 265 | 3.2 MB |

## Quick Start

### Python

```python
import pandas as pd

# Load E. coli data
df = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')

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
df <- read_tsv('outputs/sacCer_global_coords.tsv')

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

Each TSV file contains the following columns:

- `trna_id`: Unique tRNA identifier
- `source_file`: Original R2DT JSON filename
- `seq_index`: Sequential position (1-based)
- `sprinzl_index`: Sprinzl position (1-76)
- `sprinzl_label`: Sprinzl label with suffixes
- `residue`: Nucleotide base
- `sprinzl_ordinal`: Ordinal in global label order
- `sprinzl_continuous`: Continuous coordinate
- `global_index`: Equal-spaced global position (1..K)
- `region`: Structural region annotation

See [OUTPUT_FORMAT.md](../OUTPUT_FORMAT.md) for detailed documentation.

## Metadata

Detailed metadata for all datasets is available in `METADATA.json`, including:
- Generation dates and versions
- Source databases
- Processing pipeline details
- Quality notes and warnings
- Validation information

## Important Notes

### Mitochondrial tRNAs

**S. cerevisiae mitochondrial tRNAs** may require additional manual curation due to:
- Lack of fungal mitochondrial-specific models in R2DT
- Non-standard tRNA structures in mitochondria
- Potential alignment inaccuracies

If working with yeast mitochondrial tRNAs, we recommend:
1. Manual review of alignments
2. Consulting Reinsch & Garcia 2025 for validated alignments
3. Using these coordinates as a starting point, not ground truth

## Data Provenance

All sequences were obtained from [GtRNAdb](https://gtrnadb.org) and processed using:
1. **R2DT 2.0** - For Sprinzl position annotation
2. **trnas_in_space.py** - For coordinate transformation

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
