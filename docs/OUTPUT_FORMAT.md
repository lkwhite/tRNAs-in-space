# Output Format Specification

This document describes the structure and contents of the TSV files produced by `trnas_in_space.py`.

## File Format

The output is a **tab-separated values (TSV)** file with a header row followed by one row per nucleotide position across all tRNAs in the dataset.

## Column Descriptions

| Column | Data Type | Range/Format | Description | Example Values |
|--------|-----------|--------------|-------------|----------------|
| `trna_id` | string | - | Unique identifier for each tRNA, derived from the input filename | `tRNA-Ala-GGC-1-1` |
| `source_file` | string | - | Original R2DT enriched JSON filename from which this data was extracted | `tRNA-Ala-GGC-1-1-B_Ala.enriched.json` |
| `seq_index` | integer | 1..N | Sequential position within the tRNA sequence (1-based, 5' to 3') | `1`, `42`, `76` |
| `sprinzl_index` | integer | -1 or 1..76 | Sprinzl position number, or -1 if unresolvable from JSON | `1`, `20`, `34`, `-1` |
| `sprinzl_label` | string | - | Sprinzl label with any suffixes (e.g., A, B, or e-notation) | `1`, `20A`, `17a`, `47:e1` |
| `residue` | string | A/C/G/U/T/etc. | Single-letter nucleotide code | `A`, `G`, `C`, `U`, `T` |
| `sprinzl_ordinal` | float | 1..K | Ordinal position in the global label ordering (sorted Sprinzl labels) | `1.0`, `42.5`, `89.0` |
| `sprinzl_continuous` | float | 0.0..K | Continuous coordinate for this specific tRNA (interpolated between labeled positions) | `1.0`, `20.333333`, `34.0` |
| `global_index` | integer | 1..K | Equal-spaced global position shared across all tRNAs (the unified coordinate) | `1`, `45`, `147` |
| `region` | string | predefined set | Structural region annotation based on Sprinzl position | `acceptor-stem`, `anticodon-loop` |

## Column Derivation

Understanding where each column comes from helps interpret the data correctly:

### Direct R2DT Outputs

**`sprinzl_label`** (from `templateNumberingLabel`)
- Canonical Sprinzl position (1-76 numbering, Sprinzl et al. 1998)
- Functionally consistent across tRNAs: position 34 is always the first anticodon base
- R2DT assigns these based on structural covariance models
- Example: Label '18' is always canonical position 18, regardless of insertions

**`sprinzl_index`** (from `templateResidueIndex`)
- Template-specific sequential position from R2DT's covariance model
- Adjusts for insertions: a tRNA with insertion "17a" will have shifted indices downstream
- Example: In Pro tRNA with insertion at 17a, label '18' maps to index 19 (shifted by 1)

**`residue`** (from `residueName`)
- Nucleotide base at each position

### Gap-Filled Columns

**`sprinzl_index`** (enhanced)
- When R2DT doesn't provide a clear index, the script infers it using forward/backward passes
- Ensures continuity where structurally reasonable
- Unresolvable positions remain as `-1`

### Derived Columns

**`sprinzl_ordinal`**
- Global ordering of unique positions across all tRNAs (1, 2, ..., 20, 20a, 21, ...)
- Built from canonical `sprinzl_index` values plus insertion `sprinzl_label` values
- Ensures consistent position ordering across the dataset

**`sprinzl_continuous`**
- Per-tRNA fractional coordinates for smooth interpolation
- Labeled positions receive integer values from `sprinzl_ordinal`
- Unlabeled positions interpolated fractionally between neighbors

**`global_index`**
- **Primary coordinate for cross-tRNA alignment**
- Maps unique `sprinzl_continuous` values to integers (1..K)
- Ensures functionally equivalent positions align across all tRNAs

**`region`**
- Structural region annotation derived from canonical `sprinzl_index`
- Based on Type I tRNA structure (acceptor-stem, D-loop, anticodon-loop, etc.)

### Why Two Position Systems?

R2DT provides both `sprinzl_label` (canonical/functional) and `sprinzl_index` (template-specific) because they serve different purposes:

- **`sprinzl_label`**: Use for functional comparisons (e.g., "find all first anticodon bases at position 34")
- **`sprinzl_index`**: Template-specific index that may shift with insertions; useful for template alignment

The derived `global_index` is built primarily from canonical `sprinzl_label` values to ensure functional equivalence across tRNAs, enabling proper alignment in heatmaps and cross-tRNA analyses.

## Detailed Column Explanations

### `trna_id`
- The unique identifier for each tRNA sequence
- Parsed from the R2DT JSON filename by removing suffixes like `-B_Ala`
- Used to group all nucleotides belonging to the same tRNA

**Example:**
```
tRNA-Ala-GGC-1-1
tRNA-Leu-CAA-2-3
tRNA-Met-CAT-1-1
```

### `source_file`
- Full filename of the R2DT enriched JSON file
- Useful for tracing back to original R2DT output
- Multiple tRNAs may share similar filenames with different suffixes

### `seq_index`
- Simple consecutive numbering from 5' to 3' end
- Always starts at 1 for each tRNA
- Does not account for structural features
- Maximum value equals the length of the tRNA sequence

### `sprinzl_index`
- Integer Sprinzl position from R2DT annotation
- Standard tRNA positions range from 1 to 76 (canonical Type I)
- Value of `-1` indicates the position could not be determined from the JSON
- Filled by inference using forward/backward passes along the sequence

**Inference logic:**
- If R2DT provides a clear Sprinzl position, it's used directly
- If missing, the script attempts to infer from neighboring positions
- Unresolvable cases remain as `-1`

### `sprinzl_label`
- Text representation of Sprinzl position, preserving all suffixes
- More detailed than `sprinzl_index`
- May include:
  - Simple numbers: `1`, `34`, `76`
  - Letter suffixes: `20A`, `20B`, `17a`
  - Insertion notation: `47:e1`, `47:e2`
  - Dotted notation: `9.1`, `9.2`

**Examples:**
```
1       (position 1)
20A     (insertion after position 20)
47:e1   (first variable loop insertion)
```

### `residue`
- Single nucleotide base at this position
- May include modified bases if present in the FASTA
- Common values: `A`, `C`, `G`, `U`, `T` (and modified bases if annotated)

### `sprinzl_ordinal`
- Global ordering of unique Sprinzl labels across all tRNAs
- Computed by sorting all unique labels: `1 < 2 < 20 < 20A < 20B < 21 ...`
- Allows comparison of relative positions across different tRNAs
- Float type to accommodate interpolated positions

**Sorting logic:**
- Numeric labels sorted numerically: `1 < 2 < 10 < 20`
- Suffixed labels come after base: `20 < 20A < 20B`
- Unknown/empty labels sorted to end

### `sprinzl_continuous`
- Continuous coordinate specific to each tRNA
- Labeled positions get integer values from `sprinzl_ordinal`
- Unlabeled positions are interpolated fractionally between adjacent labels
- Rounded to 6 decimal places for consistency

**Interpolation:**
```
Position:   1    2    _    _    5
Continuous: 1.0  2.0  3.0  4.0  5.0

Position:   20   _    _    21
Continuous: 20.0 20.333 20.667 21.0
```

### `global_index`
- **The main output coordinate system**
- Integer index from 1 to K (where K = number of unique continuous positions)
- Every unique `sprinzl_continuous` value maps to one `global_index`
- Shared across ALL tRNAs in the dataset
- Equal-spaced: suitable for heatmaps, alignments, plotting

**Key property:**
- Missing positions show as `NA` in downstream analysis
- Enables cross-tRNA comparisons on the same X-axis

### `region`
- Structural annotation based on Sprinzl position
- Predefined categories based on canonical tRNA structure
- Assigned using the leading number from `sprinzl_label` (falling back to `sprinzl_index`)

**Valid region values:**

| Region | Sprinzl Positions | Description |
|--------|-------------------|-------------|
| `acceptor-stem` | 1-7, 66-72 | Base-paired acceptor stem |
| `acceptor-tail` | 73-76 | 3' CCA tail and terminal nucleotides |
| `D-stem` | 10-13, 22-25 | Base-paired D arm stem |
| `D-loop` | 14-21 | D arm loop |
| `anticodon-stem` | 27-31, 39-43 | Base-paired anticodon arm stem |
| `anticodon-loop` | 32-38 | Anticodon loop (includes anticodon at 34-36) |
| `variable-region` | 44-46 | Variable region (Type I, short) |
| `variable-arm` | 47-48 | Extended variable arm (Type II, long) |
| `T-stem` | 49-53, 61-65 | Base-paired T arm stem |
| `T-loop` | 54-60 | T arm loop (contains TψC motif) |
| `unknown` | - | Position cannot be assigned to a region |

**Notes:**
- Based on Type I canonical tRNA structure
- Variable loop positions 44-48 can vary significantly in length
- Mitochondrial tRNAs may have non-standard structures

## Example Rows

```tsv
trna_id	source_file	seq_index	sprinzl_index	sprinzl_label	residue	sprinzl_ordinal	sprinzl_continuous	global_index	region
tRNA-Ala-GGC-1-1	tRNA-Ala-GGC-1-1-B_Ala.enriched.json	1	1	1	G	1.0	1.0	1	acceptor-stem
tRNA-Ala-GGC-1-1	tRNA-Ala-GGC-1-1-B_Ala.enriched.json	2	2	2	G	2.0	2.0	2	acceptor-stem
tRNA-Ala-GGC-1-1	tRNA-Ala-GGC-1-1-B_Ala.enriched.json	34	34	34	G	71.0	71.0	123	anticodon-loop
tRNA-Ala-GGC-1-1	tRNA-Ala-GGC-1-1-B_Ala.enriched.json	35	35	35	G	72.0	72.0	124	anticodon-loop
tRNA-Ala-GGC-1-1	tRNA-Ala-GGC-1-1-B_Ala.enriched.json	36	36	36	C	73.0	73.0	125	anticodon-loop
```

## Usage in Analysis

### Loading the Data

**Python:**
```python
import pandas as pd

df = pd.read_csv("ecoliK12_global_coords.tsv", sep="\t")
print(f"Loaded {len(df)} positions from {df['trna_id'].nunique()} tRNAs")
```

**R:**
```r
library(readr)
df <- read_tsv("ecoliK12_global_coords.tsv")
cat(sprintf("Loaded %d positions from %d tRNAs\n",
    nrow(df), n_distinct(df$trna_id)))
```

### Creating an Alignment Matrix

**Python (pandas):**
```python
# Pivot to create alignment matrix
alignment = df.pivot_table(
    index='trna_id',
    columns='global_index',
    values='residue',
    aggfunc='first'
)

# Missing positions will be NaN
print(alignment.shape)  # (n_trnas, max_global_index)
```

### Filtering by Region

**Python:**
```python
# Analyze only anticodon loop
anticodon = df[df['region'] == 'anticodon-loop']
print(anticodon.groupby('trna_id')['residue'].apply(''.join))
```

### Finding Specific Positions

**Find all position 34 (first anticodon position):**
```python
pos34 = df[df['sprinzl_index'] == 34]
print(pos34[['trna_id', 'residue']].to_string(index=False))
```

## File Size Expectations

Approximate file sizes for reference:

| Organism | tRNAs | Nucleotides | File Size | Global Positions |
|----------|-------|-------------|-----------|------------------|
| E. coli K12 | 87 | ~6,800 | ~600 KB | ~150 |
| S. cerevisiae | 297 | ~23,000 | ~2.1 MB | ~180 |
| H. sapiens | 454 | ~35,000 | ~3.3 MB | ~190 |

## Edge Cases and Special Situations

### Mitochondrial tRNAs
- May have non-standard structures
- Variable loop can be very short or absent
- D-loop may be reduced or missing
- Region annotations may be less accurate

### Type II tRNAs
- Long variable arms (positions 47-48 extended)
- Extra insertions in variable region
- More positions in `variable-arm` region

### Missing Positions
- `sprinzl_index` of `-1` indicates unresolvable position
- `global_index` may still be assigned based on sequence position
- `region` will be `unknown`

### Insertions and Deletions
- Insertions: shown with letter suffixes (20A, 20B) or e-notation (47:e1)
- Deletions: missing Sprinzl positions, gaps in `global_index`
- Each unique insertion gets its own `global_index`

## Quality Control

### Validation Checks

**Check for complete data:**
```python
# Ensure no missing critical columns
assert df['trna_id'].notna().all()
assert df['residue'].notna().all()
assert df['seq_index'].notna().all()

# Check global_index continuity
global_positions = sorted(df['global_index'].dropna().unique())
assert global_positions[0] == 1
assert len(global_positions) == global_positions[-1]  # Should be continuous
```

**Check region assignments:**
```python
valid_regions = [
    'acceptor-stem', 'acceptor-tail', 'D-stem', 'D-loop',
    'anticodon-stem', 'anticodon-loop', 'variable-region',
    'variable-arm', 'T-stem', 'T-loop', 'unknown'
]
assert df['region'].isin(valid_regions).all()
```

## Version History

- **v1.0.0** (2025-01-01): Initial format specification
  - 10 columns
  - Region annotations
  - 6 decimal precision for continuous coordinates

## Related Documentation

- [README.md](../README.md) - Overall project documentation
- [FAQ.md](FAQ.md) - Frequently asked questions
- [test_trnas_in_space.py](../test_trnas_in_space.py) - Validation tests
