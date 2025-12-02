# tRNA Coordinate System

This document explains the unified coordinate system used for cross-tRNA alignment in this project.

## Overview

The coordinate system assigns a unique integer (`global_index`) to each position across all tRNAs, enabling:
- Cross-tRNA alignment and comparison
- Heatmap visualization with consistent X-axis
- Modification site analysis across isoacceptors

## Key Concepts

### Sprinzl Numbering

tRNA positions are named using **Sprinzl numbering** (Sprinzl et al., 1998):
- Canonical positions: 1-76 (standard tRNA)
- Letter suffixes for insertions: 20A, 20B (insertions after position 20)
- E-positions for Type II variable arm: e1-e27 (extended variable arm in Leu, Ser, Tyr)

### tRNA Types

**Type I (Standard)**: ~95% of tRNAs
- Short variable arm (positions 44-48)
- ~76 nucleotides total
- Examples: Ala, Arg, Asp, Glu, Gly, His, Ile, Lys, Met, Phe, Pro, Thr, Trp, Val

**Type II (Extended)**: Leu, Ser, Tyr
- Long variable arm with e-positions
- ~85-90 nucleotides total
- E-positions form a hairpin structure

### E-Position Biological Ordering

E-positions in Type II tRNAs form a hairpin, so they must be sorted in **biological order** (5'→3' along RNA backbone), not numerical order:

```
Numerical:  e1, e2, e3, ..., e27
Biological: e11, e12, e13, e14, e15, e16, e17, e1, e2, e3, e4, e5, e27, e26, e25, e24, e23, e22, e21
            ←───── stem 5' ────────────────→  ←─ loop ─→  ←───── stem 3' ────────────────────────→
```

This ordering ensures positions align correctly when comparing across tRNAs.

## How global_index is Built

### Step 1: Extract Sprinzl Labels

From R2DT JSON files, we extract `templateNumberingLabel` as the canonical Sprinzl position:
- "18" = canonical position 18
- "20A" = insertion after position 20
- "e12" = extended variable arm position 12

### Step 2: Sort Labels Globally

All unique labels across all tRNAs are sorted using biological ordering:
1. Numeric positions: 1 < 2 < ... < 20 < 20A < 20B < 21 < ...
2. E-positions at position 45 in biological hairpin order
3. Empty labels (insertions) sorted to end

### Step 3: Assign Ordinal Values

Each unique label gets an ordinal (integer rank):
```
Label    Ordinal
1        1
2        2
...
20       20
20A      21
20B      22
21       23
...
```

### Step 4: Generate Per-tRNA Continuous Coordinates

For each tRNA:
- Labeled positions get their ordinal value
- Unlabeled positions (insertions without canonical labels) are interpolated between neighbors

Example:
```
Position:   1    2    _    _    5
Ordinal:    1    2    -    -    5
Continuous: 1.0  2.0  3.0  4.0  5.0
```

### Step 5: Map to Integer global_index

All unique continuous values across all tRNAs are collected, sorted, and assigned integers 1..K:
```
Continuous   global_index
1.0          1
2.0          2
3.0          3
4.0          4
5.0          5
...
```

## Excluded tRNA Types

Some tRNAs are excluded from the unified coordinate system:
- **Mitochondrial tRNAs**: Non-standard structures
- **Selenocysteine (SeC)**: Unusual structure
- **Initiator Met (iMet)**: Special initiator tRNA with structural differences

## Validation

The coordinate system is validated to ensure no **collisions** occur:
- A collision is when two different Sprinzl labels map to the same global_index
- This would indicate a bug in the ordering or interpolation logic
- The script exits with an error if collisions are detected

## Example Usage

```python
import pandas as pd

# Load coordinates
df = pd.read_csv("hg38_global_coords.tsv", sep="\t")

# Create alignment matrix (tRNAs x positions)
alignment = df.pivot_table(
    index='trna_id',
    columns='global_index',
    values='residue',
    aggfunc='first'
)

# Find position 34 (anticodon wobble) across all tRNAs
wobble = df[df['sprinzl_label'] == '34']
```

## Technical Details

### Precision

Continuous coordinates are rounded to 6 decimal places before mapping to global_index. This ensures consistent grouping of identical positions across tRNAs.

### Sort Key Consistency

All sort keys return tuples with string third elements to ensure Python 3 comparison works correctly:
```python
# Good: consistent types
(20, 0, "")      # position 20
(20, 1, "A")     # position 20A
(45, 2, "005")   # e-position with bio_order 5 (zero-padded)

# Bad: mixed types (causes TypeError)
(45, 2, 5)       # integer
(45, 2, "A")     # string
```

## History

The coordinate system evolved through several iterations:

1. **Initial**: Single unified file (had collisions)
2. **Offset-type grouping**: Separate files per offset and type (worked around bugs)
3. **Current**: Single unified file with fixed bugs (no collisions)

The offset-type grouping was a workaround for two bugs:
1. Empty labels getting wrong index values
2. Sort key type inconsistencies

With these bugs fixed, a single unified coordinate file works for all organisms.
