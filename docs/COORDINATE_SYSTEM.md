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

#### Auto-Fill Missing Labels

R2DT templates sometimes fail to assign Sprinzl labels even when the position is unambiguous. The `auto_fill_missing_labels()` function detects and corrects these gaps:

**Pattern detected**: `labeled(5) → unlabeled → labeled(7)`

When a single unlabeled nucleotide sits between two labeled positions with a numeric difference of exactly 2, the missing label is filled in (position 6 in this example).

**Before fix**: 45-55% of tRNAs had missing labels at positions like 6, 13, 22, 67
**After fix**: 0% affected

This only fills unambiguous cases (gap of exactly 2, purely numeric labels).

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
- Unlabeled positions (insertions without canonical labels) are interpolated using **fixed-slot alignment**

#### Fixed-Slot Alignment for Insertions

When tRNAs have different numbers of insertions between the same pair of canonical positions, we use fixed slots to ensure alignment:

1. **Compute max insertions per gap**: Find the maximum number of insertions between each pair of canonical positions across ALL tRNAs
2. **Allocate fixed slots**: Reserve that many fractional positions for the gap
3. **Assign left-aligned**: Each tRNA's insertions map to slots starting from the left

**Example**: If one tRNA has 2 insertions and another has 3 between positions 21→22:
```
Max insertions = 3, so allocate slots: 21.25, 21.50, 21.75

tRNA A (2 insertions): 21.0 → 21.25 → 21.50 → 22.0
tRNA B (3 insertions): 21.0 → 21.25 → 21.50 → 21.75 → 22.0
```

All tRNAs with insertions at the same site share the same `global_index` columns, ensuring clean heatmap alignment.

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

### Structurally Incompatible
- **Mitochondrial tRNAs**: Non-standard structures (60-75 nt vs 76)
- **Selenocysteine (SeC)**: Unusual structure (~95 nt, extended variable arm)
- **Initiator Met (iMet)**: Special initiator tRNA with structural differences

### Poorly Annotated (R2DT Issues)
- **Missing anticodon labels**: Can't verify tRNA identity
- **Wrong anticodon**: R2DT annotation doesn't match tRNA name
- **Invalid T-loop**: Non-*TC pattern indicates shifted annotations

See [EXCLUDED_TRNAS.md](EXCLUDED_TRNAS.md) for the complete list with reasons.

## Validation

The coordinate system is validated at two levels:

### Collision Detection
- A collision is when two different Sprinzl labels map to the same global_index
- This would indicate a bug in the ordering or interpolation logic
- The script exits with an error if collisions are detected

### Biological Validation
Every tRNA must pass two biological checks:
- **Anticodon match**: Positions 34-35-36 must match the anticodon in the tRNA name
- **T-loop validity**: Positions 54-55-56 must contain TTC, TTT, or a *TC variant

tRNAs failing these checks are excluded. See [BIOLOGICAL_VALIDATION.md](BIOLOGICAL_VALIDATION.md) for details.

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
3. **Unified (v1.1)**: Single unified file with fixed bugs (no collisions)
4. **Current (v1.2)**: Fixed-slot alignment for insertions + auto-fill missing labels

The offset-type grouping was a workaround for two bugs:
1. Empty labels getting wrong index values
2. Sort key type inconsistencies

With these bugs fixed, a single unified coordinate file works for all organisms.

**v1.2 improvements:**
- Fixed-slot alignment ensures tRNAs with different insertion counts share the same global_index columns
- Auto-fill missing labels corrects R2DT template gaps at positions like 6, 13, 22, 67
- Yeast coordinates reduced from 133 to 115 unique positions (cleaner alignment)
