# R2DT Label/Index Mismatch Bug Documentation

## Summary

R2DT's enriched JSON output contains two related but distinct position identifiers:
- `sprinzl_index` (templateResidueIndex): The structural position in R2DT's alignment
- `sprinzl_label` (templateNumberingLabel): The canonical Sprinzl position label

**The Bug**: In ~64 tRNAs across E. coli, yeast, and human, these two values are inconsistent in a way that breaks cross-tRNA alignment. Specifically, when R2DT detects a **deletion** (gap in the sequence), it correctly skips the `sprinzl_index` but **fails to skip the corresponding `sprinzl_label`**.

## The Pattern

### Normal Case (no deletion)
```
seq_index=20, sprinzl_index=20, sprinzl_label="20"  → consistent
seq_index=21, sprinzl_index=21, sprinzl_label="21"  → consistent
seq_index=22, sprinzl_index=22, sprinzl_label="22"  → consistent
```

### Bug Case (deletion at position 21)
```
seq_index=20, sprinzl_index=20, sprinzl_label="20"  → OK
seq_index=21, sprinzl_index=22, sprinzl_label="21"  → BUG! index skipped 21, but label didn't
seq_index=22, sprinzl_index=23, sprinzl_label="22"  → BUG propagates
```

**What should happen** (deletion at position 21):
```
seq_index=20, sprinzl_index=20, sprinzl_label="20"  → OK
seq_index=21, sprinzl_index=22, sprinzl_label="22"  → label should also skip 21
seq_index=22, sprinzl_index=23, sprinzl_label="23"  → consistent
```

## Why This Matters

The tRNAs-in-space coordinate system uses `sprinzl_label` to align tRNAs across isodecoders. When one tRNA has label="21" at a position that structurally corresponds to position 22, it gets placed at the wrong global_index in the heatmap.

**Visual effect**: The affected tRNA appears shifted relative to other tRNAs in the same family.

## Detection Logic

The bug can be detected by looking for this pattern:
```python
# If sprinzl_index jumps (indicating deletion)
if curr_idx > prev_idx + 1:
    skipped_positions = range(prev_idx + 1, curr_idx)
    # And the label falls within the skipped range
    if int(curr_label) in skipped_positions:
        # This is the bug - label should have skipped too
```

## Affected tRNAs (64 total)

### E. coli (17 tRNAs)
- tRNA-Arg-CCU-1-1
- tRNA-Asn-GUU-1-1 through 1-4
- tRNA-His-GUG-1-1
- tRNA-Leu-CAA-1-1, CAG-1-1 through 2-1, GAG-1-1, UAA-1-1, UAG-1-1
- tRNA-Tyr-GUA-1-1, 2-1, 2-2

### Yeast (21 tRNAs)
- nuc-tRNA-Asp-GUC-2-1
- nuc-tRNA-Leu-* (multiple isodecoders)
- nuc-tRNA-Tyr-AUA-1-1

### Human (27 tRNAs)
- nuc-tRNA-Leu-AAG-*, CAA-*, UAA-*, UAG-*
- nuc-tRNA-Pro-AGG-*, CGG-*, UGG-*

## Current State

### Manual Fixes Applied
We added `LABEL_OVERRIDES` dict in `scripts/trnas_in_space.py` for:
- `nuc-tRNA-Arg-CCU-1-1` (yeast) - the original case we investigated
- `nuc-tRNA-Leu-CAA-6-1` (human)
- `nuc-tRNA-His-GUG-1-1` through `1-9` (human)

### Tests Added
In `test_trnas_in_space.py`:
- `test_label_index_consistency_within_trna()` - stub (disabled, was too strict)
- `test_no_label_index_mismatch_at_deletion_sites()` - detects the bug pattern
- `test_label_overrides_applied()` - verifies override mechanism works
- `test_arg_ccu_alignment_fixed()` - regression test for original fix

## Proposed Automated Fix

Instead of manual overrides for 64 tRNAs, detect and fix automatically:

```python
def fix_label_index_mismatch(rows):
    """
    Auto-correct R2DT's label/index mismatch at deletion sites.

    When sprinzl_index skips a position (deletion), the sprinzl_label
    should also skip. R2DT sometimes fails to do this.
    """
    rows = sorted(rows, key=lambda r: r["seq_index"])

    for i in range(1, len(rows)):
        prev_idx = rows[i-1]["sprinzl_index"]
        curr_idx = rows[i]["sprinzl_index"]
        curr_label = rows[i]["sprinzl_label"]

        # Check for deletion (gap in index)
        if prev_idx > 0 and curr_idx > 0 and curr_idx > prev_idx + 1:
            if curr_label.isdigit():
                label_num = int(curr_label)
                skipped = set(range(prev_idx + 1, curr_idx))

                # If label is in the skipped range, it's wrong
                if label_num in skipped:
                    # Correct: label should equal index (for offset=0 tRNAs)
                    # or maintain consistent offset
                    rows[i]["sprinzl_label"] = str(curr_idx)

    return rows
```

### Concerns with Automated Fix

1. **Offset handling**: Some tRNAs have consistent non-zero offset (label = index + N). The fix needs to preserve this.

2. **Cascading corrections**: If we fix position 22, do we need to fix 23, 24, etc.? Yes - all subsequent positions need adjustment.

3. **Non-numeric labels**: Labels like "20a", "e5" need special handling.

4. **False positives**: Are there cases where label < index is intentional? Need to verify.

## Files Involved

- `scripts/trnas_in_space.py`: Main pipeline, contains `LABEL_OVERRIDES` and `collect_rows_from_json()`
- `test_trnas_in_space.py`: Tests including new label/index consistency tests
- `outputs/*_global_coords_offset*.tsv`: Generated coordinate files

## Next Steps

1. Fully understand when label < index is legitimate vs. a bug
2. Implement automated detection and correction
3. Validate that correction doesn't break working tRNAs
4. Remove manual `LABEL_OVERRIDES` once automated fix works
5. Add comprehensive tests

## Key Insight

The core issue is: **R2DT's `sprinzl_index` correctly reflects structural deletions, but `sprinzl_label` doesn't always follow suit.** Our coordinate system relies on `sprinzl_label` for alignment, so we need labels to be consistent with the structural reality that `sprinzl_index` represents.

## Commands to Reproduce

```bash
# Find all affected tRNAs
python3 << 'EOF'
import pandas as pd
from pathlib import Path

for f in Path("outputs").glob("*_global_coords_offset*.tsv"):
    df = pd.read_csv(f, sep="\t")
    for trna_id, group in df.groupby("trna_id"):
        group = group.sort_values("seq_index")
        rows = group.to_dict("records")
        for i in range(1, len(rows)):
            prev_idx, curr_idx = rows[i-1]["sprinzl_index"], rows[i]["sprinzl_index"]
            curr_label = str(rows[i]["sprinzl_label"]).strip()
            if prev_idx > 0 and curr_idx > 0 and curr_idx > prev_idx + 1:
                if curr_label.isdigit() and int(curr_label) in range(prev_idx+1, curr_idx):
                    print(f"{f.name}: {trna_id} seq={rows[i]['seq_index']} idx={curr_idx} label={curr_label}")
EOF

# Run tests
python test_trnas_in_space.py

# Regenerate coordinates (after any fix)
python scripts/trnas_in_space.py outputs/sacCer_jsons outputs/sacCer_global_coords.tsv --split-by-offset-and-type
python scripts/trnas_in_space.py outputs/hg38_jsons outputs/hg38_global_coords.tsv --split-by-offset-and-type
python scripts/trnas_in_space.py outputs/ecoliK12_jsons outputs/ecoliK12_global_coords.tsv --split-by-offset-and-type
```
