# ARCHIVED: Multi-File Coordinate System Implementation

---
## Completion Status

**Completed**: 2025-11-25
**Branch**: `feature/offset-type-coordinate-groups`

### Items Completed

1. **Extensible GroupingStrategy architecture** - Added `GroupKey`, `GroupingStrategy` base class, and `OffsetTypeStrategy` to `scripts/trnas_in_space.py`

2. **Fixed `build_pref_label()`** - Now prefers `sprinzl_label` (canonical position) over `sprinzl_index` for correct position alignment

3. **Added `--split-by-offset-and-type` CLI argument** - New mode for generating offset+type grouped coordinate files

4. **Fixed OUTPUT_FORMAT.md** - Corrected swapped descriptions of templateResidueIndex vs templateNumberingLabel

5. **Generated coordinate files for all organisms**:
   - E. coli: 5 files (offset -1/0/+1 × type1/2)
   - Yeast: 4 files (offset -1/0/+1 × type1/2)
   - Human: 7 files (offset -3/-1/0/+1 × type1/2)

6. **Validation**:
   - All 16 tests pass (10 existing + 6 new)
   - No collisions in any generated file
   - Position 55 aligns correctly within each offset+type group
   - Key modification sites (8, 34, 54, 55) all align correctly

### Deferred Items

- Clover integration (to be handled on clover side)
- Visualization/plotting of alignment comparison

---

# Original Handoff Document

**Date**: 2025-11-25
**Current Branch**: `feature/offset-type-coordinate-groups` in `/Users/laurawhite/Projects/tRNAs-in-space`
**Related Branch**: `feature/global-tRNA-coordinates` in `/Users/laurawhite/Projects/clover-globalcoords` (pushed)

## Context Summary

User reported that Sprinzl position 55 (and 8/11 other modification sites) don't align between Type I and Type II tRNAs in heatmaps. After extensive investigation, discovered this is a two-dimensional problem requiring a multi-file coordinate system.

## Key Findings

### 1. Documentation Error in OUTPUT_FORMAT.md
**Lines 30-39 have the field descriptions backwards:**
- Documentation says `templateResidueIndex` is canonical → WRONG
- Documentation says `templateNumberingLabel` is template-specific → WRONG

**Reality (proven by testing):**
- `templateNumberingLabel` = canonical Sprinzl position (e.g., position '18' is always '18')
- `templateResidueIndex` = template-specific, shifts with insertions (e.g., '18' might be index 18 or 19)

**Proof:** Pro with insertion '17a':
- Label '18' has index=19 (shifted)
- Ala without insertion: label '18' has index=18
- Both identify canonical position 18 via LABEL, not INDEX

### 2. Root Cause: Two-Dimensional Heterogeneity

**Dimension 1: Labeling Offset**
Different tRNA families have different offsets between label and index:
- E. coli: 3 offset groups (-1, 0, +1)
- Yeast: 2 offset groups (0, +1)

**Dimension 2: Structural Type**
- Type I: Standard tRNAs (~76 nt)
- Type II: Extended variable arm tRNAs (~85-90 nt, with e1-e24 positions)

### 3. Why Single-Dimension Solutions Failed

**Using `templateNumberingLabel` alone**: 70+ collisions
- Even though it's the canonical coordinate, still causes collisions
- Reason: Different tRNAs have incompatible coordinate schemes

**Offset normalization**: Still causes collisions
- Adjusting indices by offset doesn't resolve fundamental incompatibility

**Dual-system mode (Type I vs Type II)**: Creates separate coordinate spaces
- Position 55: Type I=68, Type II=77 (incomparable)
- Defeats purpose of unified coordinate system

**Offset grouping alone**: Still has collisions
- E.g., offset -1 group contains Pro (Type I), Ser (Type II), Tyr (Type II)
- Different structures → collisions between positions 47-76

## Solution: Multi-File System (offset × type)

User confirmed:
- Perfect alignment NOT essential at this phase
- Multiple coordinate files approach is acceptable
- Wants to see E. coli prototype/plot to verify it works

### Implementation Plan

Generate separate coordinate files for each **offset + type** combination:

**E. coli (6 files expected):**
```
ecoliK12_offset-1_type1.tsv  (~3 Type I tRNAs: Pro family)
ecoliK12_offset-1_type2.tsv  (~8 Type II tRNAs: Ser, Tyr)
ecoliK12_offset0_type1.tsv   (~55 Type I tRNAs: majority)
ecoliK12_offset0_type2.tsv   (~5 Type II tRNAs: some Leu)
ecoliK12_offset+1_type1.tsv  (~8 Type I tRNAs: various)
ecoliK12_offset+1_type2.tsv  (~3 Type II tRNAs: various)
```

**Yeast (4 files expected):**
```
sacCer_offset0_type1.tsv
sacCer_offset0_type2.tsv
sacCer_offset+1_type1.tsv
sacCer_offset+1_type2.tsv
```

## Key Files and Locations

### tRNAs-in-space (work here next)
- **Main script**: `scripts/trnas_in_space.py`
  - `build_pref_label()` at lines 324-346 (needs modification)
  - `classify_trna_type()` at line 88 (already exists)
- **R2DT data**: `r2dt_outputs/ecoli_jsons/`, `r2dt_outputs/yeast_jsons/`
- **Current outputs**: `outputs/ecoliK12_global_coords.tsv` (has misalignment)
- **This handoff**: `NEXT_SESSION_HANDOFF.md`

### clover-globalcoords (investigation docs)
- `SPRINZL_ALIGNMENT_ISSUE.md` - Initial problem report
- `ALIGNMENT_INVESTIGATION_SUMMARY.md` - Complete investigation
- Current branch: `feature/global-tRNA-coordinates` (pushed)

### Test scripts (in /tmp)
- `/tmp/test_ecoli_alignment_fix.py` - Tested simple templateNumberingLabel fix
- `/tmp/test_offset_grouping.py` - Tested offset-only grouping
- `/tmp/offset_normalization_prototype.py` - Tested offset normalization

## Implementation Details

### Offset Calculation (from test scripts)
```python
def calculate_trna_offset(trna_rows):
    """Calculate labeling offset using positions 15-25 (conserved region)"""
    offsets = []
    for row in trna_rows:
        label = row['sprinzl_label']
        idx = row['sprinzl_index']

        if label and label.isdigit() and idx > 0:
            label_num = int(label)
            if 15 <= label_num <= 25:
                offset = label_num - idx
                offsets.append(offset)

    if offsets:
        return max(set(offsets), key=offsets.count)  # Most common
    return None
```

### Modified build_pref_label() Approach
```python
def build_pref_label(df: pd.DataFrame) -> pd.Series:
    """Prefer templateNumberingLabel (canonical coordinate)"""
    lbl = df["sprinzl_label"].astype("string").fillna("").str.strip()
    num = pd.to_numeric(df["sprinzl_index"], errors="coerce")
    num_ok = num.where((num >= 1) & (num <= 76))
    num_ok_str = num_ok.astype("Int64").astype(str).replace({"<NA>": ""})

    # Prefer label (canonical), fallback to index when empty
    pref = lbl.mask(lbl.eq(""), other=num_ok_str)
    return pref
```

### Type Classification (already exists)
- `classify_trna_type()` at line 88
- Returns 'type1' or 'type2' based on presence of extended arm (e positions)

## Expected Outcomes

### Within Each Group
- ✅ No collisions (positions don't conflict)
- ✅ Positions align (same sprinzl_label → same global_index)
- ✅ Type I and Type II can be compared within the same offset group

### Across Groups
- ⚠️ Cannot directly compare (different coordinate spaces)
- ⚠️ Position 55 in offset0_type1 ≠ position 55 in offset+1_type1
- This is acceptable per user requirement (perfect alignment not essential)

## Validation Tests

For each generated file, verify:
1. No collisions: `validate_no_global_index_collisions(df)` passes
2. Alignment within group: Position 55 (and others) have same global_index across Type I/II
3. Common positions present: Test sites (8, 13, 16, 20, 26, 34, 37, 47, 54, 55, 58)
