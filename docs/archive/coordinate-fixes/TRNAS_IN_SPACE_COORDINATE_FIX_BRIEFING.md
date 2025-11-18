# tRNAs-in-Space Coordinate System Fix - Technical Briefing

## Problem Summary

The current coordinate system in tRNAs-in-space has fundamental flaws that make cross-tRNA structural alignment impossible for visualization. A recent "fix" for extended variable arms created massive global_index collisions that broke the entire coordinate system.

## Critical Issues Discovered

### 1. Global Index Collisions (Current State)
The recent extended variable arm fix created coordinate conflicts:

```
Type I (Ala):    position 47 = global_index 63
Type II (Leu):   e12        = global_index 63  ❌ COLLISION
SeC:             47a        = global_index 65  ❌ COLLISION with Ala-49
```

**Result**: Multiple structurally different positions share the same global_index, breaking cross-tRNA alignment.

### 2. Selenocysteine Incompatibility
SeC tRNAs are **structurally incompatible** with the global coordinate approach:
- **95 nucleotides** vs. standard 76
- **17+ nucleotide extended variable arm** (47a-47r)
- **All downstream positions shifted** by ~20 nucleotides
- **Cannot be meaningfully aligned** with standard tRNAs

### 3. Space Requirements Analysis
Three extended variable arm systems compete for the same coordinate space:

| **tRNA Type** | **Variable Arm** | **Space Needed** | **Current global_index** |
|---------------|-----------------|------------------|-------------------------|
| Type I (Standard) | 47-48 (2 nt) | 2 positions | 63-64 |
| Type II (Leu/Ser/Tyr) | e1-e24 (~15 nt) | 15 positions | 63-77 |
| SeC | 47a-47r (17 nt) | 17 positions | 65-85 |
| **TOTAL** | | **34 positions** | **All conflicting** |

## Proposed Solution: Type I + Type II Coordinate System

### Core Strategy
1. **Exclude SeC tRNAs entirely** - treat as separate structural class
2. **Design coordinate system** for Type I + Type II compatibility only
3. **Reserve coordinate space** (e.g., 65-85) for Type II extended arms
4. **Maintain structural alignment** for all other regions

### Implementation Plan

#### 1. Modify `sort_key()` Function
```python
def sort_key(s):
    """Enhanced sort key that properly handles Type I/II extended arms while excluding SeC"""

    # Standard numeric positions
    if s.isdigit():
        return (int(s), 0, 0)

    # Type II extended variable arms - map to reserved space
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        e_num = int(m.group(1))
        # Map to reserved coordinate space (65-85)
        return (64, 1, e_num)  # Places after position 64, before regular 65

    # Handle other extensions...
    # NOTE: SeC positions (47a-47r) should be handled in separate coordinate system
```

#### 2. Add SeC Filtering
```python
def should_exclude_trna(trna_id):
    """Filter out structurally incompatible tRNAs"""
    # Exclude selenocysteine - incompatible structure
    if 'SeC' in trna_id or 'Sec' in trna_id:
        return True
    # Add other exclusions as needed
    return False
```

#### 3. Update Processing Logic
- **Filter SeC tRNAs** at input processing stage
- **Reserve global_index 65-85** for Type II extended arms
- **Ensure Type I tRNAs** skip reserved space (create gap)
- **Validate no collisions** in final coordinate assignments

### Expected Coordinate Layout
```
Type I:  1...64-[GAP:65-85]-86...102  (76 total positions)
Type II: 1...64-[e1-e24:65-85]-86...102  (91 total positions)
```

### Validation Requirements
1. **No global_index collisions** between any positions
2. **Type I tRNAs** have gap in 65-85 range
3. **Type II "e" positions** map to 65-85 range only
4. **Downstream positions** (T-stem, T-loop, acceptor) align correctly
5. **Cross-tRNA heatmaps** show meaningful structural alignment

## Technical Context

### Current File Structure
- `scripts/trnas_in_space.py` - Main processing script
- `outputs/*_global_coords.tsv` - Generated coordinate files
- Recent commit: "Fix extended variable arm global_index assignment bug"

### Key Functions to Modify
- `sort_key()` - Position sorting logic
- `collect_rows_from_json()` - Data processing
- Add filtering logic for incompatible tRNA types

## Success Criteria

### 1. Coordinate Integrity
- ✅ Zero global_index collisions
- ✅ Contiguous numbering within each tRNA type
- ✅ Reserved space properly utilized

### 2. Structural Alignment
- ✅ Type I/II anticodon loops align
- ✅ Type I/II T-stems align (after extended arm)
- ✅ Type I/II acceptor stems align

### 3. Visualization Compatibility
- ✅ Heatmaps show meaningful cross-tRNA patterns
- ✅ Extended variable arms appear in designated region
- ✅ No massive gaps or misalignment artifacts

## Testing Strategy

1. **Generate test coordinates** for representative tRNAs:
   - Type I: Ala-GGC, Phe-GAA, Gly-GCC
   - Type II: Leu-CAA, Ser-GCT, Tyr-GTA
   - Verify: No SeC included

2. **Validate coordinate assignments**:
   - Check global_index uniqueness
   - Verify extended arm positioning
   - Confirm downstream realignment

3. **Generate test heatmap** to verify visual alignment

## Background Context

This analysis comes from debugging broken tRNA modification heatmaps where:
- Extended variable arms appeared as isolated sections
- Cross-tRNA alignment was completely broken
- Multiple coordinate "fixes" made the problem worse
- Only exclusion-based filtering produced usable visualizations

The root cause is that the global coordinate approach has fundamental structural assumptions that don't accommodate the diversity of tRNA architectures, particularly selenocysteine tRNAs.

## Immediate Actions Needed

1. **Implement Type I/II coordinate system** (exclude SeC)
2. **Test with E. coli dataset** first
3. **Regenerate coordinate files** for all species
4. **Validate heatmap functionality**
5. **Document SeC exclusion** and rationale
6. **Consider separate SeC analysis pipeline** for future work

---

**Prompt for new Claude instance:**

"Please implement a Type I + Type II compatible coordinate system for tRNAs-in-space that excludes selenocysteine tRNAs and properly handles extended variable arms without creating global_index collisions. Focus on modifying the sort_key() function and adding appropriate filtering logic. The goal is to enable meaningful cross-tRNA structural alignment for standard and extended variable arm tRNAs while acknowledging that selenocysteine represents a fundamentally incompatible structural class."