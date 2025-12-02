# DEPRECATED: Offset-Type Coordinate Grouping

**Status**: Deprecated as of 2025-12-01
**Replaced by**: Unified coordinate system (single file per organism)

## What It Was

The offset-type grouping system generated separate coordinate files for combinations of:
- **Offset** (-3, -2, -1, 0, +1): D-loop labeling variation
- **Type** (type1, type2): Standard vs extended variable arm

Example files:
```
hg38_global_coords_offset0_type1.tsv
hg38_global_coords_offset0_type2.tsv
hg38_global_coords_offset-1_type1.tsv
...
```

## Why It Existed

The grouping was implemented as a workaround for two bugs that caused global_index collisions:

### Bug 1: Index Fallback in build_pref_label()

**Problem**: When `sprinzl_label` was empty (insertions without canonical positions), the code fell back to using `sprinzl_index`. This caused insertions to get wrong labels (e.g., "60" instead of empty), which sorted incorrectly.

**Example**:
```
Insertion after position 52 got label "60" (from sprinzl_index)
This sorted AFTER position 52, causing collision with position 60
```

### Bug 2: Mixed Types in Sort Keys

**Problem**: The third element of sort key tuples had mixed types:
- E-positions returned `int` (bio_order)
- Letter suffixes returned `str` ("A", "B")

Python 3 can't compare `int` with `str`, causing `TypeError`.

## Why It's No Longer Needed

Both bugs have been fixed:

1. **build_pref_label()**: No longer falls back to sprinzl_index. Empty labels stay empty and get interpolated coordinates.

2. **sort_key()** (and variants): All third elements are now strings, using zero-padding for numeric values.

With these fixes, the unified coordinate system produces **0 collisions** for all tested organisms:
- E. coli K12: 82 tRNAs, 136 positions
- S. cerevisiae: 268 tRNAs, 145 positions
- H. sapiens: 422 tRNAs, 253 positions

## Migration

If you have code using the old offset-type files:

**Before**:
```python
df = pd.read_csv("hg38_global_coords_offset0_type1.tsv", sep="\t")
```

**After**:
```python
df = pd.read_csv("hg38_global_coords.tsv", sep="\t")
```

All tRNA types are now in a single file with consistent global_index values.

## Removed Code

The following were removed from `trnas_in_space.py`:
- `GroupKey` dataclass
- `GroupingStrategy` base class
- `OffsetTypeStrategy` class
- `generate_grouped_coordinates()` function
- `_generate_coordinates_for_group()` function
- `--split-by-offset-and-type` CLI flag

## Historical Context

The offset-type grouping was added on 2025-11-25 to work around collisions that appeared when combining all tRNAs. At the time, the root cause (the two bugs above) was not identified.

Investigation in late November 2025 revealed:
1. The e-position ordering bug was fixed first
2. Then the empty label violation pattern was discovered
3. Finally, the index fallback bug was identified as the root cause

With the bugs fixed, offset-type grouping became unnecessary and was deprecated.

## Related Documentation

- [COORDINATE_SYSTEM.md](../COORDINATE_SYSTEM.md) - Current unified system
- [OUTPUT_FORMAT.md](../OUTPUT_FORMAT.md) - Output file format
- Plan file: `~/.claude/plans/parallel-squishing-book.md` (development history)
