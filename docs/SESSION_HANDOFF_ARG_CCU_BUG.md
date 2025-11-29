# Session Handoff: R2DT Label/Index Mismatch Bug

## Quick Context

We discovered that yeast tRNA-Arg-CCU-1-1 was misaligned in heatmaps. Investigation revealed this is a **systemic R2DT bug affecting ~64 tRNAs** across all three organisms.

## The Bug in One Sentence

When R2DT detects a structural deletion (gap in `sprinzl_index`), it correctly skips the index but **fails to skip the corresponding `sprinzl_label`**, causing alignment errors.

## Example

```
# BUG: index skips 21→22, but label doesn't
seq=20, idx=20, label="20"  ← OK
seq=21, idx=22, label="21"  ← WRONG! label should be "22"

# CORRECT behavior would be:
seq=20, idx=20, label="20"
seq=21, idx=22, label="22"  ← label also skips
```

## Current State

1. **Manual fixes applied** for 3 tRNAs via `LABEL_OVERRIDES` in `scripts/trnas_in_space.py`:
   - `nuc-tRNA-Arg-CCU-1-1` (yeast)
   - `nuc-tRNA-Leu-CAA-6-1` (human)
   - `nuc-tRNA-His-GUG-1-*` (human, 9 isodecoders)

2. **Tests added** in `test_trnas_in_space.py`:
   - Detection test (currently warn-only, doesn't fail)
   - Regression test for Arg-CCU fix

3. **Documentation** at `docs/R2DT_LABEL_INDEX_MISMATCH_BUG.md`

## What Needs to Be Done

**Option 1 (Recommended): Automated Fix**
- Detect the pattern during `collect_rows_from_json()`
- Auto-correct labels when `sprinzl_index` skips but `sprinzl_label` doesn't
- Remove manual `LABEL_OVERRIDES`

**Concerns to address:**
- Handle non-zero offset tRNAs (where label = index + N consistently)
- Handle cascading corrections (fixing one position may require fixing all subsequent)
- Don't break tRNAs that are currently working
- Handle non-numeric labels ("20a", "e5")

## Key Files

- `scripts/trnas_in_space.py` - Main pipeline, lines 30-73 have `LABEL_OVERRIDES`
- `test_trnas_in_space.py` - Tests, lines 456-507 have detection test
- `docs/R2DT_LABEL_INDEX_MISMATCH_BUG.md` - Full documentation

## Commands

```bash
# See all 64 affected tRNAs
python3 -c "
import pandas as pd
from pathlib import Path
for f in Path('outputs').glob('*_offset*.tsv'):
    df = pd.read_csv(f, sep='\t')
    for tid, g in df.groupby('trna_id'):
        g = g.sort_values('seq_index')
        rows = g.to_dict('records')
        for i in range(1, len(rows)):
            pi, ci = rows[i-1]['sprinzl_index'], rows[i]['sprinzl_index']
            cl = str(rows[i]['sprinzl_label']).strip()
            if pi > 0 and ci > 0 and ci > pi + 1 and cl.isdigit() and int(cl) in range(pi+1, ci):
                print(f'{tid}: seq={rows[i][\"seq_index\"]} idx={ci} label={cl}')
"

# Run tests
python test_trnas_in_space.py

# Regenerate after changes
python scripts/trnas_in_space.py outputs/sacCer_jsons outputs/sacCer_global_coords.tsv --split-by-offset-and-type
```

## The Fundamental Question

Is it safe to assume: **when `sprinzl_index` jumps (deletion), the `sprinzl_label` should also jump by the same amount**?

If yes → automated fix is straightforward
If no → need to understand exceptions first
