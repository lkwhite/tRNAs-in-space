# Biological Validation

This document describes the biological validation tests used to verify the quality and correctness of tRNA coordinate assignments.

## Overview

The coordinate system relies on R2DT's Sprinzl position annotations. To catch annotation errors that could produce incorrect alignments, we validate two highly conserved tRNA features:

1. **Anticodon** (positions 34-35-36): Must match the tRNA's identity
2. **T-loop** (positions 54-55-56): Must contain a valid TψC or variant sequence

These tests run automatically via pytest and enforce 100% pass rates for all included tRNAs.

---

## Anticodon Validation

### What It Checks

The anticodon at positions 34-35-36 defines which amino acid a tRNA carries. The nucleotides at these positions must match the anticodon in the tRNA name.

**Example:**
- `tRNA-Ala-AGC-1-1` must have A-G-C at positions 34-35-36
- `tRNA-Arg-CCU-1-1` must have C-C-T (T=U equivalent) at positions 34-35-36

### Why It Matters

If the anticodon doesn't match, it indicates one of:
1. **Shifted labels**: R2DT assigned wrong Sprinzl positions
2. **Annotation error**: Wrong anticodon in R2DT source data
3. **Processing bug**: Our code incorrectly modified labels

Either way, the tRNA cannot be correctly aligned and must be excluded.

### Implementation

```python
# From test_trnas_in_space.py
def test_anticodon_matches_trna_name():
    """Verify positions 34-35-36 contain the anticodon from the tRNA name."""
    # Extract expected anticodon from name (e.g., Ala-AGC -> AGC)
    # Check positions 34-35-36 match
    # Allow T/U equivalence
```

---

## T-loop Validation

### What It Checks

The T-loop (TψC loop) at positions 54-55-56 is one of the most conserved features in tRNA. It contains a characteristic T-Ψ-C sequence (shown as T-T-C in genomic/DNA notation since Ψ is pseudouridine, a modified U).

### Valid Patterns

| Pattern | Type | Prevalence |
|---------|------|------------|
| TTC / UUC | Canonical | ~97% of tRNAs |
| TTT / UUU | Common variant | ~1% of tRNAs |
| CTC, ATC, GTC | *TC variants | ~2% of tRNAs (7 human) |

### T-loop Variants Discovery

During validation, we discovered that 7 human tRNAs have non-canonical but valid *TC patterns at the T-loop:

| tRNA | T-loop | Notes |
|------|--------|-------|
| nuc-tRNA-Ile-GAU-1-1 | CTC | All 3 Ile-GAU copies |
| nuc-tRNA-Ile-GAU-1-2 | CTC | |
| nuc-tRNA-Ile-GAU-1-3 | CTC | |
| nuc-tRNA-Gly-UCC-4-1 | CTC | |
| nuc-tRNA-Lys-CUU-8-1 | CTC | |
| nuc-tRNA-Val-AAC-3-1 | ATC | |
| nuc-tRNA-Val-UAC-4-1 | GTC | |

These tRNAs:
- Have **correct anticodons** (verified against tRNA name)
- Are **correctly aligned** (Sprinzl positions match expected structure)
- Simply have a **variant T-loop sequence** (biological variation, not error)

### Invalid Patterns

Patterns that don't end in TC indicate annotation errors:

| Pattern | Indicates |
|---------|-----------|
| GCG | Shifted annotation (Trp-CCA-4-1, Trp-CCA-5-1) |
| CGA | Shifted annotation (Gln-UUG-4-1) |

tRNAs with invalid T-loop patterns are excluded from the coordinate system.

### Implementation

```python
# From test_trnas_in_space.py
def test_tloop_contains_ttc():
    """Verify positions 54-55-56 contain valid T-loop sequence."""
    # Valid: TTC, UUC, TTT, UUU, or *TC patterns
    valid_tloop = (
        tloop in ("TTC", "UUC", "TTU", "UUU", "TTT") or
        tloop.endswith("TC") or tloop.endswith("UC")
    )
```

---

## Why Strict Validation Matters

### Previous Approach (Thresholds)

Early versions used percentage thresholds (e.g., "pass if >90% of tRNAs have correct anticodons"). This masked:
- Individual tRNAs with annotation errors
- Systematic bugs affecting subsets of tRNAs

### Current Approach (100% Pass Rate)

Every tRNA in the coordinate files must pass both validations:

1. **Anticodon matches name** (positions 34-35-36)
2. **T-loop is valid** (positions 54-55-56: TTC/TTT or *TC variant)

tRNAs that fail are added to `EXCLUDED_POORLY_ANNOTATED` and excluded from coordinate generation.

### Benefits

- **Data integrity**: All included tRNAs are verified correct
- **Early detection**: Bugs show as test failures, not subtle misalignments
- **Clear exclusions**: Problematic tRNAs are documented with reasons

---

## Historical Context

### The fix_label_index_mismatch Bug

An early attempt to fix R2DT annotation issues introduced a bug:

1. **Observed**: Some tRNAs had `sprinzl_index` > `sprinzl_label` at certain positions
2. **Hypothesis**: R2DT failed to skip labels when indices skipped (deletions)
3. **"Fix"**: Implemented `fix_label_index_mismatch()` to shift labels

**The problem**: R2DT's labels were actually correct. The index/label divergence was intentional for tRNAs with structural deletions. The "fix" broke 21 yeast tRNAs by shifting their anticodons.

**Discovery**: Anticodon validation tests showed shifted sequences:
- Arg-CCU → TCC (shifted +1)
- Leu-CAA → TCA (shifted +1)

**Resolution**: Removed `fix_label_index_mismatch()`, added biological validation tests.

---

## Adding New Validations

To add additional biological checks:

1. Add test function in `test_trnas_in_space.py`
2. Use 100% pass rate (no thresholds)
3. Add failing tRNAs to `EXCLUDED_POORLY_ANNOTATED` in `trnas_in_space.py`
4. Document in this file

### Potential Future Validations

- **Acceptor stem**: Positions 1-7 should base-pair with 66-72
- **D-loop conserved G**: Position 18 is usually G
- **Variable region length**: Type I vs Type II classification

---

## Related Documentation

- [EXCLUDED_TRNAS.md](EXCLUDED_TRNAS.md) - List of excluded tRNAs with reasons
- [COORDINATE_SYSTEM.md](COORDINATE_SYSTEM.md) - How coordinates are built
- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) - Output file structure
