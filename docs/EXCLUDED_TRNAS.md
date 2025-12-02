# Excluded tRNAs

This document lists tRNAs excluded from the unified coordinate system due to poor R2DT annotation quality or incompatible structures.

## Overview

tRNAs are excluded for two main reasons:

1. **Structural incompatibility**: tRNAs with non-standard structures that can't align to the canonical coordinate system (mitochondrial, selenocysteine, initiator Met)

2. **Annotation quality**: tRNAs with R2DT annotation errors that prevent reliable alignment (missing labels, wrong anticodon, invalid T-loop)

---

## Structurally Incompatible tRNAs

These are excluded automatically based on tRNA type:

### Mitochondrial tRNAs (`mito-tRNA-*`)

Mitochondrial tRNAs have fundamentally different architecture:
- Often 60-75 nucleotides vs. 76 for nuclear tRNAs
- Can lack certain structural features (e.g., D-loop)
- Different Sprinzl numbering patterns
- Designed for mitochondrial translation system

### Selenocysteine tRNAs (`*SeC*`)

Selenocysteine tRNAs have unusual structure:
- ~95 nucleotides (vs. 76 standard)
- Extended variable arm (17+ nucleotides)
- Unique requirements for stop codon suppression
- Cannot align with Type I/II tRNA coordinates

### Initiator Methionine tRNAs (`*iMet*`, `*fMet*`)

Initiator tRNAs have modified structures:
- Special features for ribosome P-site binding
- Different base pairing patterns
- Cannot be meaningfully compared to elongator tRNAs

---

## Poorly Annotated tRNAs

These are excluded due to R2DT annotation issues. They are listed in `EXCLUDED_POORLY_ANNOTATED` in `scripts/trnas_in_space.py`.

### Current Exclusion List

| tRNA | Organism | Issue | Details |
|------|----------|-------|---------|
| nuc-tRNA-Leu-CAA-5-1 | Human | Missing labels | 59.5% empty sprinzl_labels, anticodon positions (34-36) unlabeled |
| nuc-tRNA-Arg-CCG-1-1 | Yeast | Missing labels | 20.8% empty sprinzl_labels, anticodon positions (34-36) unlabeled |
| nuc-tRNA-Tyr-AUA-1-1 | Human | Wrong anticodon | R2DT annotated anticodon as GGT instead of ATA |
| nuc-tRNA-Gln-UUG-4-1 | Human | Invalid T-loop | T-loop sequence is CGA (not *TC pattern) |
| nuc-tRNA-Trp-CCA-4-1 | Human | Invalid T-loop | T-loop sequence is GCG (not *TC pattern) |
| nuc-tRNA-Trp-CCA-5-1 | Human | Invalid T-loop | T-loop sequence is GCG (not *TC pattern) |

### Exclusion Categories

#### Missing Anticodon Labels

tRNAs where positions 34-35-36 lack sprinzl_labels cannot have their identity verified. Without knowing where the anticodon is, we can't confirm the alignment is correct.

**Affected**: Leu-CAA-5-1 (human), Arg-CCG-1-1 (yeast)

#### Wrong Anticodon Annotation

tRNAs where the anticodon sequence at positions 34-35-36 doesn't match the tRNA name. This indicates R2DT assigned wrong sprinzl_labels to the source sequence.

**Affected**: Tyr-AUA-1-1 (human) - has GGT instead of ATA

#### Invalid T-loop Pattern

tRNAs where positions 54-55-56 don't contain a valid T-loop sequence. Valid patterns are:
- TTC/UUC (canonical)
- TTT/UUU (variant)
- *TC patterns like CTC, ATC, GTC (biological variants)

Patterns like GCG or CGA indicate shifted or incorrect annotations.

**Affected**: Gln-UUG-4-1, Trp-CCA-4-1, Trp-CCA-5-1 (all human)

---

## Impact on tRNA Counts

With exclusions applied:

| Organism | Total in R2DT | Structurally Excluded | Annotation Excluded | Final Count |
|----------|---------------|----------------------|---------------------|-------------|
| E. coli K12 | 86+ | ~4 (iMet) | 0 | 82 |
| S. cerevisiae | 290+ | ~22 (mito, iMet) | 1 | 267 |
| H. sapiens | 480+ | ~35 (mito, SeC, iMet) | 5 | 416 |

---

## How to Add New Exclusions

If you discover a tRNA with annotation issues:

### 1. Identify the Issue

Run biological validation tests:
```bash
python -m pytest test_trnas_in_space.py::test_anticodon_matches_trna_name -v
python -m pytest test_trnas_in_space.py::test_tloop_contains_ttc -v
```

### 2. Verify It's an R2DT Issue

Check the raw JSON to confirm the issue is in R2DT output, not our processing:
```python
import json
with open("outputs/hg38_jsons/nuc-tRNA-XXX.enriched.json") as f:
    data = json.load(f)
mol = data["rnaComplexes"][0]["rnaMolecules"][0]
for item in mol["sequence"]:
    label = item["info"].get("templateNumberingLabel", "")
    if label in ["34", "35", "36"]:  # anticodon
        print(f"pos {label}: {item['residueName']}")
```

### 3. Add to Exclusion List

Edit `scripts/trnas_in_space.py`:
```python
EXCLUDED_POORLY_ANNOTATED = {
    # ... existing entries ...
    "nuc-tRNA-XXX-1-1",   # reason for exclusion
}
```

### 4. Regenerate Coordinates

```bash
python scripts/trnas_in_space.py outputs/hg38_jsons outputs/hg38_global_coords.tsv
```

### 5. Update This Documentation

Add the new exclusion to the table above with:
- tRNA name
- Organism
- Issue category
- Specific details

---

## T-loop Variants (NOT Excluded)

Some tRNAs have non-canonical T-loop sequences that are valid biological variants:

| tRNA | T-loop | Status |
|------|--------|--------|
| nuc-tRNA-Ile-GAU-1-1/1-2/1-3 | CTC | Included (valid *TC) |
| nuc-tRNA-Gly-UCC-4-1 | CTC | Included (valid *TC) |
| nuc-tRNA-Lys-CUU-8-1 | CTC | Included (valid *TC) |
| nuc-tRNA-Val-AAC-3-1 | ATC | Included (valid *TC) |
| nuc-tRNA-Val-UAC-4-1 | GTC | Included (valid *TC) |

These tRNAs:
- Have correct anticodons (verified against name)
- Have proper Sprinzl position assignments
- Simply have variant T-loop sequences (CTC, ATC, GTC instead of TTC)

The *TC pattern (any nucleotide + TC) is accepted as a valid T-loop variant.

---

## Related Documentation

- [BIOLOGICAL_VALIDATION.md](BIOLOGICAL_VALIDATION.md) - Validation test details
- [COORDINATE_SYSTEM.md](COORDINATE_SYSTEM.md) - How coordinates are built
- [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) - Output file structure
