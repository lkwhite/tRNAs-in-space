# Modomics Integration - Implementation Guide

**Status:** Phase 1 Setup Complete
**Next Step:** Continue with unmodified FASTA processing

---

## Overview

This document tracks the implementation of Modomics tRNA modification data integration into tRNAs-in-space. The goal is to map known modifications from Modomics mass spec data to Sprinzl positions in our global coordinate system.

## Use Cases Enabled

1. **Structure-aware analysis** - Analyze modifications by tRNA structural region (acceptor-stem, D-loop, anticodon-loop, etc.)
2. **Position-specific patterns** - Extract all modifications at specific Sprinzl positions (e.g., position 34 wobble)
3. **Cross-tRNA comparisons** - Align modifications across different tRNAs using global indices
4. **Cross-species conservation** - Compare modifications at equivalent Sprinzl positions across species
5. **Detection validation** - Compare detected modifications to known Modomics annotations by structural context

---

## Implementation Phases

### ✅ Phase 1: Modomics Parser & Modification Decoder (COMPLETED)

**What was built:**
- `scripts/modomics/modification_codes.py` - Decoder for Modomics single-character codes
- `scripts/modomics/parse_modomics.py` - Parser for Modomics FASTA files
- Directory structure: `scripts/modomics/` and `outputs/modomics/`

**Files created:**
```
scripts/modomics/
├── __init__.py
├── modification_codes.py    # Maps single-char codes → modification names
└── parse_modomics.py         # Parses Modomics FASTA files

outputs/modomics/             # For output files
```

**Capabilities:**
- Parse `modomicscodes.csv` to decode modification characters
- Parse Modomics FASTA headers to extract metadata
- Compare modified vs unmodified sequences to detect modification positions
- Export structured JSON with all modification data

---

### 🔄 Phase 2: Sequence Alignment & Position Mapping (NEXT)

**Goal:** Map Modomics sequences to gtRNAdb sequences to transfer Sprinzl positions

**Required inputs:**
- ✅ `docs/development/modified_tRNA_all_all_rna_sequences.fasta` (already present)
- ⏳ `docs/development/unmodified_tRNA_all_all_rna_sequences.fasta` (YOU NEED TO ADD THIS)
- ✅ `docs/development/modomicscodes.csv` (already present)
- ✅ `outputs/{organism}_global_coords.tsv` files (already present for 3 species)

**Implementation tasks:**
1. Run Phase 1 parser to generate `modomics_modifications.json`
2. Build species name normalization (Modomics ↔ gtRNAdb)
3. Implement tRNA sequence alignment (BioPython)
4. Map Modomics positions → gtRNAdb positions → Sprinzl positions
5. Generate `modomics_to_sprinzl_mapping.tsv`

**File to create:**
- `scripts/modomics/align_to_sprinzl.py` - Alignment and mapping logic

---

### 📋 Phase 3: Automated Pipeline (FUTURE)

**Goal:** Integrate Modomics mapping into `process_organisms.py`

**Tasks:**
- Add `--with-modomics` flag to organism processing
- Auto-generate species-specific Modomics annotation files
- Validation and QC reporting

---

### 📊 Phase 4: Analysis Tools (FUTURE)

**Goal:** Provide query functions for common analysis patterns

**Tasks:**
- Query by structural region
- Query by Sprinzl position
- Cross-species comparisons
- Detection validation helpers

---

## How to Continue This Work (NEXT SESSION)

### Prerequisites

Before starting the next session, ensure you have:
1. ✅ This repository checked out on branch `claude/modomics-species-mapping-011CV4PiWXABfn89mspa1222`
2. ⏳ **ADD THIS FILE:** `docs/development/unmodified_tRNA_all_all_rna_sequences.fasta`
   - Should have same format as modified version
   - Same headers, but sequences contain only A, C, G, U (no modification codes)

### Starting the Next Session

**Option 1: Quick continuation (if unmodified FASTA is ready)**

Use this prompt:
```
Continue the Modomics integration work from docs/development/MODOMICS_INTEGRATION.md.
The unmodified FASTA file is now available at docs/development/unmodified_tRNA_all_all_rna_sequences.fasta.

Please:
1. Run the parser to generate modomics_modifications.json
2. Start Phase 2: Build the alignment and Sprinzl mapping pipeline

Focus on E. coli first as a test case.
```

**Option 2: Test Phase 1 first**

Use this prompt:
```
Test the Modomics Phase 1 parser that was set up in the previous session.

The unmodified FASTA is at: docs/development/unmodified_tRNA_all_all_rna_sequences.fasta

Run:
python scripts/modomics/parse_modomics.py \
    --modified docs/development/modified_tRNA_all_all_rna_sequences.fasta \
    --unmodified docs/development/unmodified_tRNA_all_all_rna_sequences.fasta \
    --codes docs/development/modomicscodes.csv \
    --output outputs/modomics/modomics_modifications.json

Then show me statistics and a few example entries.
```

**Option 3: Build Phase 2 alignment pipeline**

Use this prompt:
```
I want to build Phase 2 of the Modomics integration (sequence alignment and Sprinzl mapping).

Context in: docs/development/MODOMICS_INTEGRATION.md

Start by creating scripts/modomics/align_to_sprinzl.py that:
1. Loads modomics_modifications.json
2. Normalizes species names (Modomics ↔ gtRNAdb)
3. Aligns Modomics unmodified sequences to gtRNAdb sequences
4. Maps modification positions through alignment to Sprinzl positions
5. Outputs modomics_to_sprinzl_mapping.tsv

Focus on E. coli as a test case first.
```

---

## Testing the Setup (Current Session)

You can test what's been built so far:

### Test 1: Modification Code Decoder
```bash
python scripts/modomics/modification_codes.py docs/development/modomicscodes.csv
```

Expected output: Statistics about loaded modifications and example lookups

### Test 2: Parser (once unmodified FASTA is added)
```bash
python scripts/modomics/parse_modomics.py \
    --modified docs/development/modified_tRNA_all_all_rna_sequences.fasta \
    --unmodified docs/development/unmodified_tRNA_all_all_rna_sequences.fasta \
    --codes docs/development/modomicscodes.csv
```

Expected output: `outputs/modomics/modomics_modifications.json` with all parsed data

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                        Phase 1 (DONE)                           │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  modomicscodes.csv ──────────► ModificationCodec               │
│         │                              │                        │
│         │                              ▼                        │
│         │                    [D] → "dihydrouridine"            │
│         │                    [T] → "ribothymidine"             │
│         │                              │                        │
│         ▼                              │                        │
│  modified_tRNA_*.fasta ────────┐       │                        │
│  unmodified_tRNA_*.fasta ──────┤       │                        │
│                                │       │                        │
│                                ▼       ▼                        │
│                        ModomicsParser                           │
│                                │                                │
│                                ▼                                │
│                   modomics_modifications.json                   │
│         (tRNA metadata + modification positions)                │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                     Phase 2 (TODO NEXT)                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  modomics_modifications.json ───────┐                          │
│                                     │                          │
│  gtRNAdb FASTAs (fastas/*.fa) ──────┤                          │
│                                     │                          │
│  global_coords.tsv files ───────────┤                          │
│                                     │                          │
│                                     ▼                          │
│                        Sequence Alignment                      │
│                     (BioPython pairwise2)                      │
│                                     │                          │
│                                     ▼                          │
│              Position Transfer via Alignment                   │
│         (Modomics pos → gtRNAdb pos → Sprinzl)               │
│                                     │                          │
│                                     ▼                          │
│                modomics_to_sprinzl_mapping.tsv                 │
│     (Every modification mapped to Sprinzl position)            │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                   Phase 3 & 4 (FUTURE)                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  • Automated pipeline integration                              │
│  • Species-specific output files                               │
│  • Query/analysis functions                                    │
│  • Cross-species comparison tools                              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## File Locations

### Input Files (Already Present)
- `docs/development/modified_tRNA_all_all_rna_sequences.fasta` - Modomics tRNAs with modifications
- `docs/development/modomicscodes.csv` - Modification code decoder
- `outputs/ecoliK12_global_coords.tsv` - E. coli Sprinzl mappings
- `outputs/sacCer_global_coords.tsv` - Yeast Sprinzl mappings
- `outputs/hg38_global_coords.tsv` - Human Sprinzl mappings

### Input Files (YOU NEED TO ADD)
- `docs/development/unmodified_tRNA_all_all_rna_sequences.fasta` - Clean sequences for alignment

### Output Files (Will Be Generated)
- `outputs/modomics/modomics_modifications.json` - Parsed Modomics data
- `outputs/modomics/modomics_to_sprinzl_mapping.tsv` - Master mapping file
- `outputs/modomics/species_coverage.tsv` - Species overlap report
- `outputs/modomics/alignment_qc.tsv` - Alignment quality metrics
- `outputs/{organism}_modomics_annotations.tsv` - Per-species files

---

## Species Coverage

**Modomics species:** 82 unique species (see full list in FASTA)

**gtRNAdb completed:** 3 species
- Escherichia coli
- Saccharomyces cerevisiae
- Homo sapiens

**gtRNAdb pending:** 12 additional species in `config/organisms.yaml`

**Overlap priorities:**
1. Start with the 3 completed species
2. When new species are processed through the gtRNAdb pipeline, automatically map Modomics data if available

---

## Technical Notes

### Dependencies
- **Existing:** pandas, numpy, json, csv
- **New (needed for Phase 2):** BioPython (`pip install biopython`)

### Alignment Strategy
For Phase 2, we'll use:
- **Algorithm:** Needleman-Wunsch (global alignment) via BioPython's pairwise2
- **Scoring:** Match=2, Mismatch=-1, Gap open=-2, Gap extend=-0.5
- **Quality threshold:** Require >85% identity for confident mapping
- **Handling gaps:** Transfer positions only for aligned regions (no gaps)

### Species Name Normalization
Will need to map between formats:
- Modomics: "Escherichia coli"
- gtRNAdb organism_id: "ecoliK12"
- Strategy: Create manual mapping file + fuzzy matching fallback

---

## Questions to Resolve (Future)

1. **Multiple isoforms:** How to handle when Modomics has 1 tRNA but gtRNAdb has multiple copies?
   - Proposed: Map to all matching isoforms

2. **Cellular localization:** Modomics field is often empty. How to infer?
   - Proposed: Use gtRNAdb annotations + sequence length heuristics

3. **Quality thresholds:** What alignment score is "good enough"?
   - Proposed: Start with >85% identity, flag 80-85% for review

4. **Ambiguous anticodons:** Modomics uses modified bases in anticodons (e.g., "ICG")
   - Proposed: Create conversion rules (I→G, etc.)

---

## Contact / Session Handoff

- **Branch:** `claude/modomics-species-mapping-011CV4PiWXABfn89mspa1222`
- **Key files to review:**
  - This document: `docs/development/MODOMICS_INTEGRATION.md`
  - Parser code: `scripts/modomics/parse_modomics.py`
  - Decoder code: `scripts/modomics/modification_codes.py`

**Ready to continue when you add:** `docs/development/unmodified_tRNA_all_all_rna_sequences.fasta`

---

**Last Updated:** Session 1 - Phase 1 Setup Complete
