# Modomics Integration Scripts

This directory contains scripts for integrating Modomics tRNA modification data with tRNAs-in-space global coordinate system.

## Overview

Modomics (https://iimcb.genesilico.pl/modomics/) provides experimentally validated tRNA modifications from mass spectrometry. These scripts enable mapping those modifications to Sprinzl positions in our structural coordinate system.

## Modules

### `modification_codes.py`
Decoder for Modomics single-character modification codes.

**Usage:**
```bash
python modification_codes.py path/to/modomicscodes.csv
```

**Features:**
- Bidirectional lookup: code ↔ modification name
- Get reference nucleobase for each modification
- Statistics about modification types

**Example:**
```python
from scripts.modomics.modification_codes import ModificationCodec

codec = ModificationCodec('docs/development/modomicscodes.csv')
info = codec.decode('D')  # dihydrouridine
print(info['name'])  # Full name
print(info['short_name'])  # Short name (e.g., "D")
print(info['reference_base'])  # Base nucleotide (U)
```

### `parse_modomics.py`
Parser for Modomics FASTA files (modified and unmodified versions).

**Usage:**
```bash
python parse_modomics.py \
    --modified docs/development/modified_tRNA_all_all_rna_sequences.fasta \
    --unmodified docs/development/unmodified_tRNA_all_all_rna_sequences.fasta \
    --codes docs/development/modomicscodes.csv \
    --output outputs/modomics/modomics_modifications.json
```

**Features:**
- Parse Modomics FASTA headers (species, tRNA type, anticodon, etc.)
- Compare modified/unmodified sequences to detect modification positions
- Decode modification codes using ModificationCodec
- Export structured JSON with all metadata and modifications

**Output format:**
```json
{
  "1": {
    "modomics_id": 1,
    "name": "tdbR00000010",
    "species": "Escherichia coli",
    "subtype": "Ala",
    "anticodon": "VGC",
    "modified_sequence": "GGGGCUAUAGCUCAGCDGGG...",
    "unmodified_sequence": "GGGGCUAUAGCUCAGCUGGG...",
    "modifications": [
      {
        "position": 17,
        "modified_char": "D",
        "unmodified_char": "U",
        "modification_name": "dihydrouridine",
        "short_name": "D",
        "reference_base": "U"
      }
    ]
  }
}
```

### `align_to_sprinzl.py` (TODO - Phase 2)
Aligns Modomics sequences to gtRNAdb sequences and maps modifications to Sprinzl positions.

Will implement:
- Species name normalization
- Sequence alignment (BioPython)
- Position transfer through alignment
- Quality control and validation

## Data Flow

```
modomicscodes.csv ──────────► ModificationCodec
                                      │
                                      ▼
modified_FASTA ──────────┐       decode codes
unmodified_FASTA ────────┤            │
                         ▼            ▼
                    ModomicsParser ───┘
                         │
                         ▼
           modomics_modifications.json
                         │
                         ▼
                [Phase 2: Alignment]
                         │
                         ▼
          modomics_to_sprinzl_mapping.tsv
```

## Configuration

See `config/modomics.yaml` for configuration options including:
- File paths
- Species name mappings
- Alignment parameters
- Quality thresholds

## Testing

Test the modification codec:
```bash
python modification_codes.py docs/development/modomicscodes.csv
```

Expected output: Statistics and example lookups

## Dependencies

- Python 3.8+
- Standard library: csv, json, pathlib, re, logging, dataclasses
- Future (Phase 2): biopython

## Documentation

See `docs/development/MODOMICS_INTEGRATION.md` for:
- Overall implementation plan
- Phase descriptions
- How to continue in next session
- Data flow diagrams

## Status

- ✅ Phase 1: Parser and decoder (COMPLETE)
- ⏳ Phase 2: Alignment to Sprinzl (TODO)
- ⏳ Phase 3: Automated pipeline (TODO)
- ⏳ Phase 4: Analysis tools (TODO)
