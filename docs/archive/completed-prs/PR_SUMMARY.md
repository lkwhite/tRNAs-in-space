# PR Summary: Infrastructure for Tier 1 Model Organisms Pre-calculation

## Overview

This PR adds the complete infrastructure needed to expand pre-calculated tRNA global coordinates from 3 organisms to 14 organisms (adding 11 Tier 1 model organisms).

## What's Included

### 1. Organism Configuration (`config/organisms.yaml`)
- Complete metadata for all 14 organisms (3 existing + 11 new)
- Assembly versions, taxonomy IDs, expected tRNA counts
- Quality filter requirements (high confidence mature tRNAs)
- Processing pipeline configuration

**New Tier 1 Organisms:**
1. Mus musculus (mouse) - mm39
2. Drosophila melanogaster (fruit fly) - dm6
3. Caenorhabditis elegans (nematode) - ce11
4. Danio rerio (zebrafish) - danRer11
5. Arabidopsis thaliana (thale cress) - araTha1
6. Schizosaccharomyces pombe (fission yeast) - schiPomb972h
7. Bacillus subtilis - bacSubt168
8. Pseudomonas aeruginosa - pseAerPAO1
9. Candida albicans - canAlbSC5314
10. Xenopus tropicalis (frog) - xenTro9
11. Rattus norvegicus (rat) - rn7

### 2. Automated Processing Pipeline (`scripts/process_organisms.py`)
- Batch processing of multiple organisms
- Integrated R2DT annotation + coordinate generation + validation
- Detailed logging and error handling
- Flexible: can process all pending or specific organisms
- Dry-run mode for testing

**Usage:**
```bash
# Process all pending organisms
python scripts/process_organisms.py

# Process specific organisms
python scripts/process_organisms.py --organisms mm39,dm6,ce11

# Skip R2DT if JSONs exist
python scripts/process_organisms.py --skip-r2dt

# Dry run
python scripts/process_organisms.py --dry-run
```

### 3. Download Helper (`scripts/download_gtrnadb_fastas.py`)
- Attempts automated download (currently blocked by GtRNAdb 403 errors)
- Provides exact URLs for manual download
- Validates downloaded files
- Tracks download progress

**Usage:**
```bash
# Show download URLs
python scripts/download_gtrnadb_fastas.py --manual

# Attempt automated download
python scripts/download_gtrnadb_fastas.py
```

### 4. Comprehensive Download Guide (`docs/DOWNLOAD_GUIDE.md`)
- Step-by-step manual download instructions
- High confidence filter criteria (scores ≥ 50/70/10 bits)
- Checklist for all 11 organisms with exact filenames
- Troubleshooting and verification steps
- Expected tRNA counts for validation

## Data Quality Requirements

All new organisms must use **high confidence mature tRNAs** from GtRNAdb:
- Overall score ≥ 50 bits
- Isotype-specific model score ≥ 70 bits (eukaryotes)
- Secondary structure score ≥ 10 bits
- Exclude pseudogenes

This ensures consistency with existing datasets and filters out low-confidence predictions.

## Files Added

```
config/
  └── organisms.yaml              (234 lines) - Organism configuration and metadata

docs/
  └── DOWNLOAD_GUIDE.md          (213 lines) - Manual download instructions

scripts/
  ├── process_organisms.py       (373 lines) - Automated bulk processing
  └── download_gtrnadb_fastas.py (268 lines) - Download helper tool
```

## Commits

1. **fd4c83d** - Add infrastructure for processing Tier 1 model organisms
   - Created organism configuration with all metadata
   - Built automated processing pipeline
   - Added download helper script
   - Initial download guide

2. **3542eca** - Specify high confidence mature tRNA requirement
   - Clarified quality filter requirements
   - Updated download guide with filter criteria
   - Documented tRNAscan-SE 2.0 confidence levels

## Next Steps (Future PRs)

The infrastructure is complete and ready. To finish the expansion:

1. **Manual download** - Download 11 FASTA files from GtRNAdb following `docs/DOWNLOAD_GUIDE.md`
2. **Process organisms** - Run `python scripts/process_organisms.py`
3. **Update metadata** - Update `outputs/METADATA.json` with new organisms
4. **Update documentation** - Update README, FAQ with new organism list
5. **Create example** - Add cross-species comparison notebook

## Testing

All scripts have been validated:
- ✅ `process_organisms.py` - Dry run tested
- ✅ `download_gtrnadb_fastas.py` - Manual mode tested
- ✅ Configuration file - Valid YAML syntax
- ✅ Download guide - Complete and accurate

## Impact

- **User benefit**: Access to pre-calculated coordinates for 11 major model organisms
- **Scalability**: Easy to add more organisms using the same infrastructure
- **Maintainability**: Well-documented, automated pipeline
- **Consistency**: Standardized quality criteria across all organisms

## Breaking Changes

None. This is purely additive infrastructure.

## Documentation

All new functionality is fully documented:
- Inline code comments
- Comprehensive usage examples
- Detailed download guide
- Configuration file is self-documenting

---

**Ready for review and merge!** Once merged, the manual download process can begin.
