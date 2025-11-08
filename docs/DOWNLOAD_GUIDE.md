# GtRNAdb FASTA Download Guide

GtRNAdb blocks automated downloads with 403 Forbidden errors. You'll need to manually download FASTA files for each organism.

## Quick Start

For each organism in Tier 1, follow these steps:

### Method 1: Direct Download from GtRNAdb (Recommended)

1. **Visit GtRNAdb**: http://gtrnadb.ucsc.edu/
2. **Navigate to organism**:
   - Click on the appropriate kingdom (Eukaryota, Bacteria, etc.)
   - Find and click on your organism
3. **Download FASTA**:
   - Look for "Download" or "Sequences" section
   - Download the mature tRNA sequences in FASTA format
   - For eukaryotes, get both nuclear and mitochondrial if available
4. **Save to fastas/ directory** with the naming convention below

### Method 2: UCSC Table Browser

Some organisms may be available through the UCSC Table Browser:

1. Visit: https://genome.ucsc.edu/cgi-bin/hgTables
2. Select your organism and assembly
3. Group: "Genes and Gene Predictions"
4. Track: "tRNAs" or search for GtRNAdb
5. Output format: "sequence"
6. Get output

## Tier 1 Organisms - Download Checklist

Save all files to the `fastas/` directory with these exact filenames:

### Mammals
- [ ] **Mouse (Mus musculus)** - Assembly: mm39/GRCm39
  - Filename: `mm39-tRNAs.fa`
  - GtRNAdb: Eukaryota → Mus musculus → mm39

- [ ] **Rat (Rattus norvegicus)** - Assembly: rn7/mRatBN7.2
  - Filename: `rn7-tRNAs.fa`
  - GtRNAdb: Eukaryota → Rattus norvegicus → rn7

### Model Organisms
- [ ] **Fruit fly (Drosophila melanogaster)** - Assembly: dm6/BDGP6
  - Filename: `dm6-tRNAs.fa`
  - GtRNAdb: Eukaryota → Drosophila melanogaster → dm6

- [ ] **Nematode (C. elegans)** - Assembly: ce11/WBcel235
  - Filename: `ce11-tRNAs.fa`
  - GtRNAdb: Eukaryota → Caenorhabditis elegans → ce11

- [ ] **Zebrafish (Danio rerio)** - Assembly: danRer11/GRCz11
  - Filename: `danRer11-tRNAs.fa`
  - GtRNAdb: Eukaryota → Danio rerio → danRer11

- [ ] **Frog (Xenopus tropicalis)** - Assembly: xenTro9
  - Filename: `xenTro9-tRNAs.fa`
  - GtRNAdb: Eukaryota → Xenopus tropicalis → xenTro9

### Yeast and Fungi
- [ ] **Fission yeast (S. pombe)** - Assembly: ASM294v2
  - Filename: `schiPomb972h-tRNAs.fa`
  - GtRNAdb: Eukaryota → Schizosaccharomyces pombe → schiPomb972h

- [ ] **C. albicans** - Assembly: SC5314
  - Filename: `canAlbSC5314-tRNAs.fa`
  - GtRNAdb: Eukaryota → Candida albicans → canAlbSC5314

### Plant
- [ ] **Arabidopsis (A. thaliana)** - Assembly: TAIR10
  - Filename: `araTha1-tRNAs.fa`
  - GtRNAdb: Eukaryota → Arabidopsis thaliana → araTha1

### Bacteria
- [ ] **B. subtilis** - Assembly: 168
  - Filename: `bacSubt168-tRNAs.fa`
  - GtRNAdb: Bacteria → Bacillus subtilis → 168

- [ ] **P. aeruginosa** - Assembly: PAO1
  - Filename: `pseAerPAO1-tRNAs.fa`
  - GtRNAdb: Bacteria → Pseudomonas aeruginosa → PAO1

## File Naming Convention

The naming pattern is:
```
{gtrnadb_id}-tRNAs.fa
```

For organisms with mitochondrial tRNAs, you may see:
```
{gtrnadb_id}-mito-and-nuclear-tRNAs.fa
```

Either naming is fine - the processing script will detect both.

## Verification

After downloading, verify the files:

```bash
# Check that files are present
ls -lh fastas/

# Verify they are FASTA format (should start with '>')
head -n 1 fastas/mm39-tRNAs.fa

# Count sequences (should match expected tRNA counts)
grep -c "^>" fastas/mm39-tRNAs.fa
```

## Troubleshooting

### File not found on GtRNAdb
- Try searching for an older assembly version (e.g., mm10 instead of mm39)
- Check if the organism is listed under a different name
- Some assemblies may not be available yet - note this in issues

### Wrong file format
- Make sure you're downloading "mature tRNA sequences"
- Avoid "tRNA genes" with genomic context
- FASTA format should have header lines starting with '>' followed by sequence

### Expected sequence counts

Approximate expected tRNA counts (for validation):

| Organism | Expected tRNAs |
|----------|----------------|
| Mouse (mm39) | 400-500 |
| Rat (rn7) | 400-500 |
| Fly (dm6) | 280-320 |
| Worm (ce11) | 590-630 |
| Zebrafish (danRer11) | 500-600 |
| Frog (xenTro9) | 400-500 |
| S. pombe | 150-180 |
| C. albicans | 90-110 |
| Arabidopsis | 600-700 |
| B. subtilis | 80-90 |
| P. aeruginosa | 60-70 |

If your counts are very different, double-check you have the right file.

## Alternative: tRNAscan-SE

If a FASTA file is not available on GtRNAdb, you can run tRNAscan-SE yourself:

```bash
# Install tRNAscan-SE (requires separate installation)
# See: http://lowelab.ucsc.edu/tRNAscan-SE/

# Run on genome FASTA
tRNAscan-SE -o output.txt -f output.ss -m stats.txt genome.fa

# Extract mature tRNA sequences (you'll need to write a custom script)
```

This is more advanced and should only be used if GtRNAdb doesn't have your organism.

## Next Steps

Once you have downloaded the FASTA files to `fastas/`, run:

```bash
# Process all organisms
python scripts/process_organisms.py

# Or process specific organisms
python scripts/process_organisms.py --organisms mm39,dm6,ce11
```

## Need Help?

If you encounter issues:
1. Check that the file is valid FASTA format
2. Verify the sequence count is reasonable
3. Check file size (should be >10 KB for most organisms)
4. Open an issue with details about which organism failed
