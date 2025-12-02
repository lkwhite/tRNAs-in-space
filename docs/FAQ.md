# Frequently Asked Questions

## Using Pre-computed Coordinates

### Which organisms are available?

Currently: E. coli K12, S. cerevisiae, and H. sapiens. See `outputs/` directory.

Each organism has:
- One unified coordinate file for nuclear tRNAs (`{species}_global_coords.tsv`)
- Optionally, a separate file for mitochondrial tRNAs (`{species}_mito_global_coords.tsv`)

### Which coordinate file should I use?

For nuclear tRNAs, use the unified file (e.g., `ecoliK12_global_coords.tsv`). This includes both Type I (standard) and Type II (extended variable arm) tRNAs in a single coordinate space.

For mitochondrial tRNAs, use the mito file (e.g., `hg38_mito_global_coords.tsv`). Mito tRNAs have different structural architecture and require a separate coordinate system.

### What are Type I and Type II tRNAs?

- **Type I**: Standard tRNAs with short variable loops (Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Lys, Met, Phe, Pro, Thr, Trp, Val)
- **Type II**: tRNAs with extended variable arms (Leu, Ser, Tyr) — these have additional e1-e24 positions

Both types are included in the unified coordinate files.

### How do I load the data?

**Python:**
```python
import pandas as pd
# Load nuclear tRNAs
df = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')

# Load mitochondrial tRNAs (if available)
df_mito = pd.read_csv('outputs/hg38_mito_global_coords.tsv', sep='\t')
```

**R:**
```r
library(readr)
df <- read_tsv('outputs/ecoliK12_global_coords.tsv')
```

### What do the columns mean?

See [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md) for detailed descriptions. Key columns:
- `global_index`: The standardized position (use this for alignment)
- `region`: Structural region (acceptor-stem, anticodon-loop, etc.)
- `residue`: The nucleotide base

### Which tRNA is the best tRNA?

Phenylalanine GAA. Obviously.

(It was the [first tRNA structure ever solved](https://doi.org/10.1038/246031a0) by crystallography in 1974, earning a Nobel Prize. Also it's just cool.)

## Generating New Coordinates

### How do I process my own tRNAs?

1. Run R2DT on your FASTA file (requires Docker):
   ```bash
   docker run --rm -v "$(pwd):/data" rnacentral/r2dt \
     r2dt.py gtrnadb draw /data/fastas/my_trnas.fa /data/outputs/my_jsons
   ```

2. Generate coordinates:
   ```bash
   # Nuclear tRNAs
   python scripts/trnas_in_space.py outputs/my_jsons outputs/my_coords.tsv

   # Mitochondrial tRNAs (separate coordinate file)
   python scripts/trnas_in_space.py outputs/my_jsons outputs/my_mito_coords.tsv --mito
   ```

### R2DT fails on my sequences. What should I do?

- Make sure sequences are mature tRNAs (no adapters)
- Check that FASTA headers are simple (no special characters)
- Try with a smaller test set first
- See [R2DT documentation](https://docs.r2dt.bio/)

### Can I use this without Docker?

R2DT requires Docker. If you can't use Docker, you can:
- Use the pre-computed datasets
- Use R2DT's web interface and download the JSON files
- Contact us about processing your sequences

## Understanding the Output

### What is `global_index`?

The standardized position number (1, 2, 3...) that's consistent across all tRNAs within the same coordinate file. Use this as your X-axis for alignments and heatmaps.

### Why are some `global_index` values missing for certain tRNAs?

tRNAs vary in length, especially in the variable loop. Missing values mean that tRNA doesn't have a nucleotide at that position.

### What's the difference between `sprinzl_index` and `global_index`?

- `sprinzl_index`: Traditional Sprinzl numbering (1-76, with gaps and letter suffixes)
- `global_index`: Continuous numbering (1, 2, 3...) for easy plotting

## Mitochondrial tRNAs

### How do I generate mito tRNA coordinates?

Use the `--mito` flag:
```bash
python scripts/trnas_in_space.py outputs/hg38_jsons outputs/hg38_mito_global_coords.tsv --mito
```

Mito tRNAs have separate coordinate files because they have different structural architecture (60-75 nt vs 76 nt for nuclear).

### Are mitochondrial tRNAs accurate?

For yeast mitochondrial tRNAs: **use with caution**. R2DT lacks fungal mitochondrial-specific models, so alignments may need verification. See Reinsch & Garcia 2025 for curated yeast mito alignments.

For human mitochondrial tRNAs: Generally good quality. Some tRNAs have R2DT label shifts that are automatically corrected.

### How can I tell which are mitochondrial?

Check the `trna_id` - mitochondrial tRNAs have "mito-" prefix in the unified output (e.g., `mito-tRNA-Leu-UAA-1-1`).

## Common Issues

### My heatmap has too many gaps

This is normal for the variable loop region (around positions 44-48). Different tRNAs have different variable loop lengths.

### How do I extract the anticodon positions?

The anticodon is always at **Sprinzl positions 34-36** and maps to the same `global_index` for all tRNAs within a coordinate file.

To extract anticodons:
```python
anticodons = df[df['sprinzl_index'].isin([34, 35, 36])]
```

### How do I cite this?

See the Citation section in the [README](README.md#citation) for BibTeX format.

## Still Have Questions?

[Open an issue](https://github.com/lkwhite/tRNAs-in-space/issues) on GitHub!
