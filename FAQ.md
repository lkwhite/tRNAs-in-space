# Frequently Asked Questions

## Using Pre-computed Coordinates

### Which organisms are available?

Currently: E. coli K12, S. cerevisiae, and H. sapiens. See `outputs/` directory.

### How do I load the data?

**Python:**
```python
import pandas as pd
df = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')
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

## Generating New Coordinates

### How do I process my own tRNAs?

1. Run R2DT on your FASTA file (requires Docker):
   ```bash
   docker run --rm -v "$(pwd):/data" rnacentral/r2dt \
     r2dt.py gtrnadb draw /data/fastas/my_trnas.fa /data/outputs/my_jsons
   ```

2. Generate coordinates:
   ```bash
   python scripts/trnas_in_space.py outputs/my_jsons outputs/my_coords.tsv
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

The standardized position number (1, 2, 3...) that's consistent across all tRNAs. Use this as your X-axis for alignments and heatmaps.

### Why are some `global_index` values missing for certain tRNAs?

tRNAs vary in length, especially in the variable loop. Missing values mean that tRNA doesn't have a nucleotide at that position.

### What's the difference between `sprinzl_index` and `global_index`?

- `sprinzl_index`: Traditional Sprinzl numbering (1-76, with gaps and letter suffixes)
- `global_index`: Continuous numbering (1, 2, 3...) for easy plotting

## Mitochondrial tRNAs

### Are mitochondrial tRNAs accurate?

For yeast mitochondrial tRNAs: **use with caution**. R2DT lacks fungal mitochondrial-specific models, so alignments may be inaccurate. See README note about Reinsch & Garcia 2025.

For human mitochondrial tRNAs: Generally good quality.

### How can I tell which are mitochondrial?

Check the `trna_id` - mitochondrial tRNAs typically have "mito" or "MT-" in their identifier.

## Common Issues

### My heatmap has too many gaps

This is normal for the variable loop region (around positions 44-48). Different tRNAs have different variable loop lengths.

### The anticodon isn't at the same `global_index` for all tRNAs

Correct! The anticodon is always at **Sprinzl positions 34-36**, but may map to different `global_index` values depending on insertions/deletions earlier in the sequence.

To extract anticodons:
```python
anticodons = df[df['sprinzl_index'].isin([34, 35, 36])]
```

### How do I cite this?

See the Citation section in the [README](README.md#citation) for BibTeX format.

## Still Have Questions?

[Open an issue](https://github.com/lkwhite/tRNAs-in-space/issues) on GitHub!
