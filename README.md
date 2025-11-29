# 🚀🍀 tRNAs in space 🍀🚀

[![Tests](https://github.com/lkwhite/tRNAs-in-space/workflows/Tests/badge.svg)](https://github.com/lkwhite/tRNAs-in-space/actions)
[![Build](https://github.com/lkwhite/tRNAs-in-space/workflows/Build/badge.svg)](https://github.com/lkwhite/tRNAs-in-space/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

*A standardized approach to generate shared tRNA coordinates for plotting.*

## Quick Start

**Using pre-computed coordinates:**
```python
import pandas as pd

# Load E. coli Type I tRNAs with standard labeling (most common)
df = pd.read_csv('outputs/ecoliK12_global_coords_offset0_type1.tsv', sep='\t')

# Create alignment matrix
alignment = df.pivot_table(
    index='trna_id',
    columns='global_index',
    values='residue'
)
```

> **Note:** Coordinate files are grouped by structural characteristics (offset and type) to ensure proper position alignment. See [Coordinate File Organization](#coordinate-file-organization) for details on choosing the right file.

**Installation:**
```bash
# Clone the repository
git clone https://github.com/lkwhite/tRNAs-in-space.git
cd tRNAs-in-space

# Install as a package
pip install -e .

# Or install with visualization tools
pip install -e ".[viz]"
```

**Generate coordinates from your own data:**
```bash
# After running R2DT on your FASTA files
python scripts/trnas_in_space.py ./r2dt_output_dir/ my_output.tsv --split-by-offset-and-type
```

This generates multiple coordinate files grouped by offset and type (e.g., `my_output_offset0_type1.tsv`).

See [examples/01_basic_visualization.ipynb](examples/01_basic_visualization.ipynb) for detailed usage examples.

## Modomics Modification Annotations

Pre-computed modification data from [MODOMICS](https://genesilico.pl/modomics/) is included, mapped to the global coordinate system. This enables comparison of experimental modification detection against known reference modifications from mass spectrometry data.

**Available species:**
- *E. coli*: 261 modification positions across 14 tRNAs (12 modification types)
- *S. cerevisiae*: 131 modification positions across 10 tRNAs (18 modification types)
- *H. sapiens*: 43 modification positions across 16 tRNAs (16 modification types)

**Usage example:**
```python
import pandas as pd

# Load Modomics annotations
mods = pd.read_csv('outputs/modomics/modomics_to_sprinzl_mapping.tsv', sep='\t')

# Join with your global coordinates (using offset0/type1 as example)
coords = pd.read_csv('outputs/ecoliK12_global_coords_offset0_type1.tsv', sep='\t')
annotated = coords.merge(
    mods[['gtRNAdb_trna_id', 'position_gtRNAdb', 'modification_short_name']],
    left_on=['trna_id', 'seq_index'],
    right_on=['gtRNAdb_trna_id', 'position_gtRNAdb'],
    how='left'
)
# Now 'annotated' includes modification_short_name for known modifications
```

See [`docs/archive/modomics-integration/MODOMICS_INTEGRATION.md`](docs/archive/modomics-integration/MODOMICS_INTEGRATION.md) for implementation details and alignment methodology.

## Global Coordinate System

This project provides a **standardized coordinate system for nuclear elongator tRNAs** that enables comparative structural analysis across different tRNA families. The coordinate system supports both Type I (standard) and Type II (extended variable arm) tRNAs, allowing researchers to perform analyses that were not previously possible with individual tRNA studies.

### Research Capabilities

**Cross-tRNA Comparative Analysis:**
- Compare modification patterns across amino acid families
- Analyze structural domain conservation (acceptor stem, anticodon loop, T-arm)
- Study extended variable arm differences between Leu/Ser/Tyr and other tRNAs
- Generate multi-tRNA heatmaps and statistical comparisons

**Position-Specific Studies:**
- Map modification frequencies to standardized structural positions
- Identify hotspots of evolutionary conservation or variation
- Correlate structural features with experimental modification data

**Type I vs Type II Analysis:**
- Compare standard tRNAs (76 nt) with extended variable arm tRNAs (~90 nt)
- Study structural adaptations in Leucine, Serine, and Tyrosine tRNAs
- Analyze how extended arms affect surrounding structural regions

### System Scope

**✅ Supported tRNA Types:**
- **Nuclear elongator tRNAs**: Standard cytoplasmic tRNAs used in protein synthesis
- **Type I**: Alanine, Phenylalanine, Glycine, and most other amino acids (standard structure)
- **Type II**: Leucine, Serine, Tyrosine (extended variable arms with e1-e24 positions)

**🚫 Excluded by Design:**
- **Selenocysteine tRNAs**: Structurally incompatible (~95 nt with unique binding requirements)
- **Mitochondrial tRNAs**: Different architecture (60-75 nt, missing structural features)
- **Initiator methionine tRNAs**: Modified structure for specialized ribosome binding

### Usage Guidelines

**High-Confidence Analyses:**
- Cross-amino-acid modification comparisons
- Structural domain analysis (acceptor, anticodon, T-regions)
- Type I vs Type II extended variable arm studies

**Moderate-Confidence Analyses:**
- Position-specific studies (validation recommended for critical positions)
- Inter-species comparative analysis

**Alternative Approaches Recommended:**
- Fine-grained analysis within single tRNA families (use individual tRNA coordinates)
- Studies requiring selenocysteine or mitochondrial tRNAs (specialized analysis needed)

For detailed analysis guidelines, see [ANALYSIS_GUIDELINES.md](ANALYSIS_GUIDELINES.md). For technical implementation details, see [docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SCOPE.md](docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SCOPE.md).

## Coordinate File Organization

tRNAs are grouped into separate coordinate files based on two structural characteristics to ensure that positions align correctly across all tRNAs within each file.

### Why Grouping?

Different tRNAs have structural variations that affect how Sprinzl positions map to sequence positions. When combined into a single coordinate space, these variations cause **position collisions**—where the same global index maps to different structural positions in different tRNAs. Grouping by offset and type eliminates these collisions.

### Grouping Dimensions

**Offset** (labeling offset from -3 to +1): Reflects variation in how the D-loop region is annotated. Most tRNAs have offset 0 (standard labeling), but some have insertions or deletions that shift the numbering.

**Type** (Type I or Type II):
- **Type I**: Standard tRNAs with short variable loops (most amino acids)
- **Type II**: tRNAs with extended variable arms (Leucine, Serine, Tyrosine)

### Available Coordinate Files

| Organism | Files | Offset Range |
|----------|-------|--------------|
| E. coli K12 | 5 | -1, 0, +1 |
| S. cerevisiae | 4 | -1, 0, +1 |
| H. sapiens | 7 | -3, -1, 0, +1 |

**File naming:** `{species}_global_coords_offset{N}_type{1|2}.tsv`

**Examples:**
- `ecoliK12_global_coords_offset0_type1.tsv` — Standard E. coli Type I tRNAs (most common)
- `hg38_global_coords_offset0_type2.tsv` — Human Leu/Ser/Tyr tRNAs with standard labeling
- `sacCer_global_coords_offset+1_type1.tsv` — Yeast Type I tRNAs with +1 labeling offset

### Choosing a File

For most analyses, start with `offset0_type1` files—these contain the majority of tRNAs with standard structure. If you're studying Leucine, Serine, or Tyrosine tRNAs (which have extended variable arms), use the corresponding `type2` files.

To find which file contains a specific tRNA, check the tRNA's amino acid:
- **Type I** (short variable loop): Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Lys, Met, Phe, Pro, Thr, Trp, Val
- **Type II** (extended variable arm): Leu, Ser, Tyr

---

This README documents how to go from tRNA reference sequences → grouped coordinate files for plotting and cross‑isodecoder comparisons.

Pre-computed coordinate files are available in the `outputs/` directory for *E. coli* K12, *S. cerevisiae*, and *H. sapiens* nuclear elongator tRNAs. Note that *S. cerevisiae* mitochondrial tRNAs may need additional hand-curation for accurate alignment due to [the unavailability of models specific to fungal mitochondria](https://github.com/r2dt-bio/R2DT/issues/197#issuecomment-3201887161). If this is relevant to your work we recommend the alignments in Reinsch and Garcia 2025 (see References).

The documentation and code in this repository can be used to generate coordinate files for your own tRNA sequences.

## The Problem

<img src="https://github.com/user-attachments/assets/f62efd38-594d-437f-9dac-2ef8c92167e1" alt="image" width="443" height="375"/>

tRNA biologists have classically used the **Sprinzl positions** pictured above (named after M. Sprinzl, see *References*) instead of consecutive numbering within each isodecoder\[1\]. This system ensures that homologous structural features line up across different tRNAs. For instance, the anticodon is always assigned to positions 34-36 regardless of whether a particular tRNA sequence is longer or shorter.

This convention is biologically meaningful, but introduces problems for data integration:

1.  **Non-contiguous across isodecoders:** not every tRNA contains every Sprinzl position, so some positions are absent depending on sequence length or loop structure

2.  **Unequal spacing:** gaps in Sprinzl numbering create irregular axes, making it difficult to generate heatmaps or plots that assume equally spaced positions

3.  **Non-integer labels** like 17a, 20a, 20b and the e-notations in the variable loop further complicate use of Sprinzl as a common coordinate system

[R2DT 2.0](https://github.com/r2dt-bio/r2dt) partly addresses this by embedding Sprinzl numbering in its structural templates, allowing researchers to annotate secondary structure images using positional annotations relevant to tRNA biology. This is very useful for RNA structure visualization. But for downstream analysis, a more unified coordinate system is needed.

## A Global Index

For tRNA sequencing (or other positionally anchored assays), it is often more useful to work in a **global coordinate system**:

-   Each nucleotide position is assigned a consecutive integer index (1,2,3...).

-   The index is consistent across all tRNAs, ensuring every position in a heatmap corresponds to an equal-spaced axis.

-   Missing Sprinzl positions can be interpolated or left blank without breaking the regular grid.

Here's an example from [our own work](https://pubmed.ncbi.nlm.nih.gov/39091754/) where we attempted to align nuclear and mitochondrial tRNAs from budding yeast using Sprinzl coordinates. You'll note some positions appear as "missing" (gray), with the large grey region between Sprinzl positions 48 and 49 reflecting variable loop length, where none of the tRNAs displayed contains sequence covering the full set of variable loop positions. <img src="https://github.com/user-attachments/assets/ed9f8001-b2d8-44e2-a044-9c8017b0a89f" alt="Screenshot 2025-08-20 at 5 13 38 PM" width="585" height="220"/>

However, there are still a few issues with the heatmap above.

-   The distribution of "missing" positions in the variable loop don't line up correctly with their Sprinzl annotations, because the Sprinzl annotations along the X axis are actually only added during the plotting step

-   Not all tRNAs are pictured, because these structural alignments were generated from a non-comprehensive `.afa` file

**Note:** This repository now includes properly aligned Modomics modification data mapped to the global coordinate system (see `outputs/modomics/`), which resolves these alignment issues for downstream analysis

By introducing a global index, we eliminate spacing irregularities and enable cross-isodecoder comparison in a clean, standardized coordinate space.

## Implementation

**Goal:** To convert heterogeneous Sprinzl-style labels from [R2DT](https://docs.r2dt.bio/en/latest/index.html) output into a unified coordinate system we need to:

-   Keep per‑nucleotide sequence order (5′→3′).

-   Preserve canonical Sprinzl labels (e.g., 20, 20A).

-   Fill unlabeled residues deterministically with fractional positions.

-   Generate a **global_index** (1..K) so all tRNAs plot on the same x‑axis; missing positions show as NA.

**Inputs** A FASTA file of mature tRNA sequences used for alignment/reference. For consistency, trim adapters out of these sequences if present.

**Outputs**:

-   Multiple TSV files (grouped by offset and type) with per-base fields:

    -   `trna_id`, `seq_index`, `sprinzl_index`, `sprinzl_label`, `residue`

    -   `sprinzl_ordinal`, `sprinzl_continuous`, `global_index`, `region`

## Prerequisites

-   [Docker](https://www.docker.com) (For [R2DT](https://github.com/r2dt-bio/r2dt))

-   Python 3.9+ with `pandas`.

-   Your tRNA reference fasta

-   `trnas_in_space.py` (from this repository): extracts per-nucleotide indices/labels from R2DT-produced `.enriched.json` files, fills Sprinzl gaps, builds global label order, assigns fractional/global indices, and annotates structural regions.

### Step 1: Run R2DT

> *Note:* Make sure Docker Desktop (or another Docker engine) is running before you start. On Mac/Windows you should see the 🐳 whale icon in your menu bar/system tray. You can test with `docker ps` — if it prints a table (even empty), you’re good.

From the project root (with your FASTA files in fastas/), run:

```         
docker run --rm \
  -v "$(pwd):/data" \
  rnacentral/r2dt \
  r2dt.py gtrnadb draw /data/fastas/yourtRNAreference.fa /data/outputs/yourprefix_jsons
```

This runs R2DT in `gtrnadb draw` mode, using covariance models and tRNAscan-SE outputs to annotate tRNAs with structural information. It creates a folder like `outputs/yourprefix_jsons/` containing R2DT `.enriched.json` files that include these fields:

-   `templateResidueIndex` = plain numeric Sprinzl positions

-   `templateNumberingLabel` = full Sprinzl label as a string (numbers + any special suffixes)

### Step 2: Build global coordinates

You can then extract information from the above by running the following script on your R2DT output directory:

```
python trnas_in_space.py ./output ecoliK12_global_coords.tsv --split-by-offset-and-type
```

This generates multiple coordinate files grouped by offset and type (e.g., `ecoliK12_global_coords_offset0_type1.tsv`). The script fills missing `sprinzl_index` values using neighboring positions, and assigns structural regions such as `anticodon-loop`, `acceptor-stem`, etc. Unresolvable cases retain `sprinzl_index` of `-1` and a `region` value of `unknown`.

## Documentation

- **[OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md)** - Detailed specification of output TSV columns
- **[FAQ.md](docs/FAQ.md)** - Frequently asked questions and practical tips
- **[examples/01_basic_visualization.ipynb](examples/01_basic_visualization.ipynb)** - Interactive visualization tutorial
- **[CHANGELOG.md](CHANGELOG.md)** - Version history and release notes
- **[ANALYSIS_GUIDELINES.md](ANALYSIS_GUIDELINES.md)** - Guidelines for using the coordinate system in research
- **[docs/archive/](docs/archive/)** - Historical documentation and completed development notes

## Citation

If you use tRNAs in space in your research, please cite:

**BibTeX:**
```bibtex
@software{trnas_in_space,
  author = {White, Laura K.},
  title = {tRNAs in space: Standardized coordinates for tRNA analysis},
  year = {2025},
  url = {https://github.com/lkwhite/tRNAs-in-space},
  license = {MIT}
}
```

**Related publication:**
```bibtex
@article{white2024comparative,
  author = {White, Laura K. and Dobson, K. and Del Pozo, S. and others},
  title = {Comparative analysis of 43 distinct RNA modifications by nanopore tRNA sequencing},
  journal = {bioRxiv},
  year = {2024},
  doi = {10.1101/2024.07.23.604651}
}
```

## Contributing

Contributions are welcome! Please feel free to:
- Report issues or bugs
- Suggest new features or improvements
- Submit pull requests

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Footnotes

-   Isodecoders: tRNAs that share the same anticodon

-   Isoacceptors: tRNAs charged by the same amino acid

## Implementation History

- **Previous**: Single unified coordinate file per organism (`{species}_global_coords.tsv`). This approach had position collisions when combining tRNAs with different structural characteristics.
- **Current (November 2025)**: Grouped files by offset and type (`{species}_global_coords_offset{N}_type{1|2}.tsv`). Separating tRNAs by structural characteristics eliminates collisions and ensures proper position alignment. Unified files may be supported in future versions with a different approach.

## References

Cappannini A., Ray A., Purta E., Mukherjee S., Boccaletto P., Moafinejad S.N., Lechner A., Barchet C., Klaholz B.P., Stefaniak F., Bujnicki J.M. (2023). MODOMICS: a database of RNA modifications and related information. 2023 update. Nucleic Acids Research 51(D1):D155–D163. <https://doi.org/10.1093/nar/gkad1083> **Resource**: <https://genesilico.pl/modomics/>

Chan P.P., Lowe T.M. (2016). GtRNAdb 2.0: an expanded database of transfer RNA genes identified in complete and draft genomes. Nucleic Acids Research 44(Database issue):D184–D189. <https://doi.org/10.1093/nar/gkv1309> **Resource**: <https://gtrnadb.org>

Chan P.P., Lin B.Y., Mak A.J., Lowe T.M. (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Research 49(16):9077–9096. <https://doi.org/10.1093/nar/gkab688>

McCann H., Meade C.D., Williams L.D., Petrov A.S., Johnson P.Z., Simon A.E., Hoksza D., Nawrocki E.P., Chan P.P., Lowe T.M., Ribas C.E., Sweeney B.A., Madeira F., Anyango S., Appasamy S.D., Deshpande M., Varadi M., Velankar S., Zirbel C.L., Naiden A., Jossinet F., Petrov A.I. (2025). R2DT: a comprehensive platform for visualizing RNA secondary structure. Nucleic Acids Research 53(4):gkaf032. <https://doi.org/10.1093/nar/gkaf032> **Resource**: <https://r2dt.bio>

Reinsch J.L., Garcia D.M. (2025). Concurrent detection of chemically modified bases in yeast mitochondrial tRNAs by nanopore direct RNA sequencing. bioRxiv \[Preprint\]. 2025 May 9:2025.05.09.653160. <https://doi.org/10.1101/2025.05.09.653160>

Sprinzl M., Horn C., Brown M., Ioudovitch A., Steinberg S. (1998). Compilation of tRNA sequences and sequences of tRNA genes. Nucleic Acids Research 26(1):148–153. [https://doi.org/10.1093/nar/26.1.148](url)

White L.K., Dobson K., Del Pozo S., Bilodeaux J.M., Andersen S.E., Baldwin A., Barrington C., Körtel N., Martinez-Seidel F., Strugar S.M., Watt K.E.N., Mukherjee N., Hesselberth J.R. (2024). Comparative analysis of 43 distinct RNA modifications by nanopore tRNA sequencing. bioRxiv \[Preprint\]. 2024 Jul 24:2024.07.23.604651. [https://doi.org/10.1101/2024.07.23.604651](url)
