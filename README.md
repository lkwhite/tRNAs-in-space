# 🚀🍀 tRNAs in space 🍀🚀

*A standardized approach to generate shared tRNA coordinate space for plotting.*

This README documents how to go from adapter‑trimmed tRNA reference sequences → R2DT drawings → a single shared, equal‑spaced coordinate axis for plotting and cross‑isodecoder comparisons.

In this repository, you can also find pre-computed tables containing these indexed coordinates (`outputs/*continuous.tsv`) for tRNAs from common model organisms as we build them. The documentation below and scripts in this repository can be further used to build such tables from scratch.

## The Problem

<img width="443" height="375" alt="image" src="https://github.com/user-attachments/assets/f62efd38-594d-437f-9dac-2ef8c92167e1" />

tRNA biologists have classically used **Sprinzl positions** (named after M. Sprinzl, see *References*) instead of consecutive numbering within each isodecoder\[1\]. This system ensures that homologous structural features line up across different tRNAs. For instance, the anticodon is always assigned to positions 34-36 regardless of whether a particular tRNA sequence is longer or shorter.

This convention is biologically meaningful, but introduces problems for data integration:

1.  **Non-contiguous across isodecoders:** not every tRNA contains every Sprinzl position, so some positions are absent depending on sequence length or loop structure

2.  **Unequal spacing:** gaps in Sprinzl numbering create irregular axes, making it difficult to generate heatmaps or plots that assume equally spaced positions

[R2DT 2.0](https://github.com/r2dt-bio/r2dt) partly addresses this by embedding Sprinzl numbering in its structural templates, allowing researchers to annotate secondary structure images using positional annotations relavant to tRNA biology. But for downstream analysis, a more unified coordinate system is needed.

## A Global Index

For tRNA sequencing (or other positionally anchored assays), it is often more useful to work in a **global coordinate system**:

-   Each nucleotide position is assigned a consecutive integer index (1,2,3...).

-   The index is consistent across all tRNAs, ensuring every position in a heatmap corresponds to an equal-spaced axis.

-   Missing Sprinzl positions can be interpolated or left blank without breaking the regular grid.

Here's an example from [our own work](https://pubmed.ncbi.nlm.nih.gov/39091754/) where we attempted to align nuclear and mitochondrial tRNAs from budding yeast using Sprinzl coordinates. You'll note some positions appear as "missing" (gray), with the large grey region between Sprinzl positions 48 and 49 reflecting variable loop length, where none of the tRNAs displayed contains sequence covering the full set of variable loop positions. <img src="https://github.com/user-attachments/assets/ed9f8001-b2d8-44e2-a044-9c8017b0a89f" alt="Screenshot 2025-08-20 at 5 13 38 PM" width="585" height="220"/>

However, there are still a few issues with the heatmap above.

-   The distribution of "missing" positions in the variable loop don't line up correctly with their Sprinzl annotations, because

-   the Sprinzl annotations along the X axis are actually only added during the plotting step, and

-   not all tRNAs are pictured, because these structural alignments were generated from a non-comprehensive `.afa` file downloaded from [Modomics](%5Burl%5D(https://genesilico.pl/modomics/rnafamilies/rf00005/))

By introducing a global index, we eliminate spacing irregularities and enable cross-isodecoder comparison in a clean, standardized coordinate space.

## Implementation

**Goal:** To convert heterogeneous Sprinzl-style labels from [R2DT](https://docs.r2dt.bio/en/latest/index.html) output into a unified coordinate system we need to:

-   Keep per‑nucleotide sequence order (5′→3′).

-   Preserve canonical Sprinzl labels (e.g., 20, 20A).

-   Fill unlabeled residues deterministically with fractional positions.

-   Generate a **global_index** (1..K) so all tRNAs plot on the same x‑axis; missing positions show as NA.

**Inputs** A FASTA file of mature tRNA sequences used for alignment/reference. For consistency, trim adapters out of these sequences if present.

**Outputs**:

-   Per‑tRNA JSON drawings from R2DT (e.g., `.enriched.json`).

-   A combined TSV with per‑base fields: `trna_id`, `seq_index`, `sprinzl_index`, `sprinzl_label`, `residue`.

-   A second TSV that adds: `sprinzl_ordinal`, `sprinzl_continuous`, `global_index`.

## Prerequisites

-   [Docker](https://www.docker.com) (For [R2DT](https://github.com/r2dt-bio/r2dt))

-   Python 3.9+ with `pandas`.

-   Your tRNA reference fasta

-   Files from this repository:

    -   `r2dt_collect_sprinzl.py` – collects per‑nucleotide indices/labels from R2DT-produced `.enriched.json` files into one TSV.

    -   `make_sprinzl_continuous.py` – builds global label order for tRNAs in the tRNA input, fills fractional positions, and emits global_index.

### Step 1: Run R2DT to generate drawings

From the directory containing your FASTA reference, run:

```         
docker run --rm \
-v "$(pwd):/data" \
rnacentral/r2dt \
r2dt.py gtrnadb draw /data/yourtRNAreference.fa /data/output
```

This runs R2DT in `gtrnadb draw` mode, using covariance models and tRNAscan-SE outputs to annotate tRNAs with structural information, generating `.json` files and secondary structure images for each tRNA in the reference.

### Step 2: Extract per-base indices from R2DT JSON → combined TSV

We pull information from all `*.enriched.json` files in `/output` using the script below and writes them to a shared `tsv` file. Example:

`python r2dt_collect_sprinzl.py ./output ecoli_tRNAs_sprinzl.tsv`

Note that this script overwrites the values of `sprinzl_index` with a best-effort fill fill based on neighboring positions. Unresolvable cases remain `-1`.

### Step 3: Build shared coordinates and equal-spaced global index

This step derives a cohort-wide order of labels for all tRNAs in the file (e.g., `… 20 < 20A < 20B < 21 …` and assigns fractional positions between neighboring labels within each tRNA. It then converts all unique fractional coordinates in the file to a `global_index` for plotting.

`python make_sprinzl_continuous.py ecoli_tRNAs_sprinzl.tsv ecoli_tRNAs_continuous.tsv`

Adds columns:

-   `sprinzl_ordinal` – integer bin from label order (for tick annotations).

-   `sprinzl_continuous` – per‑tRNA fractional coordinate (integers at labeled bins; fractions fill gaps).

-   `global_index` – equal‑spaced 1..K from sorted unique `sprinzl_continuous` across all tRNAs (ideal x‑axis for heatmaps/line plots).

*Example workflow in full:*

```         
# 1) Run R2DT on adapter‑trimmed tRNAs
mkdir -p output
docker run --rm -v "$(pwd):/data" rnacentral/r2dt \
r2dt.py gtrnadb draw /data/ecoliK12MG1655-mature-tRNAs.fa /data/output

# 2) Extract per‑base labels into one TSV
python r2dt_collect_sprinzl.py ./output ecoli_tRNAs_sprinzl.tsv

# 3) Build continuous positions and shared global index
python make_sprinzl_continuous.py ecoli_tRNAs_sprinzl.tsv ecoli_tRNAs_continuous.tsv
```

## Footnotes

\[1\] Isodecoders: tRNAs that share the same anticodon; isoacceptors: tRNAs charged by the same amino acid.

## References

Cappannini A., Ray A., Purta E., Mukherjee S., Boccaletto P., Moafinejad S.N., Lechner A., Barchet C., Klaholz B.P., Stefaniak F., Bujnicki J.M. (2023). MODOMICS: a database of RNA modifications and related information. 2023 update. Nucleic Acids Research 51(D1):D155–D163. <https://doi.org/10.1093/nar/gkad1083> **Resource**: <https://genesilico.pl/modomics/>

Chan P.P., Lowe T.M. (2016). GtRNAdb 2.0: an expanded database of transfer RNA genes identified in complete and draft genomes. Nucleic Acids Research 44(Database issue):D184–D189. <https://doi.org/10.1093/nar/gkv1309> **Resource**: <https://gtrnadb.org>

Chan P.P., Lin B.Y., Mak A.J., Lowe T.M. (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Research 49(16):9077–9096. <https://doi.org/10.1093/nar/gkab688>

McCann H., Meade C.D., Williams L.D., Petrov A.S., Johnson P.Z., Simon A.E., Hoksza D., Nawrocki E.P., Chan P.P., Lowe T.M., Ribas C.E., Sweeney B.A., Madeira F., Anyango S., Appasamy S.D., Deshpande M., Varadi M., Velankar S., Zirbel C.L., Naiden A., Jossinet F., Petrov A.I. (2025). R2DT: a comprehensive platform for visualizing RNA secondary structure. Nucleic Acids Research 53(4):gkaf032. <https://doi.org/10.1093/nar/gkaf032> **Resource**: <https://r2dt.bio>

Sprinzl M., Horn C., Brown M., Ioudovitch A., Steinberg S. (1998). Compilation of tRNA sequences and sequences of tRNA genes. Nucleic Acids Research 26(1):148–153. [https://doi.org/10.1093/nar/26.1.148](url)

White L.K., Dobson K., Del Pozo S., Bilodeaux J.M., Andersen S.E., Baldwin A., Barrington C., Körtel N., Martinez-Seidel F., Strugar S.M., Watt K.E.N., Mukherjee N., Hesselberth J.R. (2024). Comparative analysis of 43 distinct RNA modifications by nanopore tRNA sequencing. bioRxiv \[Preprint\]. 2024 Jul 24:2024.07.23.604651. [https://doi.org/10.1101/2024.07.23.604651](url)
