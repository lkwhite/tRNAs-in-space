# tRNAs-in-space
_A standardized approach to generate shared tRNA coordinate space for plotting._

This README documents how to go from adapter‑trimmed tRNA reference sequences → R2DT drawings → a single shared, equal‑spaced coordinate axis for plotting and cross‑isodecoder comparisons.

## The Problem
tRNA biologists have classifically used **Sprinzl positions** (named after M. Sprinzl, see _References_) instead of consecutive numbering within each  isodecoder[1]. This system ensures that homologous structural features line up across different tRNA, for instance, always assigning the anticodon to positions 34-36 regardless of whether a particular tRNA sequence is longer or shorter.

This convention is biologically meaningful, but introduces problems for data integration: 
* **Non-contiguous across isodecoders:** not every tRNA contains every Sprinzl position, so some positions are absent depending on sequence length or loop structure
* **Unequal spacing:** gaps in Sprinzl numbering create irregular axes, making it difficult to generate heatmaps or plots that assume equally spaced positions

R2DT 2.0 partly addresses this by embedding Sprinzl numbering in its structural templates, allowing researchers to annotate secondary structure images using positional annotations relavant to tRNA biology. But for downstream analysis, a more unified coordinate system is needed.

## A Global Index
For tRNA sequencing (or other positionally anchored assays), it is often more useful to work in a **global coordinate system**:
* Each nucleotide position is assigned a consecutive integer index (1,2,3...).
* The index is consistent across all tRNAs, ensuring every position in a heatmap corresponds to an equal-spaced axis.
* Missing Sprinzl positions can be interpolated or left blank without breaking the regular grid.

Here's an example from [our own work](https://pubmed.ncbi.nlm.nih.gov/39091754/) where we attempted to align nuclear and mitochondrial tRNAs from budding yeast using Sprinzl coordinates. You'll note some positions appear as "missing" (gray), with the large grey box between Sprinzl positions 48 and 49 reflecting variable loop length, where none of the tRNAs displayed contains sequence covering the full set of variable loop positions. 
<img width="585" height="220" alt="Screenshot 2025-08-20 at 5 13 38 PM" src="https://github.com/user-attachments/assets/ed9f8001-b2d8-44e2-a044-9c8017b0a89f" />

However, there are still a few issues with the heatmap above.
* The distribution of "missing" positions in the variable loop don't line up correctly with their Sprinzl annotations, because
* the Sprinzl annotations along the X axis are actually only added during the plotting step, and 
* not all tRNAs are pictured, because these structural alignments were generated from a non-comprehensive `.afa` file downloaded from [Modomics]([url](https://genesilico.pl/modomics/rnafamilies/rf00005/))

By introducing a global index, we eliminate spacing irregularities and enable cross-isodecoder comparison in a clean, standardized coordinate space.



## Footnotes
[1] Isodecoders: tRNAs that share the same anticodon; isoacceptors: tRNAs charged by the same amino acid.

## References
Sprinzl M., Horn C., Brown M., Ioudovitch A., Steinberg S. (1998). Compilation of tRNA sequences and sequences of tRNA genes. Nucleic Acids Research 26(1):148–153. [https://doi.org/10.1093/nar/26.1.148](url)

White L.K., Dobson K., Del Pozo S., Bilodeaux J.M., Andersen S.E., Baldwin A., Barrington C., Körtel N., Martinez-Seidel F., Strugar S.M., Watt K.E.N., Mukherjee N., Hesselberth J.R. (2024). Comparative analysis of 43 distinct RNA modifications by nanopore tRNA sequencing. bioRxiv [Preprint]. 2024 Jul 24:2024.07.23.604651. [https://doi.org/10.1101/2024.07.23.604651](url)

