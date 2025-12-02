#!/usr/bin/env python3
"""
Visualization script for exploring the tRNAs-in-space alignment system.

Generates multiple visualizations showing how the coordinate system aligns
tRNAs with different structural variations.

Usage:
    python scripts/visualize_alignment_system.py

Outputs PNG files to outputs/figures/
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from pathlib import Path


# Color schemes - Okabe-Ito colorblind-friendly palette
# https://jfly.uni-koeln.de/color/
OKABE_ITO = {
    'orange': '#E69F00',
    'skyblue': '#56B4E9',
    'green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'black': '#000000',
}

NUCLEOTIDE_COLORS = {
    'A': OKABE_ITO['green'],      # green
    'C': OKABE_ITO['blue'],       # blue
    'G': OKABE_ITO['orange'],     # orange
    'T': OKABE_ITO['vermillion'], # vermillion/red
    'U': OKABE_ITO['vermillion'], # vermillion/red
}
REGION_COLORS = {
    'acceptor-stem': '#e74c3c',
    'D-stem': '#9b59b6',
    'D-loop': '#8e44ad',
    'anticodon-stem': '#3498db',
    'anticodon-loop': '#2980b9',
    'variable-region': '#27ae60',
    'variable-arm': '#16a085',
    'T-stem': '#f39c12',
    'T-loop': '#d35400',
    'acceptor-tail': '#c0392b',
    'unknown': '#95a5a6',
}


def load_data(filepath: str) -> pd.DataFrame:
    """Load the global coordinates TSV file."""
    df = pd.read_csv(filepath, sep='\t')
    # Convert empty strings to NaN for sprinzl_label
    df['sprinzl_label'] = df['sprinzl_label'].replace('', np.nan)
    return df


def classify_label(label: str) -> str:
    """Classify a Sprinzl label by type."""
    if pd.isna(label) or label == '':
        return 'unlabeled'
    label = str(label)
    if label.startswith('e'):
        return 'e-position'
    if any(c.isalpha() for c in label) and any(c.isdigit() for c in label):
        # Has both letters and numbers like "20a"
        if label[0].isdigit():
            return 'insertion'
    return 'standard'


def viz_01_label_mapping(df: pd.DataFrame, output_dir: Path):
    """
    Visualization 1: Label-to-Coordinate Mapping
    Shows how Sprinzl labels map to global_index coordinates.
    """
    # Get unique label-to-index mappings
    label_map = df[df['sprinzl_label'].notna()].groupby('sprinzl_label').agg({
        'global_index': 'first',
        'region': 'first'
    }).reset_index()

    # Classify each label
    label_map['label_type'] = label_map['sprinzl_label'].apply(classify_label)

    # Sort by global_index
    label_map = label_map.sort_values('global_index')

    # Color mapping for label types
    type_colors = {
        'standard': '#3498db',
        'insertion': '#e74c3c',
        'e-position': '#27ae60',
        'unlabeled': '#95a5a6'
    }

    fig, ax = plt.subplots(figsize=(16, 10))

    # Plot each label
    colors = [type_colors[t] for t in label_map['label_type']]
    ax.scatter(label_map['global_index'], range(len(label_map)),
               c=colors, s=50, alpha=0.7)

    # Add label text
    for i, (_, row) in enumerate(label_map.iterrows()):
        ax.annotate(row['sprinzl_label'],
                   (row['global_index'], i),
                   xytext=(5, 0), textcoords='offset points',
                   fontsize=7, alpha=0.8)

    # Highlight gaps in global_index
    all_indices = set(range(1, int(df['global_index'].max()) + 1))
    used_indices = set(label_map['global_index'].astype(int))
    gaps = sorted(all_indices - used_indices)
    for gap in gaps:
        ax.axvline(gap, color='red', alpha=0.1, linewidth=1)

    ax.set_xlabel('Global Index', fontsize=12)
    ax.set_ylabel('Sprinzl Labels (sorted by global_index)', fontsize=12)
    ax.set_title('S. cerevisiae Type II tRNAs: Sprinzl Label → Global Index Mapping\n'
                 'Red = insertions (20a), Green = e-positions, Blue = standard', fontsize=12)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=t) for t, c in type_colors.items() if t != 'unlabeled']
    ax.legend(handles=legend_patches, loc='upper left')

    ax.set_yticks([])
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / '01_label_to_coord_mapping.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / '01_label_to_coord_mapping.png'}")


def viz_02_heatmap(df: pd.DataFrame, output_dir: Path):
    """
    Visualization 2: Cross-tRNA Alignment Heatmap
    Shows multiple tRNAs aligned by global_index.
    """
    # Get unique tRNAs and select a sample
    trna_ids = df['trna_id'].unique()
    # Take first 25 or all if fewer
    sample_trnas = trna_ids[:25]

    df_sample = df[df['trna_id'].isin(sample_trnas)]

    # Create pivot table
    alignment = df_sample.pivot_table(
        index='trna_id',
        columns='global_index',
        values='residue',
        aggfunc='first'
    )

    # Encode nucleotides as numbers
    nuc_to_num = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
    alignment_numeric = alignment.map(lambda x: nuc_to_num.get(x, 0) if pd.notna(x) else 0)

    # Create figure
    fig, ax = plt.subplots(figsize=(20, 8))

    # Custom colormap: white for gaps, then A/C/G/T colors
    cmap = ListedColormap(['white', '#2ecc71', '#3498db', '#f39c12', '#e74c3c'])

    im = ax.imshow(alignment_numeric.values, aspect='auto', cmap=cmap, vmin=0, vmax=4)

    # Labels
    ax.set_yticks(range(len(alignment.index)))
    ax.set_yticklabels([tid.replace('nuc-tRNA-', '') for tid in alignment.index], fontsize=6)

    # X-axis: show every 10th position
    xticks = range(0, len(alignment.columns), 10)
    ax.set_xticks(xticks)
    ax.set_xticklabels([alignment.columns[i] for i in xticks], fontsize=8)

    ax.set_xlabel('Global Index', fontsize=12)
    ax.set_ylabel('tRNA', fontsize=12)
    ax.set_title('S. cerevisiae Type II tRNAs: Aligned by Global Coordinate\n'
                 'Green=A, Blue=C, Orange=G, Red=T/U, White=gap', fontsize=12)

    # Add region annotations at top
    region_map = df.groupby('global_index')['region'].first()
    prev_region = None
    region_starts = []
    for i, col in enumerate(alignment.columns):
        region = region_map.get(col, 'unknown')
        if region != prev_region:
            region_starts.append((i, region))
            prev_region = region

    # Draw region bars
    for i, (start, region) in enumerate(region_starts):
        end = region_starts[i+1][0] if i < len(region_starts)-1 else len(alignment.columns)
        color = REGION_COLORS.get(region, '#95a5a6')
        ax.axvspan(start-0.5, end-0.5, ymin=1.0, ymax=1.05,
                   color=color, alpha=0.7, clip_on=False)

    plt.tight_layout()
    plt.savefig(output_dir / '02_alignment_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / '02_alignment_heatmap.png'}")


def viz_03_coverage(df: pd.DataFrame, output_dir: Path):
    """
    Visualization 3: Position Coverage
    Shows which positions are present in how many tRNAs.
    """
    n_trnas = df['trna_id'].nunique()

    # Count tRNAs at each position
    coverage = df.groupby('global_index').agg({
        'trna_id': 'nunique',
        'region': 'first'
    }).reset_index()
    coverage.columns = ['global_index', 'count', 'region']

    fig, ax = plt.subplots(figsize=(16, 6))

    # Color bars by region
    colors = [REGION_COLORS.get(r, '#95a5a6') for r in coverage['region']]

    ax.bar(coverage['global_index'], coverage['count'], color=colors, width=1.0, alpha=0.8)

    # Add line at total tRNAs
    ax.axhline(n_trnas, color='black', linestyle='--', alpha=0.5,
               label=f'Total tRNAs: {n_trnas}')

    ax.set_xlabel('Global Index', fontsize=12)
    ax.set_ylabel('Number of tRNAs', fontsize=12)
    ax.set_title('S. cerevisiae Type II tRNAs: Position Coverage\n'
                 'Lower bars = variable positions (insertions, e-region)', fontsize=12)

    # Legend for regions
    unique_regions = coverage['region'].unique()
    legend_patches = [mpatches.Patch(color=REGION_COLORS.get(r, '#95a5a6'), label=r)
                      for r in unique_regions]
    ax.legend(handles=legend_patches, loc='lower right', ncol=3, fontsize=8)

    ax.set_xlim(0, coverage['global_index'].max() + 1)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / '03_position_coverage.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / '03_position_coverage.png'}")


def viz_04_ruler_tracks(df: pd.DataFrame, output_dir: Path):
    """
    Visualization 4: Ruler/Tracks View
    Shows horizontal tracks comparing several tRNAs.
    """
    # Select 4 different tRNAs (try to get Leu, Ser, Tyr)
    trna_ids = df['trna_id'].unique()
    sample_trnas = []
    for aa in ['Leu', 'Ser', 'Tyr']:
        matches = [t for t in trna_ids if aa in t]
        if matches:
            sample_trnas.append(matches[0])
    # Add one more if needed
    if len(sample_trnas) < 4:
        for t in trna_ids:
            if t not in sample_trnas:
                sample_trnas.append(t)
                break

    sample_trnas = sample_trnas[:4]

    fig, axes = plt.subplots(len(sample_trnas), 1, figsize=(18, 3*len(sample_trnas)),
                              sharex=True)
    if len(sample_trnas) == 1:
        axes = [axes]

    max_global = int(df['global_index'].max())

    for ax, trna_id in zip(axes, sample_trnas):
        trna_df = df[df['trna_id'] == trna_id].sort_values('global_index')

        # Plot each position as a colored box
        for _, row in trna_df.iterrows():
            gi = row['global_index']
            residue = row['residue']
            label = row['sprinzl_label'] if pd.notna(row['sprinzl_label']) else '?'
            region = row['region']

            # Box color by nucleotide
            color = NUCLEOTIDE_COLORS.get(residue, '#95a5a6')

            # Draw box
            rect = mpatches.Rectangle((gi-0.4, 0.2), 0.8, 0.6,
                                       facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)

            # Add residue letter
            ax.text(gi, 0.5, residue, ha='center', va='center', fontsize=7, fontweight='bold')

            # Add Sprinzl label above
            label_type = classify_label(label)
            label_color = 'red' if label_type == 'insertion' else 'green' if label_type == 'e-position' else 'black'
            ax.text(gi, 1.0, str(label), ha='center', va='bottom', fontsize=5,
                   rotation=90, color=label_color)

        # Highlight gaps
        present_indices = set(trna_df['global_index'].astype(int))
        for gi in range(1, max_global + 1):
            if gi not in present_indices:
                ax.axvspan(gi-0.4, gi+0.4, color='lightgray', alpha=0.3)

        ax.set_xlim(0, max_global + 1)
        ax.set_ylim(0, 1.8)
        ax.set_yticks([])
        short_name = trna_id.replace('nuc-tRNA-', '')
        ax.set_ylabel(short_name, fontsize=9, rotation=0, ha='right', va='center')
        ax.set_title(f'{short_name} ({len(trna_df)} positions)', fontsize=10, loc='left')

    axes[-1].set_xlabel('Global Index', fontsize=12)
    fig.suptitle('S. cerevisiae Type II tRNAs: Track Comparison\n'
                 'Labels above: black=standard, red=insertion, green=e-position',
                 fontsize=12, y=1.02)

    plt.tight_layout()
    plt.savefig(output_dir / '04_ruler_tracks.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / '04_ruler_tracks.png'}")


def viz_05_arrow_schematic(df: pd.DataFrame, output_dir: Path):
    """
    Visualization 5: Arrow Schematic
    Shows coordinate transformation with arrows for one tRNA.
    """
    # Pick one tRNA with interesting features
    trna_id = df['trna_id'].iloc[0]
    trna_df = df[df['trna_id'] == trna_id].sort_values('seq_index')

    fig, ax = plt.subplots(figsize=(20, 8))

    # Two rows: top = seq_index, bottom = global_index
    top_y = 2
    bottom_y = 0

    # Scale factor for x positions
    scale = 1.0

    # Draw top row (sequence positions)
    for _, row in trna_df.iterrows():
        si = row['seq_index'] * scale
        gi = row['global_index'] * scale
        residue = row['residue']
        label = row['sprinzl_label'] if pd.notna(row['sprinzl_label']) else '?'

        color = NUCLEOTIDE_COLORS.get(residue, '#95a5a6')

        # Top box (seq_index)
        rect_top = mpatches.FancyBboxPatch((si-0.3, top_y), 0.6, 0.6,
                                            boxstyle="round,pad=0.02",
                                            facecolor=color, edgecolor='black')
        ax.add_patch(rect_top)
        ax.text(si, top_y+0.3, residue, ha='center', va='center', fontsize=8, fontweight='bold')
        ax.text(si, top_y+0.8, str(int(row['seq_index'])), ha='center', va='bottom', fontsize=6)

        # Bottom box (global_index)
        rect_bot = mpatches.FancyBboxPatch((gi-0.3, bottom_y), 0.6, 0.6,
                                            boxstyle="round,pad=0.02",
                                            facecolor=color, edgecolor='black')
        ax.add_patch(rect_bot)
        ax.text(gi, bottom_y+0.3, residue, ha='center', va='center', fontsize=8, fontweight='bold')
        ax.text(gi, bottom_y-0.2, str(label), ha='center', va='top', fontsize=6,
               color='red' if classify_label(label) in ['insertion', 'e-position'] else 'black')

        # Arrow connecting them
        arrow_color = 'red' if si != gi else 'gray'
        alpha = 0.8 if si != gi else 0.3
        ax.annotate('', xy=(gi, bottom_y+0.6), xytext=(si, top_y),
                   arrowprops=dict(arrowstyle='->', color=arrow_color, alpha=alpha, lw=0.5))

    # Labels
    ax.text(-2, top_y+0.3, 'seq_index\n(5\'→3\')', ha='right', va='center', fontsize=10)
    ax.text(-2, bottom_y+0.3, 'global_index\n(aligned)', ha='right', va='center', fontsize=10)

    max_x = max(trna_df['seq_index'].max(), trna_df['global_index'].max()) * scale
    ax.set_xlim(-5, max_x + 2)
    ax.set_ylim(-1, 3.5)
    ax.set_aspect('equal')
    ax.axis('off')

    short_name = trna_id.replace('nuc-tRNA-', '')
    ax.set_title(f'Coordinate Transformation: {short_name}\n'
                 'Top: sequence position → Bottom: aligned global_index\n'
                 'Red arrows = position shifts, Gray arrows = no change', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_dir / '05_arrow_schematic.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / '05_arrow_schematic.png'}")


def viz_06_text_alignment(df: pd.DataFrame, output_dir: Path,
                          suffix: str = None, title: str = None,
                          max_trnas: int = 5):
    """
    Visualization 6: Full Coordinate Path Alignment View
    Shows all 4 columns: global_index, sprinzl_label (shared), seq_index, residue (per-tRNA).

    Parameters:
        df: DataFrame with tRNA alignment data
        output_dir: Directory to save output
        suffix: Optional suffix for filename (e.g., 'offset0_type2')
        title: Optional custom title
        max_trnas: Maximum number of tRNAs to show (default 5)
    """
    # Select a sample of tRNAs - auto-detect what amino acids are present
    trna_ids = df['trna_id'].unique()

    # Extract amino acid from tRNA IDs and find representative samples
    aa_to_trnas = {}
    for trna_id in trna_ids:
        # Extract amino acid from tRNA ID (e.g., 'nuc-tRNA-Leu-CAA-1-1' -> 'Leu')
        parts = trna_id.replace('nuc-tRNA-', '').split('-')
        if parts:
            aa = parts[0]
            if aa not in aa_to_trnas:
                aa_to_trnas[aa] = []
            aa_to_trnas[aa].append(trna_id)

    # Select diverse samples: prioritize variety across amino acids
    sample_trnas = []
    for aa in sorted(aa_to_trnas.keys()):
        if len(sample_trnas) >= max_trnas:
            break
        # Take first tRNA of each amino acid type
        sample_trnas.append(aa_to_trnas[aa][0])

    # If we have fewer amino acids than max, add more from same aa
    if len(sample_trnas) < max_trnas:
        for aa in sorted(aa_to_trnas.keys()):
            for trna in aa_to_trnas[aa][1:]:
                if len(sample_trnas) >= max_trnas:
                    break
                if trna not in sample_trnas:
                    sample_trnas.append(trna)
            if len(sample_trnas) >= max_trnas:
                break

    sample_trnas = sample_trnas[:max_trnas]

    df_sample = df[df['trna_id'].isin(sample_trnas)]

    # Create alignment matrices for residue and seq_index
    alignment_residue = df_sample.pivot_table(
        index='trna_id',
        columns='global_index',
        values='residue',
        aggfunc='first'
    ).fillna('-')

    alignment_seqidx = df_sample.pivot_table(
        index='trna_id',
        columns='global_index',
        values='seq_index',
        aggfunc='first'
    )  # Keep NaN for gaps

    # Get labels for header
    label_map = df.groupby('global_index')['sprinzl_label'].first()

    # Build the alignment as strings for cleaner rendering
    cols = alignment_residue.columns.tolist()
    n_cols = len(cols)
    n_trnas = len(alignment_residue.index)

    # Figure size based on content - each tRNA needs 2 rows (seq_index + residue)
    fig_width = max(18, n_cols * 0.20)
    fig_height = max(6, (n_trnas * 1.8 + 4) * 0.5 + 3)  # condensed rows per tRNA + headers + legend
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Use a tight grid layout
    row_height = 1.0
    col_width = 0.8
    name_width = 14  # characters for tRNA name

    # Starting positions
    x_start = name_width * col_width + 1
    y_start = (n_trnas * 1.8 + 3) * row_height + 4  # offset for legend space (condensed)

    # --- SHARED HEADER ROWS ---
    y = y_start

    # Draw global_index header (angled 45 degrees) - FIRST (top)
    ax.text(0, y, 'global_index:', fontsize=10, fontfamily='monospace', fontweight='bold', va='center')
    for i, col in enumerate(cols):
        x = x_start + i * col_width
        ax.text(x, y, str(col), fontsize=8, fontfamily='monospace', ha='left', va='bottom',
                color=OKABE_ITO['purple'], rotation=45)

    # Draw sprinzl_label header (angled 45 degrees) - SECOND (tighter spacing)
    y -= row_height * 0.6
    ax.text(0, y, 'sprinzl_label:', fontsize=10, fontfamily='monospace', fontweight='bold', va='center')
    for i, col in enumerate(cols):
        x = x_start + i * col_width
        raw_label = label_map.get(col, None)
        # Convert None/NaN to '-', and format numeric labels as integers
        if raw_label is None or pd.isna(raw_label) or str(raw_label) in ('None', 'nan', ''):
            label = '-'
        elif isinstance(raw_label, float) and raw_label == int(raw_label):
            # Numeric label stored as float - convert to int string
            label = str(int(raw_label))
        else:
            label = str(raw_label)
        # Color by label type (Okabe-Ito)
        label_type = classify_label(label)
        if label_type == 'insertion':
            color = OKABE_ITO['vermillion']
        elif label_type == 'e-position':
            color = OKABE_ITO['skyblue']
        else:
            color = OKABE_ITO['black']
        ax.text(x, y, label, fontsize=9, fontfamily='monospace', ha='left', va='bottom',
                color=color, fontweight='bold', rotation=45)

    # Separator line
    y -= row_height * 0.5
    ax.axhline(y, color='black', linewidth=1.5, xmin=0, xmax=1)

    # --- PER-tRNA ROWS ---
    for trna_id in alignment_residue.index:
        # tRNA name header
        y -= row_height * 0.6
        short_name = trna_id.replace('nuc-tRNA-', '')
        ax.text(0, y, short_name, fontsize=10, fontfamily='monospace', va='center', fontweight='bold')

        # seq_index row (tighter)
        y -= row_height * 0.45
        ax.text(0, y, '  seq_index:', fontsize=8, fontfamily='monospace', va='center', color='gray')
        seq_idx_row = alignment_seqidx.loc[trna_id]
        for i, col in enumerate(cols):
            x = x_start + i * col_width
            val = seq_idx_row[col]
            if pd.isna(val):
                ax.text(x, y, '-', fontsize=8, fontfamily='monospace', ha='center', va='center', color='lightgray')
            else:
                ax.text(x, y, str(int(val)), fontsize=8, fontfamily='monospace', ha='center', va='center', color='gray')

        # residue row (tighter)
        y -= row_height * 0.45
        ax.text(0, y, '  residue:', fontsize=8, fontfamily='monospace', va='center', color='gray')
        residue_row = alignment_residue.loc[trna_id]
        for i, col in enumerate(cols):
            x = x_start + i * col_width
            residue = residue_row[col]
            if residue == '-':
                ax.text(x, y, '·', fontsize=10, fontfamily='monospace', ha='center', va='center', color='lightgray')
            else:
                color = NUCLEOTIDE_COLORS.get(residue, 'black')
                ax.text(x, y, residue, fontsize=10, fontfamily='monospace', ha='center', va='center',
                        color=color, fontweight='bold')

        # Small gap between tRNAs
        y -= row_height * 0.2

    # Store bottom of alignment for legend positioning
    alignment_bottom = y - row_height * 0.5

    # --- LEGEND SECTION ---
    legend_y = alignment_bottom - row_height * 1.0
    ax.axhline(legend_y + row_height * 0.5, color='gray', linewidth=0.5, linestyle='--', xmin=0, xmax=1)

    # Legend title
    legend_y -= row_height * 0.3
    ax.text(0, legend_y, 'LEGEND', fontsize=11, fontweight='bold', va='top')

    # Two columns for legend
    col1_x = 0
    col2_x = x_start + n_cols * col_width * 0.35

    # Column 1: Shared alignment coordinates
    legend_y -= row_height * 0.7
    ax.text(col1_x, legend_y, 'Shared (alignment coordinates)', fontsize=10, fontweight='bold', va='top')

    legend_y -= row_height * 0.5
    ax.text(col1_x, legend_y, 'global_index:', fontsize=9, fontweight='bold', va='top', style='italic',
            color=OKABE_ITO['purple'])
    ax.text(col1_x + 8, legend_y, 'unified coordinate (computed)', fontsize=9, va='top')

    legend_y -= row_height * 0.4
    ax.text(col1_x, legend_y, 'sprinzl_label:', fontsize=9, fontweight='bold', va='top', style='italic')
    ax.text(col1_x + 8, legend_y, 'structural position (from R2DT)', fontsize=9, va='top')

    legend_y -= row_height * 0.5
    # sprinzl_label color key
    for label_type, color, lbl, xoff in [
        ('standard', OKABE_ITO['black'], 'Standard (1-76)', 0),
        ('insertion', OKABE_ITO['vermillion'], 'Insertion (20a)', 12),
        ('e-position', OKABE_ITO['skyblue'], 'e-position (e1-e24)', 24)
    ]:
        rect = mpatches.Rectangle((col1_x + xoff, legend_y - 0.12), 0.25, 0.25, facecolor=color)
        ax.add_patch(rect)
        ax.text(col1_x + xoff + 0.4, legend_y, lbl, fontsize=8, va='center')

    # Column 2: Per-tRNA data
    legend_y2 = alignment_bottom - row_height * 1.7
    ax.text(col2_x, legend_y2, 'Per-tRNA', fontsize=10, fontweight='bold', va='top')

    legend_y2 -= row_height * 0.5
    ax.text(col2_x, legend_y2, 'seq_index:', fontsize=9, fontweight='bold', va='top', style='italic', color='gray')
    ax.text(col2_x + 6, legend_y2, "position in sequence 5'→3' (1-based)", fontsize=9, va='top')

    legend_y2 -= row_height * 0.4
    ax.text(col2_x, legend_y2, 'residue:', fontsize=9, fontweight='bold', va='top', style='italic')
    ax.text(col2_x + 6, legend_y2, 'nucleotide (from input FASTA)', fontsize=9, va='top')

    legend_y2 -= row_height * 0.5
    # Nucleotide colors (Okabe-Ito)
    for nuc, color, offset in [('A', NUCLEOTIDE_COLORS['A'], 0),
                                ('C', NUCLEOTIDE_COLORS['C'], 2),
                                ('G', NUCLEOTIDE_COLORS['G'], 4),
                                ('T/U', NUCLEOTIDE_COLORS['T'], 6)]:
        rect = mpatches.Rectangle((col2_x + offset, legend_y2 - 0.12), 0.25, 0.25, facecolor=color)
        ax.add_patch(rect)
        ax.text(col2_x + offset + 0.35, legend_y2, nuc, fontsize=8, va='center', fontweight='bold')

    legend_y2 -= row_height * 0.4
    ax.text(col2_x, legend_y2, '· or - = gap (position not present in this tRNA)', fontsize=8, va='center', color='gray')

    # Set axis limits
    final_bottom = min(legend_y, legend_y2) - row_height * 1.0
    ax.set_xlim(-1, x_start + n_cols * col_width + 1)
    ax.set_ylim(final_bottom, y_start + row_height * 1.5)
    ax.axis('off')

    # Dynamic title based on suffix or custom title
    if title:
        fig_title = title
    elif suffix:
        # Parse suffix for readable title
        fig_title = f'S. cerevisiae tRNAs ({suffix}): Full Coordinate Alignment Path'
    else:
        fig_title = 'S. cerevisiae Type II tRNAs: Full Coordinate Alignment Path'

    ax.set_title(fig_title, fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()

    # Dynamic filename based on suffix
    if suffix:
        output_filename = f'sacCer_{suffix}_alignment.png'
    else:
        output_filename = '06_text_alignment.png'

    plt.savefig(output_dir / output_filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / output_filename}")


def generate_all_yeast_alignments(output_dir: Path, outputs_dir: Path):
    """
    Generate alignment visualizations for all yeast offset groups.

    Parameters:
        output_dir: Directory to save figures
        outputs_dir: Directory containing the global_coords TSV files
    """
    # All yeast data files with their suffixes and descriptive titles
    YEAST_FILES = [
        ('sacCer_global_coords_offset-1_type2.tsv', 'offset-1_type2',
         'S. cerevisiae Type II tRNAs (offset -1): Full Coordinate Alignment Path'),
        ('sacCer_global_coords_offset+1_type1.tsv', 'offset+1_type1',
         'S. cerevisiae Type I tRNAs (offset +1): Full Coordinate Alignment Path'),
        ('sacCer_global_coords_offset0_type1.tsv', 'offset0_type1',
         'S. cerevisiae Type I tRNAs (offset 0): Full Coordinate Alignment Path'),
        ('sacCer_global_coords_offset0_type2.tsv', 'offset0_type2',
         'S. cerevisiae Type II tRNAs (offset 0): Full Coordinate Alignment Path'),
    ]

    print("Generating alignment visualizations for all yeast offset groups...")
    print()

    for filename, suffix, title in YEAST_FILES:
        filepath = outputs_dir / filename
        if not filepath.exists():
            print(f"  Skipping {filename} (file not found)")
            continue

        print(f"Processing: {filename}")
        df = load_data(filepath)
        n_trnas = df['trna_id'].nunique()
        print(f"  Loaded {len(df)} rows, {n_trnas} tRNAs")

        # Determine max_trnas based on group size
        max_trnas = min(5, n_trnas)

        viz_06_text_alignment(df, output_dir, suffix=suffix, title=title, max_trnas=max_trnas)
        print()

    print("All yeast alignment visualizations complete.")


def main():
    """Generate all visualizations."""
    import argparse

    parser = argparse.ArgumentParser(description='Generate tRNA alignment visualizations')
    parser.add_argument('--yeast-all', action='store_true',
                        help='Generate alignment visualizations for all yeast offset groups')
    args = parser.parse_args()

    # Paths
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
    outputs_dir = project_dir / 'outputs'
    output_dir = outputs_dir / 'figures'

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.yeast_all:
        # Generate only the alignment visualizations for all yeast offset groups
        generate_all_yeast_alignments(output_dir, outputs_dir)
    else:
        # Default: generate all visualizations for offset0_type2
        input_file = outputs_dir / 'sacCer_global_coords_offset0_type2.tsv'

        print(f"Loading data from: {input_file}")
        df = load_data(input_file)
        print(f"Loaded {len(df)} rows, {df['trna_id'].nunique()} tRNAs")
        print()

        # Generate all visualizations
        print("Generating visualizations...")
        viz_01_label_mapping(df, output_dir)
        viz_02_heatmap(df, output_dir)
        viz_03_coverage(df, output_dir)
        viz_04_ruler_tracks(df, output_dir)
        viz_05_arrow_schematic(df, output_dir)
        viz_06_text_alignment(df, output_dir)

        print()
        print(f"All visualizations saved to: {output_dir}")


if __name__ == '__main__':
    main()
