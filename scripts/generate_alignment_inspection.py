#!/usr/bin/env python3
"""
Generate a transformed TSV view for manual inspection of tRNA alignments.
Pivots the data to show global_index as rows and tRNA IDs as columns.
"""

import pandas as pd
import sys

def generate_alignment_inspection(input_file, output_file):
    """
    Transform the global coordinates TSV for alignment inspection.

    Args:
        input_file: Path to the input TSV file
        output_file: Path to the output TSV file
    """
    # Read the input TSV
    df = pd.read_csv(input_file, sep='\t')

    # Create a pivot table
    # Rows: global_index, Columns: trna_id, Values: sprinzl_label
    pivot = df.pivot_table(
        index='global_index',
        columns='trna_id',
        values='sprinzl_label',
        aggfunc='first'  # In case of duplicates, take the first
    )

    # Sort columns alphabetically
    pivot = pivot[sorted(pivot.columns)]

    # Fill missing values with 'NA'
    pivot = pivot.fillna('NA')

    # Ensure all global_index values from min to max are present
    min_idx = int(pivot.index.min())
    max_idx = int(pivot.index.max())
    full_index = range(min_idx, max_idx + 1)
    pivot = pivot.reindex(full_index, fill_value='NA')

    # Reset index to make global_index a column
    pivot = pivot.reset_index()
    pivot = pivot.rename(columns={'index': 'global_index'})

    # Save to TSV
    pivot.to_csv(output_file, sep='\t', index=False)

    print(f"Generated alignment inspection file: {output_file}")
    print(f"Dimensions: {len(pivot)} rows (global_index) × {len(pivot.columns)-1} columns (tRNA IDs)")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python generate_alignment_inspection.py <input_tsv> <output_tsv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    generate_alignment_inspection(input_file, output_file)
