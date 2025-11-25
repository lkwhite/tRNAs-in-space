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
    # Rows: trna_id, Columns: global_index, Values: sprinzl_label
    pivot = df.pivot_table(
        index='trna_id',
        columns='global_index',
        values='sprinzl_label',
        aggfunc='first'  # In case of duplicates, take the first
    )

    # Ensure all global_index values from min to max are present
    min_idx = int(pivot.columns.min())
    max_idx = int(pivot.columns.max())
    full_columns = range(min_idx, max_idx + 1)
    pivot = pivot.reindex(columns=full_columns, fill_value='NA')

    # Fill missing values with 'NA'
    pivot = pivot.fillna('NA')

    # Sort rows (tRNA IDs) alphabetically
    pivot = pivot.sort_index()

    # Reset index to make trna_id a column
    pivot = pivot.reset_index()
    pivot = pivot.rename(columns={'trna_id': 'tRNA'})

    # Save to TSV
    pivot.to_csv(output_file, sep='\t', index=False)

    print(f"Generated alignment inspection file: {output_file}")
    print(f"Dimensions: {len(pivot)} rows (tRNAs) × {len(pivot.columns)-1} columns (global_index positions)")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python generate_alignment_inspection.py <input_tsv> <output_tsv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    generate_alignment_inspection(input_file, output_file)
