#!/usr/bin/env python3
"""
investigate_3prime_variation.py

Investigate structural variation at the 3' end of tRNAs (positions 73-76).
This helps understand why positions 73 and 76 map to multiple global_index values.

Usage:
    python scripts/investigate_3prime_variation.py outputs/ecoliK12_global_coords.tsv
"""

import sys
import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict


def analyze_3prime_structure(filepath: str):
    """Analyze 3' end structural variation across tRNAs."""
    df = pd.read_csv(filepath, sep='\t')
    species = Path(filepath).stem

    print(f"\n{'='*80}")
    print(f"3' END STRUCTURAL VARIATION ANALYSIS: {species}")
    print(f"{'='*80}\n")

    print(f"Total tRNAs: {df['trna_id'].nunique()}")
    print(f"Total positions: {len(df)}\n")

    # Analyze each tRNA's 3' structure
    trna_3prime_structures = {}

    for trna_id in df['trna_id'].unique():
        trna_df = df[df['trna_id'] == trna_id].sort_values('seq_index')

        # Get all positions from 72 onwards
        tail_positions = trna_df[trna_df['sprinzl_index'] >= 72]

        # What positions does this tRNA have?
        has_positions = {}
        for pos in [72, 73, 74, 75, 76]:
            pos_exists = (tail_positions['sprinzl_label'] == str(pos)).any()
            has_positions[pos] = pos_exists

        # Check for insertions (positions with letters like 73a, 74a)
        insertions = tail_positions[
            tail_positions['sprinzl_label'].str.contains('[a-z]', case=False, na=False, regex=True)
        ]['sprinzl_label'].unique()

        # Determine structure pattern
        last_pos = max([p for p, exists in has_positions.items() if exists], default=72)

        structure_key = f"ends_at_{last_pos}"
        if len(insertions) > 0:
            structure_key += f"_insertions_{'-'.join(sorted(insertions))}"

        trna_3prime_structures[trna_id] = {
            'has_73': has_positions[73],
            'has_74': has_positions[74],
            'has_75': has_positions[75],
            'has_76': has_positions[76],
            'last_position': last_pos,
            'insertions': list(insertions),
            'structure_pattern': structure_key,
            'length': len(trna_df)
        }

    # Group tRNAs by structure pattern
    structure_patterns = defaultdict(list)
    for trna_id, info in trna_3prime_structures.items():
        structure_patterns[info['structure_pattern']].append(trna_id)

    print(f"{'='*80}")
    print("3' STRUCTURE PATTERNS")
    print(f"{'='*80}\n")

    print(f"Found {len(structure_patterns)} distinct 3' structure patterns:\n")

    for pattern, trnas in sorted(structure_patterns.items(), key=lambda x: -len(x[1])):
        print(f"\n{pattern}:")
        print(f"  Count: {len(trnas)} tRNAs")
        print(f"  Examples: {', '.join(trnas[:5])}")

        # Show global_index mapping for positions 73, 76 in this pattern
        example_trna = trnas[0]
        example_df = df[df['trna_id'] == example_trna]

        for pos in [73, 76]:
            pos_row = example_df[example_df['sprinzl_label'] == str(pos)]
            if len(pos_row) > 0:
                gidx = pos_row['global_index'].values[0]
                print(f"    Position {pos} -> global_index {gidx}")

    # Analyze positions 73 and 76 specifically
    print(f"\n{'='*80}")
    print("POSITION-SPECIFIC ANALYSIS")
    print(f"{'='*80}\n")

    for pos in [73, 76]:
        pos_df = df[df['sprinzl_label'] == str(pos)]

        if len(pos_df) == 0:
            print(f"Position {pos}: NOT FOUND in any tRNA\n")
            continue

        print(f"Position {pos}:")
        print(f"  Present in: {len(pos_df['trna_id'].unique())} tRNAs")
        print(f"  Absent from: {df['trna_id'].nunique() - len(pos_df['trna_id'].unique())} tRNAs\n")

        # Global index mapping
        gidx_mapping = pos_df.groupby('global_index')['trna_id'].apply(list).to_dict()
        print(f"  Maps to {len(gidx_mapping)} different global_index values:")

        for gidx, trnas in sorted(gidx_mapping.items()):
            print(f"    global_index {gidx}: {len(trnas)} tRNAs")

            # What structure pattern do these tRNAs have?
            patterns_in_group = defaultdict(int)
            for trna in trnas:
                pattern = trna_3prime_structures[trna]['structure_pattern']
                patterns_in_group[pattern] += 1

            for pattern, count in sorted(patterns_in_group.items(), key=lambda x: -x[1]):
                print(f"      {pattern}: {count}")

        print()

    # Check correlation with tRNA length
    print(f"{'='*80}")
    print("CORRELATION WITH tRNA LENGTH")
    print(f"{'='*80}\n")

    length_to_trnas = defaultdict(list)
    for trna_id, info in trna_3prime_structures.items():
        length_to_trnas[info['length']].append(trna_id)

    print(f"tRNA length distribution:")
    for length in sorted(length_to_trnas.keys()):
        trnas = length_to_trnas[length]
        print(f"  {length} nt: {len(trnas)} tRNAs")

        # Show global_index for position 73 in this length group
        if len(trnas) > 0:
            example = trnas[0]
            example_df = df[df['trna_id'] == example]
            pos73 = example_df[example_df['sprinzl_label'] == '73']
            if len(pos73) > 0:
                gidx = pos73['global_index'].values[0]
                print(f"    Position 73 -> global_index {gidx}")

    # Recommendations
    print(f"\n{'='*80}")
    print("RECOMMENDATIONS")
    print(f"{'='*80}\n")

    if len(structure_patterns) == 1:
        print("✅ All tRNAs have identical 3' structure - variation is NOT structural")
        print("   Issue is likely algorithmic (sprinzl_continuous rounding)")
    elif len(structure_patterns) <= 3:
        print("⚠️  Small number of structure patterns - could normalize")
        print("   Consider padding shorter tRNAs to match longest pattern")
    else:
        print("❌ High diversity in 3' structure patterns")
        print("   Positions 73/76 variation reflects biological reality")
        print("   Options:")
        print("   1. Document variation as expected (biological)")
        print("   2. Create '3prime_structure_class' metadata column")
        print("   3. Allow users to filter by structure class if needed")


def main():
    parser = argparse.ArgumentParser(
        description="Investigate 3' end structural variation in tRNA coordinates"
    )
    parser.add_argument(
        'file',
        help='Coordinate TSV file to analyze'
    )

    args = parser.parse_args()

    analyze_3prime_structure(args.file)


if __name__ == "__main__":
    main()
