#!/usr/bin/env python3
"""
Add global coordinate columns to MODOMICS mapping files.

Enriches MODOMICS TSV files with:
- offset: labeling offset (-1, 0, +1, etc.)
- structural_type: type1 or type2
- global_index: position in the unified coordinate space

Joins on gtRNAdb_trna_id + sprinzl_label to look up the global_index
from the appropriate offset×type coordinate file.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Tuple


def parse_coord_filename(filename: str) -> Optional[str]:
    """
    Parse organism from coordinate filename.

    Examples:
        sacCer_global_coords.tsv -> 'sacCer'
        sacCer_mito_global_coords.tsv -> 'sacCer'
        hg38_global_coords.tsv -> 'hg38'
    """
    # Match both nuclear and mito coordinate files
    match = re.match(r'(.+?)(?:_mito)?_global_coords\.tsv', filename)
    if match:
        return match.group(1)
    return None


def load_coordinate_files(output_dir: Path) -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Load unified coordinate files and build lookup tables.

    Returns:
        Dict mapping organism_id -> {trna_id -> {sprinzl_label: global_index}}
    """
    # Structure: organism -> trna_id -> {sprinzl_label: global_index}
    coord_lookup: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(
        lambda: defaultdict(dict)
    )

    # Load both nuclear and mito coord files
    for coord_file in output_dir.glob("*_global_coords.tsv"):
        organism = parse_coord_filename(coord_file.name)
        if not organism:
            continue

        with open(coord_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                trna_id = row['trna_id']
                sprinzl_label = row.get('sprinzl_label', '')
                global_index_str = row.get('global_index', '')

                # Skip rows with missing sprinzl_label or global_index
                if not sprinzl_label or not global_index_str:
                    continue

                global_index = int(global_index_str)

                coord_lookup[organism][trna_id][sprinzl_label] = global_index

    return coord_lookup


# Map MODOMICS species names to organism IDs used in coordinate files
SPECIES_TO_ORGANISM = {
    'Escherichia coli': 'ecoliK12',
    'Saccharomyces cerevisiae': 'sacCer',
    'Homo sapiens': 'hg38',
}


def enrich_modomics_file(
    input_path: Path,
    output_path: Path,
    coord_lookup: Dict,
):
    """
    Add global coordinate columns to a MODOMICS mapping file.
    """
    with open(input_path, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        input_fieldnames = list(reader.fieldnames)

        # Remove existing coordinate columns if present (for idempotent re-runs)
        # Note: offset and structural_type are no longer used in unified format
        old_cols = ['offset', 'structural_type', 'global_index']
        base_fieldnames = [f for f in input_fieldnames if f not in old_cols]

        # Add global_index column (offset/type no longer needed)
        output_fieldnames = base_fieldnames + ['global_index']

        rows = []
        stats = {'matched': 0, 'unmatched': 0, 'unknown_species': 0}

        for row in reader:
            species = row['species']
            trna_id = row['gtRNAdb_trna_id']
            sprinzl_label = row['sprinzl_label']

            # Clean up old columns if present
            row.pop('offset', None)
            row.pop('structural_type', None)

            # Look up organism
            organism = SPECIES_TO_ORGANISM.get(species)
            if not organism:
                row['global_index'] = ''
                stats['unknown_species'] += 1
                rows.append(row)
                continue

            # Look up global_index
            if (organism in coord_lookup and
                trna_id in coord_lookup[organism] and
                sprinzl_label in coord_lookup[organism][trna_id]):

                global_index = coord_lookup[organism][trna_id][sprinzl_label]
                row['global_index'] = str(global_index)
                stats['matched'] += 1
            else:
                row['global_index'] = ''
                stats['unmatched'] += 1

            rows.append(row)

    # Write output
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Add global coordinate columns to MODOMICS mapping files'
    )
    parser.add_argument(
        '--modomics-dir',
        default='outputs/modomics',
        help='Directory containing MODOMICS TSV files'
    )
    parser.add_argument(
        '--coords-dir',
        default='outputs',
        help='Directory containing offset×type coordinate files'
    )
    parser.add_argument(
        '--in-place',
        action='store_true',
        help='Modify files in place (otherwise creates .enriched.tsv files)'
    )

    args = parser.parse_args()

    modomics_dir = Path(args.modomics_dir)
    coords_dir = Path(args.coords_dir)

    # Load coordinate files
    print(f"Loading coordinate files from {coords_dir}...")
    coord_lookup = load_coordinate_files(coords_dir)

    organisms_loaded = list(coord_lookup.keys())
    print(f"  Loaded coordinates for: {', '.join(organisms_loaded)}")

    # Process MODOMICS files - both species-specific and combined
    modomics_files = list(modomics_dir.glob("*_modomics_to_sprinzl.tsv"))

    # Also include the combined file
    combined = modomics_dir / "modomics_to_sprinzl_mapping.tsv"
    if combined.exists():
        modomics_files.append(combined)

    if not modomics_files:
        print(f"No MODOMICS files found in {modomics_dir}")
        return 1

    print(f"\nProcessing {len(modomics_files)} MODOMICS file(s)...")

    for input_path in modomics_files:
        if args.in_place:
            output_path = input_path
        else:
            # Insert .enriched before .tsv
            output_path = input_path.with_suffix('.enriched.tsv')
            if input_path.suffix == '.tsv':
                output_path = input_path.parent / (input_path.stem + '.with_coords.tsv')

        print(f"\n  {input_path.name}:")
        stats = enrich_modomics_file(input_path, output_path, coord_lookup)

        print(f"    Matched: {stats['matched']}")
        print(f"    Unmatched: {stats['unmatched']}")
        if stats['unknown_species'] > 0:
            print(f"    Unknown species: {stats['unknown_species']}")
        print(f"    Output: {output_path.name}")

    print("\nDone!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
