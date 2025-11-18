#!/usr/bin/env python3
"""
validate_coordinate_system.py

Comprehensive validation suite for tRNA global coordinate system.
Tests all success criteria from the assessment plan:
1. Shared coordinate system for Type I and Type II tRNAs
2. Retention of sprinzl_label values
3. Region boundary alignment with Sprinzl positions
4. Consistent global_index mapping for positions 1, 73, 76

Usage:
    python scripts/validate_coordinate_system.py outputs/ecoliK12_global_coords.tsv
    python scripts/validate_coordinate_system.py outputs/*.tsv  # validate all
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Any


class CoordinateValidator:
    """Validator for tRNA global coordinate system."""

    def __init__(self, filepath: str):
        """Initialize validator with coordinate file."""
        self.filepath = Path(filepath)
        self.df = pd.read_csv(filepath, sep='\t')
        self.species_name = self.filepath.stem
        self.results = {}

    def run_all_tests(self) -> Dict[str, Any]:
        """Run all validation tests and return results."""
        print(f"\n{'='*80}")
        print(f"VALIDATION REPORT: {self.species_name}")
        print(f"{'='*80}\n")

        self.results['file'] = str(self.filepath)
        self.results['n_trnas'] = self.df['trna_id'].nunique()
        self.results['n_rows'] = len(self.df)
        self.results['global_index_range'] = (
            int(self.df['global_index'].min()),
            int(self.df['global_index'].max())
        )

        print(f"File: {self.filepath.name}")
        print(f"Total tRNAs: {self.results['n_trnas']}")
        print(f"Total positions: {self.results['n_rows']}")
        print(f"Global index range: {self.results['global_index_range'][0]} - {self.results['global_index_range'][1]}")

        # Run all tests
        self.test_excluded_types()
        self.test_type_sharing()
        self.test_sprinzl_retention()
        self.test_region_boundaries()
        self.test_critical_positions()
        self.test_e_position_ordering()
        self.test_collisions()

        # Overall score
        self.calculate_quality_score()

        return self.results

    def test_excluded_types(self):
        """Test 0: Check for excluded tRNA types that should be filtered."""
        print(f"\n{'-'*80}")
        print("TEST 0: Excluded tRNA Types")
        print(f"{'-'*80}")

        excluded_types = {
            'SeC': self.df['trna_id'].str.contains('SeC|Sec', case=False, na=False).sum(),
            'Mitochondrial': self.df['trna_id'].str.contains('mito', case=False, na=False).sum(),
            'iMet': self.df['trna_id'].str.contains('iMet', case=False, na=False).sum(),
            'fMet': self.df['trna_id'].str.contains('fMet', case=False, na=False).sum(),
        }

        total_excluded = sum(excluded_types.values())

        if total_excluded == 0:
            print("✅ PASS: No excluded tRNA types found")
            self.results['excluded_types'] = {'status': 'PASS', 'count': 0}
        else:
            print(f"❌ FAIL: Found {total_excluded} excluded tRNA types:")
            for type_name, count in excluded_types.items():
                if count > 0:
                    examples = self.df[self.df['trna_id'].str.contains(
                        type_name, case=False, na=False
                    )]['trna_id'].unique()[:3]
                    print(f"  - {type_name}: {count} tRNAs (e.g., {', '.join(examples)})")
            self.results['excluded_types'] = {
                'status': 'FAIL',
                'count': total_excluded,
                'details': excluded_types
            }

    def test_type_sharing(self):
        """Test 1: Verify Type I and Type II can share global coordinate system."""
        print(f"\n{'-'*80}")
        print("TEST 1: Type I & Type II Coordinate Sharing")
        print(f"{'-'*80}")

        # Identify Type II amino acids (Leu, Ser, Tyr with extended variable arms)
        self.df['amino_acid'] = self.df['trna_id'].str.extract(r'tRNA-([A-Za-z]+)-')[0]

        type2_aas = ['Leu', 'Ser', 'Tyr']
        has_e_positions = self.df['sprinzl_label'].str.match(r'^e\d+$', na=False).any()

        type2_trnas = self.df[self.df['amino_acid'].isin(type2_aas)]['trna_id'].unique()
        type1_trnas = self.df[~self.df['amino_acid'].isin(type2_aas)]['trna_id'].unique()

        print(f"Type I tRNAs (standard): {len(type1_trnas)}")
        print(f"Type II tRNAs (Leu/Ser/Tyr): {len(type2_trnas)}")
        print(f"Extended arm positions (e-values) present: {has_e_positions}")

        if has_e_positions:
            e_positions = self.df[self.df['sprinzl_label'].str.match(r'^e\d+$', na=False)]
            e_global_range = (e_positions['global_index'].min(), e_positions['global_index'].max())
            print(f"  e-positions occupy global_index: {e_global_range[0]} - {e_global_range[1]}")

            # Check if e-positions are in reserved space (should be roughly 65-85)
            if e_global_range[0] >= 60 and e_global_range[1] <= 90:
                print("  ✅ e-positions in expected reserved space")
                status = 'PASS'
            else:
                print("  ⚠️  e-positions outside expected reserved space (60-90)")
                status = 'PARTIAL'
        else:
            print("  ℹ️  No Type II extended arm positions in dataset")
            status = 'PASS'

        self.results['type_sharing'] = {
            'status': status,
            'type1_count': len(type1_trnas),
            'type2_count': len(type2_trnas),
            'has_e_positions': has_e_positions
        }

    def test_sprinzl_retention(self):
        """Test 2: Verify sprinzl_label values are retained."""
        print(f"\n{'-'*80}")
        print("TEST 2: Sprinzl Label Retention")
        print(f"{'-'*80}")

        total_rows = len(self.df)
        non_null_labels = self.df['sprinzl_label'].notna().sum()
        non_empty_labels = (self.df['sprinzl_label'].astype(str) != '').sum()

        retention_rate = non_empty_labels / total_rows * 100

        print(f"Total rows: {total_rows}")
        print(f"Non-null sprinzl_label: {non_null_labels} ({non_null_labels/total_rows*100:.1f}%)")
        print(f"Non-empty sprinzl_label: {non_empty_labels} ({retention_rate:.1f}%)")

        if retention_rate >= 95:
            print(f"✅ PASS: {retention_rate:.1f}% retention (≥95% threshold)")
            status = 'PASS'
        elif retention_rate >= 90:
            print(f"⚠️  PARTIAL: {retention_rate:.1f}% retention (90-95%)")
            status = 'PARTIAL'
        else:
            print(f"❌ FAIL: {retention_rate:.1f}% retention (<90%)")
            status = 'FAIL'

        self.results['sprinzl_retention'] = {
            'status': status,
            'retention_rate': retention_rate
        }

    def test_region_boundaries(self):
        """Test 3: Verify regions align with expected Sprinzl position ranges."""
        print(f"\n{'-'*80}")
        print("TEST 3: Region Boundary Alignment")
        print(f"{'-'*80}")

        # Expected region assignments based on Sprinzl positions
        expected_regions = {
            '1-7': 'acceptor-stem',
            '8-13': 'D-stem',
            '14-21': 'D-loop',
            '22-25': 'D-stem',
            '26-31': 'anticodon-stem',
            '32-38': 'anticodon-loop',
            '39-43': 'anticodon-stem',
            '44-48': 'variable-region',
            '49-53': 'T-stem',
            '54-60': 'T-loop',
            '61-65': 'T-stem',
            '66-72': 'acceptor-stem',
            '73-76': 'acceptor-tail'
        }

        # Check critical positions
        test_positions = {
            1: 'acceptor-stem',
            8: 'D-stem',  # Was falling through to unknown
            9: 'D-stem',  # Was falling through to unknown
            26: 'anticodon-stem',  # Was falling through to unknown
            34: 'anticodon-loop',
            45: 'variable-region',
            54: 'T-loop',
            73: 'acceptor-tail'
        }

        mismatches = []
        for pos, expected_region in test_positions.items():
            pos_str = str(pos)
            pos_df = self.df[self.df['sprinzl_label'] == pos_str]
            if len(pos_df) > 0:
                actual_regions = pos_df['region'].unique()
                if expected_region not in actual_regions:
                    mismatches.append((pos, expected_region, actual_regions))
                    print(f"  ❌ Position {pos}: expected '{expected_region}', got {actual_regions}")
                else:
                    print(f"  ✅ Position {pos}: '{expected_region}' ✓")

        # Check for unknown regions in positions 1-76
        unknown_positions = self.df[
            (self.df['region'] == 'unknown') &
            (self.df['sprinzl_index'].between(1, 76))
        ]
        unknown_count = len(unknown_positions)

        if unknown_count > 0:
            unknown_labels = unknown_positions['sprinzl_label'].unique()[:10]
            print(f"\n  ⚠️  {unknown_count} positions (1-76) marked as 'unknown': {list(unknown_labels)}")

        # Check e-positions get "variable-arm"
        e_positions = self.df[self.df['sprinzl_label'].str.match(r'^e\d+$', na=False)]
        if len(e_positions) > 0:
            e_regions = e_positions['region'].unique()
            if len(e_regions) == 1 and e_regions[0] == 'variable-arm':
                print(f"\n  ✅ All e-positions assigned to 'variable-arm'")
            else:
                print(f"\n  ❌ e-positions have incorrect regions: {e_regions}")
                mismatches.append(('e*', 'variable-arm', e_regions))

        if len(mismatches) == 0 and unknown_count == 0:
            print(f"\n✅ PASS: All region boundaries correctly aligned")
            status = 'PASS'
        elif len(mismatches) <= 2:
            print(f"\n⚠️  PARTIAL: {len(mismatches)} mismatches, {unknown_count} unknown")
            status = 'PARTIAL'
        else:
            print(f"\n❌ FAIL: {len(mismatches)} mismatches, {unknown_count} unknown")
            status = 'FAIL'

        self.results['region_boundaries'] = {
            'status': status,
            'mismatches': len(mismatches),
            'unknown_count': unknown_count
        }

    def test_critical_positions(self):
        """Test 4: Verify positions 1, 73, 76 map consistently to global_index."""
        print(f"\n{'-'*80}")
        print("TEST 4: Critical Position Consistency (1, 73, 76)")
        print(f"{'-'*80}")

        critical_positions = ['1', '73', '76']
        position_results = {}

        for pos in critical_positions:
            pos_df = self.df[self.df['sprinzl_label'] == pos]

            if len(pos_df) == 0:
                print(f"\nSprinzl {pos}: NOT FOUND")
                position_results[pos] = {'status': 'MISSING', 'unique_indices': 0}
                continue

            unique_indices = pos_df['global_index'].unique()
            n_unique = len(unique_indices)
            n_trnas = len(pos_df['trna_id'].unique())

            print(f"\nSprinzl position {pos}:")
            print(f"  Appears in {len(pos_df)} rows across {n_trnas} tRNAs")
            print(f"  Maps to {n_unique} unique global_index value(s): {sorted(unique_indices)}")

            if n_unique == 1:
                print(f"  ✅ PASS: Maps to single global_index ({unique_indices[0]})")
                position_results[pos] = {'status': 'PASS', 'unique_indices': 1}
            else:
                print(f"  ❌ FAIL: Maps to {n_unique} different global_index values")
                # Show which tRNAs cause variation
                for gidx in sorted(unique_indices):
                    trnas = pos_df[pos_df['global_index'] == gidx]['trna_id'].unique()
                    print(f"    global_index {gidx}: {len(trnas)} tRNAs (e.g., {trnas[0]})")
                position_results[pos] = {
                    'status': 'FAIL',
                    'unique_indices': n_unique,
                    'index_values': list(unique_indices)
                }

        # Overall status
        statuses = [r['status'] for r in position_results.values()]
        if all(s == 'PASS' for s in statuses):
            overall = 'PASS'
        elif position_results.get('1', {}).get('status') == 'PASS':
            overall = 'PARTIAL'  # Position 1 is critical
        else:
            overall = 'FAIL'

        self.results['critical_positions'] = {
            'status': overall,
            'positions': position_results
        }

    def test_e_position_ordering(self):
        """Test 5: Verify e-positions sort in correct order."""
        print(f"\n{'-'*80}")
        print("TEST 5: Extended Variable Arm (e-position) Ordering")
        print(f"{'-'*80}")

        e_positions = self.df[self.df['sprinzl_label'].str.match(r'^e\d+$', na=False)]

        if len(e_positions) == 0:
            print("ℹ️  No e-positions in dataset (Type I tRNAs only)")
            self.results['e_position_ordering'] = {'status': 'N/A', 'count': 0}
            return

        # Get unique e-labels and their global_index
        e_mapping = e_positions.groupby('sprinzl_label')['global_index'].first().sort_index()

        print(f"Found {len(e_mapping)} unique e-positions:")
        for label, gidx in e_mapping.items():
            print(f"  {label}: global_index {gidx}")

        # Check if global_index values are monotonically increasing
        indices = list(e_mapping.values)
        is_ordered = all(indices[i] < indices[i+1] for i in range(len(indices)-1))

        if is_ordered:
            print(f"\n✅ PASS: e-positions have monotonically increasing global_index")
            status = 'PASS'
        else:
            print(f"\n❌ FAIL: e-positions not properly ordered by global_index")
            status = 'FAIL'

        self.results['e_position_ordering'] = {
            'status': status,
            'count': len(e_mapping)
        }

    def test_collisions(self):
        """Test 6: Check for global_index collisions."""
        print(f"\n{'-'*80}")
        print("TEST 6: Global Index Collision Detection")
        print(f"{'-'*80}")

        # Find global_index values with multiple different sprinzl_labels
        collision_check = self.df.groupby('global_index')['sprinzl_label'].apply(
            lambda x: [lbl for lbl in x.unique() if pd.notna(lbl) and lbl != '']
        )
        collisions = collision_check[collision_check.apply(len) > 1]

        if len(collisions) == 0:
            print("✅ PASS: No global_index collisions detected")
            self.results['collisions'] = {'status': 'PASS', 'count': 0}
        else:
            print(f"❌ FAIL: Found {len(collisions)} global_index collision points:")
            # Show first 5 examples
            for i, (gidx, labels) in enumerate(collisions.items()):
                if i >= 5:
                    print(f"  ... and {len(collisions) - 5} more collisions")
                    break
                print(f"\n  global_index {gidx} has {len(labels)} labels: {labels}")
                for lbl in labels[:2]:  # Show 2 examples per collision
                    trnas = self.df[
                        (self.df['global_index'] == gidx) &
                        (self.df['sprinzl_label'] == lbl)
                    ]['trna_id'].unique()
                    if len(trnas) > 0:
                        print(f"    '{lbl}': {len(trnas)} tRNAs (e.g., {trnas[0]})")

            self.results['collisions'] = {
                'status': 'FAIL',
                'count': len(collisions)
            }

    def calculate_quality_score(self):
        """Calculate overall quality score (0-100)."""
        print(f"\n{'='*80}")
        print("OVERALL QUALITY ASSESSMENT")
        print(f"{'='*80}\n")

        scores = {
            'excluded_types': 15,      # Critical
            'type_sharing': 10,
            'sprinzl_retention': 15,   # Critical
            'region_boundaries': 20,   # Critical
            'critical_positions': 25,  # Most critical
            'e_position_ordering': 5,
            'collisions': 10
        }

        total_score = 0
        for test_name, max_points in scores.items():
            if test_name not in self.results:
                continue

            status = self.results[test_name].get('status', 'N/A')

            if status == 'PASS':
                points = max_points
            elif status == 'PARTIAL':
                points = max_points * 0.5
            elif status == 'N/A':
                points = max_points  # Don't penalize for N/A
            else:
                points = 0

            total_score += points
            print(f"{test_name:25s}: {points:4.1f}/{max_points} ({status})")

        print(f"\n{'-'*80}")
        print(f"TOTAL QUALITY SCORE: {total_score:.1f}/100")
        print(f"{'-'*80}")

        if total_score >= 90:
            quality = "EXCELLENT - Ready for production use"
        elif total_score >= 75:
            quality = "GOOD - Minor issues to address"
        elif total_score >= 50:
            quality = "FAIR - Significant issues need fixing"
        else:
            quality = "POOR - Major remediation required"

        print(f"\nQuality Level: {quality}\n")

        self.results['quality_score'] = {
            'total': total_score,
            'max': 100,
            'level': quality
        }


def main():
    parser = argparse.ArgumentParser(
        description="Validate tRNA global coordinate system files"
    )
    parser.add_argument(
        'files',
        nargs='+',
        help='Coordinate TSV file(s) to validate'
    )
    parser.add_argument(
        '--output',
        help='Output validation report to file (markdown format)',
        default=None
    )

    args = parser.parse_args()

    all_results = []

    for filepath in args.files:
        validator = CoordinateValidator(filepath)
        results = validator.run_all_tests()
        all_results.append(results)

    # Summary if multiple files
    if len(all_results) > 1:
        print(f"\n{'='*80}")
        print(f"SUMMARY: {len(all_results)} files validated")
        print(f"{'='*80}\n")
        for result in all_results:
            filename = Path(result['file']).name
            score = result['quality_score']['total']
            print(f"{filename:40s}: {score:.1f}/100")

    return 0 if all(r['quality_score']['total'] >= 75 for r in all_results) else 1


if __name__ == "__main__":
    sys.exit(main())
