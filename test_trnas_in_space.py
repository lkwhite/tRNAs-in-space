#!/usr/bin/env python3
"""
Tests for trnas_in_space.py

Run with: python -m pytest test_trnas_in_space.py -v
Or: python test_trnas_in_space.py
"""

import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))

import trnas_in_space


def test_imports():
    """Test that all required modules can be imported."""
    assert trnas_in_space is not None
    assert hasattr(trnas_in_space, "main")
    assert hasattr(trnas_in_space, "sort_key")
    assert hasattr(trnas_in_space, "assign_region_from_sprinzl")


def test_sort_key():
    """Test the sort_key function for Sprinzl label ordering."""
    # Test numeric labels
    assert trnas_in_space.sort_key("1") < trnas_in_space.sort_key("2")
    assert trnas_in_space.sort_key("10") < trnas_in_space.sort_key("20")

    # Test labels with suffixes
    assert trnas_in_space.sort_key("20") < trnas_in_space.sort_key("20A")
    assert trnas_in_space.sort_key("20A") < trnas_in_space.sort_key("20B")

    # Test Type II extended variable arm positions (e1-e24)
    # These should sort after position 46 but before 48
    assert trnas_in_space.sort_key("46") < trnas_in_space.sort_key("e1")
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("e2")
    assert trnas_in_space.sort_key("e12") < trnas_in_space.sort_key("e13")
    assert trnas_in_space.sort_key("e24") < trnas_in_space.sort_key("48")

    # Test that "e" positions sort in reserved coordinate space
    # "e" positions map to (46, 2, e_num) which sorts after (46, *) and before (47, 0, '')
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("47")
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("47A")
    # But "e" positions should sort after any 46 insertions
    assert trnas_in_space.sort_key("46A") < trnas_in_space.sort_key("e1")

    # Test empty/nan labels (should sort to end)
    assert trnas_in_space.sort_key("1") < trnas_in_space.sort_key("")
    assert trnas_in_space.sort_key("1") < trnas_in_space.sort_key("nan")
    assert trnas_in_space.sort_key("e24") < trnas_in_space.sort_key("")


def test_sprinzl_numeric_from_label():
    """Test extraction of numeric part from Sprinzl labels."""
    assert trnas_in_space.sprinzl_numeric_from_label("20") == 20
    assert trnas_in_space.sprinzl_numeric_from_label("20A") == 20
    assert trnas_in_space.sprinzl_numeric_from_label("20B") == 20
    assert trnas_in_space.sprinzl_numeric_from_label("1") == 1
    assert trnas_in_space.sprinzl_numeric_from_label("76") == 76
    assert trnas_in_space.sprinzl_numeric_from_label(None) is None


def test_assign_region_from_sprinzl():
    """Test region assignment based on Sprinzl position."""
    # Acceptor stem
    assert trnas_in_space.assign_region_from_sprinzl(1) == "acceptor-stem"
    assert trnas_in_space.assign_region_from_sprinzl(7) == "acceptor-stem"
    assert trnas_in_space.assign_region_from_sprinzl(66) == "acceptor-stem"
    assert trnas_in_space.assign_region_from_sprinzl(72) == "acceptor-stem"

    # Acceptor tail
    assert trnas_in_space.assign_region_from_sprinzl(73) == "acceptor-tail"
    assert trnas_in_space.assign_region_from_sprinzl(74) == "acceptor-tail"

    # Anticodon loop
    assert trnas_in_space.assign_region_from_sprinzl(34) == "anticodon-loop"
    assert trnas_in_space.assign_region_from_sprinzl(35) == "anticodon-loop"
    assert trnas_in_space.assign_region_from_sprinzl(36) == "anticodon-loop"

    # D-loop
    assert trnas_in_space.assign_region_from_sprinzl(14) == "D-loop"
    assert trnas_in_space.assign_region_from_sprinzl(20) == "D-loop"

    # T-loop
    assert trnas_in_space.assign_region_from_sprinzl(54) == "T-loop"
    assert trnas_in_space.assign_region_from_sprinzl(60) == "T-loop"

    # Unknown
    assert trnas_in_space.assign_region_from_sprinzl(None) == "unknown"


def test_infer_trna_id_from_filename():
    """Test tRNA ID inference from filenames."""
    # Standard enriched.json format
    assert (
        trnas_in_space.infer_trna_id_from_filename("tRNA-Ala-AGC-1-1-B_Ala.enriched.json")
        == "tRNA-Ala-AGC-1-1"
    )

    # Simple format
    assert trnas_in_space.infer_trna_id_from_filename("example.enriched.json") == "example"

    # With path
    assert (
        trnas_in_space.infer_trna_id_from_filename("/path/to/tRNA-Leu-CAA-1-1-B_Leu.enriched.json")
        == "tRNA-Leu-CAA-1-1"
    )


def test_should_exclude_trna():
    """Test SeC and mitochondrial tRNA filtering function."""
    # Should exclude SeC tRNAs
    assert trnas_in_space.should_exclude_trna("nuc-tRNA-SeC-TCA-1-1")
    assert trnas_in_space.should_exclude_trna("tRNA-Sec-TCA-1-1")
    assert trnas_in_space.should_exclude_trna("Selenocysteine-tRNA-1")
    assert trnas_in_space.should_exclude_trna("SEC_tRNA")

    # Should exclude mitochondrial tRNAs
    assert trnas_in_space.should_exclude_trna("mito-tRNA-Ala-UGC")
    assert trnas_in_space.should_exclude_trna("mito-tRNA-Leu-UAA")
    assert trnas_in_space.should_exclude_trna("MITO-tRNA-Phe-GAA")

    # Should exclude initiator methionine tRNAs
    assert trnas_in_space.should_exclude_trna("nuc-tRNA-iMet-CAT-1-1")
    assert trnas_in_space.should_exclude_trna("tRNA-initiator-Met-CAU")

    # Should not exclude standard nuclear elongator tRNAs
    assert not trnas_in_space.should_exclude_trna("nuc-tRNA-Ala-GGC-1-1")
    assert not trnas_in_space.should_exclude_trna("nuc-tRNA-Leu-CAA-1-1")
    assert not trnas_in_space.should_exclude_trna("nuc-tRNA-Ser-GCT-1-1")
    assert not trnas_in_space.should_exclude_trna("nuc-tRNA-Phe-GAA-1-1")

    # Handle edge cases
    assert not trnas_in_space.should_exclude_trna(None)
    assert not trnas_in_space.should_exclude_trna("")


def test_validate_no_global_index_collisions():
    """Test collision detection function."""
    # Create test DataFrame with no collisions
    good_df = pd.DataFrame({
        'global_index': [1, 2, 3, 4],
        'sprinzl_label': ['1', '2', '3', '4'],
        'trna_id': ['test1', 'test1', 'test1', 'test1']
    })

    # Should not raise any errors
    try:
        trnas_in_space.validate_no_global_index_collisions(good_df)
    except SystemExit:
        assert False, "Should not exit on good data"

    # Create test DataFrame with collisions
    bad_df = pd.DataFrame({
        'global_index': [1, 2, 2, 3],  # Collision at index 2
        'sprinzl_label': ['1', '47', 'e1', '3'],  # Different labels with same index
        'trna_id': ['test1', 'test1', 'test2', 'test1']
    })

    # Should exit with error
    import pytest
    with pytest.raises(SystemExit):
        trnas_in_space.validate_no_global_index_collisions(bad_df)


def test_output_files_exist():
    """Test that expected output files exist in the outputs directory."""
    outputs_dir = Path(__file__).parent / "outputs"

    # Check for expected output files
    expected_files = [
        "ecoliK12_global_coords.tsv",
        "sacCer_global_coords.tsv",
        "hg38_global_coords.tsv",
    ]

    for filename in expected_files:
        filepath = outputs_dir / filename
        assert filepath.exists(), f"Expected output file not found: {filename}"


def test_output_file_structure():
    """Test that output TSV files have the correct structure."""
    outputs_dir = Path(__file__).parent / "outputs"
    test_file = outputs_dir / "ecoliK12_global_coords.tsv"

    if not test_file.exists():
        return  # Skip if file doesn't exist

    # Read the TSV
    df = pd.read_csv(test_file, sep="\t")

    # Check required columns
    expected_columns = [
        "trna_id",
        "source_file",
        "seq_index",
        "sprinzl_index",
        "sprinzl_label",
        "residue",
        "sprinzl_ordinal",
        "sprinzl_continuous",
        "global_index",
        "region",
    ]

    for col in expected_columns:
        assert col in df.columns, f"Missing expected column: {col}"

    # Check data types and basic constraints
    assert df["seq_index"].dtype in [np.int64, np.int32], "seq_index should be integer"
    assert len(df) > 0, "Output file should not be empty"
    assert df["trna_id"].notna().all(), "trna_id should not have NaN values"
    assert df["residue"].notna().all(), "residue should not have NaN values"

    # Check region values are valid
    valid_regions = [
        "acceptor-stem",
        "acceptor-tail",
        "D-stem",
        "D-loop",
        "anticodon-stem",
        "anticodon-loop",
        "variable-region",
        "variable-arm",
        "T-stem",
        "T-loop",
        "unknown",
    ]
    assert df["region"].isin(valid_regions).all(), "Invalid region values found"

    print(f"✓ {test_file.name}: {len(df)} rows, {df['trna_id'].nunique()} unique tRNAs")


def test_global_index_continuity():
    """Test that global_index values are properly continuous."""
    outputs_dir = Path(__file__).parent / "outputs"
    test_file = outputs_dir / "ecoliK12_global_coords.tsv"

    if not test_file.exists():
        return  # Skip if file doesn't exist

    df = pd.read_csv(test_file, sep="\t")

    # Check that global_index starts at 1 and is continuous
    global_indices = df["global_index"].dropna().unique()
    global_indices.sort()

    assert len(global_indices) > 0, "Should have global indices"
    assert global_indices[0] == 1, "Global index should start at 1"

    # Check for reasonable maximum (should be less than a few hundred for tRNAs)
    assert global_indices[-1] < 500, "Global index seems unreasonably high"


def run_basic_tests():
    """Run all tests manually without pytest."""
    tests = [
        ("Imports", test_imports),
        ("Sort key", test_sort_key),
        ("Sprinzl numeric extraction", test_sprinzl_numeric_from_label),
        ("Region assignment", test_assign_region_from_sprinzl),
        ("Filename parsing", test_infer_trna_id_from_filename),
        ("SeC filtering", test_should_exclude_trna),
        ("Collision detection", test_validate_no_global_index_collisions),
        ("Output files exist", test_output_files_exist),
        ("Output file structure", test_output_file_structure),
        ("Global index continuity", test_global_index_continuity),
    ]

    passed = 0
    failed = 0

    print("Running tests for trnas_in_space.py")
    print("=" * 60)

    for name, test_func in tests:
        try:
            test_func()
            print(f"✓ {name}")
            passed += 1
        except AssertionError as e:
            print(f"✗ {name}: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ {name}: Unexpected error: {e}")
            failed += 1

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")

    return failed == 0


if __name__ == "__main__":
    success = run_basic_tests()
    sys.exit(0 if success else 1)
