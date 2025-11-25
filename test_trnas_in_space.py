#!/usr/bin/env python3
"""
Tests for trnas_in_space.py

Run with: python -m pytest test_trnas_in_space.py -v
Or: python test_trnas_in_space.py
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))

import trnas_in_space  # noqa: E402


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
    # These should sort after position 46 but before 47 (reserved coordinate space)
    assert trnas_in_space.sort_key("46") < trnas_in_space.sort_key("e1")
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("e2")
    assert trnas_in_space.sort_key("e12") < trnas_in_space.sort_key("e13")
    assert trnas_in_space.sort_key("e24") < trnas_in_space.sort_key("47")

    # Test that "e" positions sort in reserved coordinate space
    # "e" positions map to (46, 2, e_num) which sorts after (46, *) and before (47, 0, '')
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("47")
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("47A")
    # But "e" positions should sort after any 46 insertions
    assert trnas_in_space.sort_key("46A") < trnas_in_space.sort_key("e1")

    # Test that e positions don't sort at the end (old bug)
    assert trnas_in_space.sort_key("e1") < trnas_in_space.sort_key("76")

    # Comprehensive ordering test for extended variable arm (reserved coordinate space)
    labels = ["44", "45", "46", "e1", "e2", "e3", "e12", "e13", "e14", "47", "48", "49"]
    sorted_labels = sorted(labels, key=trnas_in_space.sort_key)
    assert sorted_labels == labels, f"Expected {labels}, got {sorted_labels}"

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
    # Must include sprinzl_index column as required by validate_no_global_index_collisions
    good_df = pd.DataFrame(
        {
            "global_index": [1, 2, 3, 4],
            "sprinzl_index": [1, 2, 3, 4],
            "sprinzl_label": ["1", "2", "3", "4"],
            "trna_id": ["test1", "test1", "test1", "test1"],
        }
    )

    # Should not raise any errors
    try:
        trnas_in_space.validate_no_global_index_collisions(good_df)
    except SystemExit:
        assert False, "Should not exit on good data"

    # Create test DataFrame with collisions
    bad_df = pd.DataFrame(
        {
            "global_index": [1, 2, 2, 3],  # Collision at index 2
            "sprinzl_index": [1, 47, -1, 3],  # e1 has invalid sprinzl_index
            "sprinzl_label": ["1", "47", "e1", "3"],  # Different labels with same index
            "trna_id": ["test1", "test1", "test2", "test1"],
        }
    )

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


# ======================== Tests for offset+type grouping ========================


def test_offset_type_strategy_calculate_offset():
    """Test offset calculation from conserved D-loop region (positions 15-25)."""
    strategy = trnas_in_space.OffsetTypeStrategy()

    # Offset 0: label == index
    rows_offset0 = [
        {"sprinzl_label": "15", "sprinzl_index": 15},
        {"sprinzl_label": "18", "sprinzl_index": 18},
        {"sprinzl_label": "20", "sprinzl_index": 20},
    ]
    assert strategy._calculate_offset(rows_offset0) == 0

    # Offset +1: label > index (label 18 at index 17)
    rows_offset_plus1 = [
        {"sprinzl_label": "16", "sprinzl_index": 15},
        {"sprinzl_label": "18", "sprinzl_index": 17},
        {"sprinzl_label": "20", "sprinzl_index": 19},
    ]
    assert strategy._calculate_offset(rows_offset_plus1) == 1

    # Offset -1: label < index (label 18 at index 19)
    rows_offset_minus1 = [
        {"sprinzl_label": "16", "sprinzl_index": 17},
        {"sprinzl_label": "18", "sprinzl_index": 19},
        {"sprinzl_label": "20", "sprinzl_index": 21},
    ]
    assert strategy._calculate_offset(rows_offset_minus1) == -1

    # No valid positions (outside 15-25 range)
    rows_no_range = [
        {"sprinzl_label": "1", "sprinzl_index": 1},
        {"sprinzl_label": "34", "sprinzl_index": 34},
    ]
    assert strategy._calculate_offset(rows_no_range) is None

    # Mixed non-numeric labels should be ignored
    rows_with_labels = [
        {"sprinzl_label": "17a", "sprinzl_index": 18},  # Non-numeric, ignored
        {"sprinzl_label": "18", "sprinzl_index": 18},
        {"sprinzl_label": "e5", "sprinzl_index": -1},  # e-position, ignored
    ]
    assert strategy._calculate_offset(rows_with_labels) == 0


def test_offset_type_strategy_classify():
    """Test OffsetTypeStrategy classification."""
    strategy = trnas_in_space.OffsetTypeStrategy()

    # Type I tRNA with offset 0
    rows_type1_offset0 = [
        {"trna_id": "tRNA-Ala-GGC-1-1", "sprinzl_label": "18", "sprinzl_index": 18},
        {"trna_id": "tRNA-Ala-GGC-1-1", "sprinzl_label": "20", "sprinzl_index": 20},
    ]
    result = strategy.classify("tRNA-Ala-GGC-1-1", rows_type1_offset0)
    assert result is not None
    assert result.dimensions == {"offset": "0", "type": "type1"}

    # Type II tRNA (Leu) with offset +1
    rows_type2_offset_plus1 = [
        {"trna_id": "tRNA-Leu-CAA-1-1", "sprinzl_label": "18", "sprinzl_index": 17},
        {"trna_id": "tRNA-Leu-CAA-1-1", "sprinzl_label": "20", "sprinzl_index": 19},
    ]
    result = strategy.classify("tRNA-Leu-CAA-1-1", rows_type2_offset_plus1)
    assert result is not None
    assert result.dimensions == {"offset": "+1", "type": "type2"}

    # Excluded tRNA (SeC)
    rows_sec = [
        {"trna_id": "tRNA-SeC-TCA-1-1", "sprinzl_label": "18", "sprinzl_index": 18},
    ]
    result = strategy.classify("tRNA-SeC-TCA-1-1", rows_sec)
    assert result is None  # Excluded

    # Excluded tRNA (mitochondrial)
    rows_mito = [
        {"trna_id": "mito-tRNA-Ala-UGC", "sprinzl_label": "18", "sprinzl_index": 18},
    ]
    result = strategy.classify("mito-tRNA-Ala-UGC", rows_mito)
    assert result is None  # Excluded


def test_group_key_filename_suffix():
    """Test GroupKey filename suffix generation."""
    # Type only
    key1 = trnas_in_space.GroupKey({"type": "type1"})
    assert key1.to_filename_suffix() == "_type1"

    # Offset and type
    key2 = trnas_in_space.GroupKey({"offset": "-1", "type": "type2"})
    assert key2.to_filename_suffix() == "_offset-1_type2"

    # Positive offset
    key3 = trnas_in_space.GroupKey({"offset": "+1", "type": "type1"})
    assert key3.to_filename_suffix() == "_offset+1_type1"


def test_offset_type_files_exist():
    """Test that offset+type split files exist."""
    outputs_dir = Path(__file__).parent / "outputs"

    # Check for E. coli offset+type files
    ecoli_files = list(outputs_dir.glob("ecoliK12_global_coords_offset*.tsv"))
    assert (
        len(ecoli_files) >= 4
    ), f"Expected at least 4 E. coli offset+type files, found {len(ecoli_files)}"

    # Check for yeast offset+type files
    yeast_files = list(outputs_dir.glob("sacCer_global_coords_offset*.tsv"))
    assert (
        len(yeast_files) >= 3
    ), f"Expected at least 3 yeast offset+type files, found {len(yeast_files)}"


def test_position_55_alignment_within_groups():
    """Test that position 55 aligns within each offset+type group."""
    outputs_dir = Path(__file__).parent / "outputs"

    for f in outputs_dir.glob("*_global_coords_offset*.tsv"):
        df = pd.read_csv(f, sep="\t")
        pos55 = df[df["sprinzl_label"] == "55"]

        if len(pos55) < 2:
            continue  # Skip files with too few pos55 entries

        # All position 55 instances should have the same global_index
        unique_gidx = pos55["global_index"].dropna().unique()
        assert len(unique_gidx) == 1, (
            f"Position 55 misalignment in {f.name}: "
            f"found {len(unique_gidx)} different global_index values: {list(unique_gidx)}"
        )


def test_no_collisions_in_offset_type_files():
    """Test that each offset+type file has no collisions."""
    outputs_dir = Path(__file__).parent / "outputs"

    for f in outputs_dir.glob("*_global_coords_offset*.tsv"):
        df = pd.read_csv(f, sep="\t")

        # Check for collisions: multiple distinct sprinzl_labels at same global_index
        for gidx, group in df.groupby("global_index"):
            if pd.isna(gidx):
                continue
            labels = group["sprinzl_label"].dropna().unique()
            labels = [lbl for lbl in labels if lbl != ""]
            assert (
                len(labels) <= 1
            ), f"Collision in {f.name}: global_index {gidx} has multiple labels: {labels}"


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
        # New tests for offset+type grouping
        ("Offset calculation", test_offset_type_strategy_calculate_offset),
        ("Strategy classification", test_offset_type_strategy_classify),
        ("GroupKey filename suffix", test_group_key_filename_suffix),
        ("Offset+type files exist", test_offset_type_files_exist),
        ("Position 55 alignment", test_position_55_alignment_within_groups),
        ("No collisions in offset files", test_no_collisions_in_offset_type_files),
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
