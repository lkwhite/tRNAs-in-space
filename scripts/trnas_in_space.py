#!/usr/bin/env python3
"""
trnas_in_space.py

Basic pipeline:
- Collect per-base Sprinzl labels from R2DT *.enriched.json
- Fill missing sprinzl_index values by inference along each tRNA's sequence
- Build shared global label order -> sprinzl_ordinal
- Compute sprinzl_continuous (per-tRNA fractional coords to 6 decimals)
- Map to equal-spaced global_index (1..K)
- Add tRNA region annotation from Sprinzl

Usage:
  python trnas_in_space.py /path/to/r2dt_output_dir out.tsv
"""

import argparse
import json
import os
import re
import sys
from glob import glob

import numpy as np
import pandas as pd

# ----------------------------- config -----------------------------
PRECISION = 6  # fixed rounding for sprinzl_continuous before uniquing

# Manual label overrides - kept as fallback for any edge cases where R2DT
# provides incorrect sprinzl_labels that need manual correction.
# Format: {"trna_id": {seq_index: "correct_label", ...}, ...}
LABEL_OVERRIDES = {}

# tRNAs excluded due to poor R2DT annotation quality.
# These cannot be reliably aligned because key positions lack sprinzl_labels
# or have incorrect T-loop sequences (indicating shifted annotations).
EXCLUDED_POORLY_ANNOTATED = {
    # ===== Nuclear tRNAs =====
    # Missing anticodon labels (positions 34-36) - can't verify tRNA identity
    "nuc-tRNA-Leu-CAA-5-1",   # human, 59.5% empty labels
    # Known R2DT annotation error - anticodon annotated as GGT instead of ATA
    "nuc-tRNA-Tyr-AUA-1-1",   # human
    # Non-standard T-loop sequences that are NOT *TC variants (R2DT errors)
    # Note: *TC variants (CTC, ATC, GTC) are valid biological variants, not errors
    "nuc-tRNA-Gln-UUG-4-1",   # human, T-loop=CGA (not *TC pattern)
    "nuc-tRNA-Trp-CCA-4-1",   # human, T-loop=GCG (not *TC pattern)
    "nuc-tRNA-Trp-CCA-5-1",   # human, T-loop=GCG (not *TC pattern)
    # Broken R2DT alignment - missing sprinzl 1-8 and 65-72, has unique insertions
    "nuc-tRNA-Asp-GUC-2-1",   # yeast, 14.9% empty labels, missing acceptor stem
    # ===== Mitochondrial tRNAs =====
    # Large gaps in R2DT annotation - positions 13-22 have no sprinzl_labels
    "mito-tRNA-Asn-GUU",      # yeast, 40% empty labels (30/75 positions unlabeled)
}

# R2DT label offset corrections for mito tRNAs with shifted Sprinzl labels.
# These tRNAs have correctly sequenced bases but the R2DT template alignment
# assigned incorrect Sprinzl position labels. The offset is added to numeric
# sprinzl_labels to correct the alignment.
# Positive correction = R2DT labels are too low (add to shift up)
# Negative correction = R2DT labels are too high (subtract to shift down)
# Keys must match the trna_id as returned by infer_trna_id_from_filename()
# Human mito tRNAs have -1-1 suffix, yeast mito tRNAs have no copy number suffix
MITO_LABEL_OFFSET_CORRECTIONS = {
    # ===== Human Mito tRNAs (hg38) =====
    # Anticodon found at labels 35-36-37 instead of 34-35-36 (need to subtract 1)
    "mito-tRNA-Leu-UAA-1-1": -1,
    "mito-tRNA-Leu-UAG-1-1": -1,
    "mito-tRNA-Ser-GCU-1-1": -1,
    # Anticodon found at labels 33-34-35 instead of 34-35-36 (need to add 1)
    "mito-tRNA-Lys-UUU-1-1": +1,
    # ===== Yeast Mito tRNAs (sacCer, no -1-1 suffix) =====
    # Anticodon found at labels 36-37-38 instead of 34-35-36 (need to subtract 2)
    # Note: Human Glu-UUC-1-1 is NOT in this dict because it doesn't need correction
    "mito-tRNA-Glu-UUC": -2,
    # Anticodon found at labels 35-36-37 instead of 34-35-36 (need to subtract 1)
    # Note: Human Leu-UAA-1-1 IS in this dict (different key due to -1-1 suffix)
    "mito-tRNA-Leu-UAA": -1,
    # Anticodon found at labels 33-34-35 instead of 34-35-36 (need to add 1)
    "mito-tRNA-Lys-UUU": +1,
}

# Biological order for variable loop hairpin (Type II tRNAs)
# Complete set e1-e27 in 5'→3' order along the RNA backbone:
# - Enter hairpin after position 45
# - Ascend one stem (e11→e17)
# - Cross loop apex (e1→e5)
# - Descend other stem (e27→e21)
# - Additional positions e6-e10 and e18-e20 placed in logical positions
E_POSITION_BIOLOGICAL_ORDER = [
    'e11', 'e12', 'e13', 'e14', 'e15', 'e16', 'e17',  # ascending stem
    'e18', 'e19', 'e20',                               # upper ascending (if present)
    'e1', 'e2', 'e3', 'e4', 'e5',                      # loop apex
    'e6', 'e7', 'e8', 'e9', 'e10',                     # lower descending (if present)
    'e27', 'e26', 'e25', 'e24', 'e23', 'e22', 'e21',  # descending stem
]
E_POSITION_ORDER_MAP = {label: idx for idx, label in enumerate(E_POSITION_BIOLOGICAL_ORDER)}


# ------------------------- helpers: files -------------------------


def infer_trna_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    name = base
    for suf in (".enriched.json", ".json"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    m = re.match(r"^(.*?)-[A-Z]_[A-Za-z0-9]+$", name)  # strip "-B_His" style suffixes if present
    trna_id = m.group(1) if m else name
    # Convert DNA notation (T) to RNA notation (U) in anticodon portion
    # Anticodon is after amino acid: nuc-tRNA-Ala-TGC-1-1 -> nuc-tRNA-Ala-UGC-1-1
    # Also handles: tRNA-Ile2-CAT-1-1, tRNA-fMet-CAT-1-1, etc.
    # Pattern: tRNA-<amino>-<anticodon>- where amino can include digits (Ile2) or lowercase (fMet, iMet)
    def replace_anticodon_t_with_u(match):
        return match.group(1) + match.group(2).replace("T", "U") + match.group(3)
    trna_id = re.sub(r"(tRNA-[A-Za-z0-9]+-)([ACGTU]{3})(-)", replace_anticodon_t_with_u, trna_id, flags=re.IGNORECASE)
    return trna_id


def is_mitochondrial_trna(trna_id: str) -> bool:
    """Check if a tRNA is mitochondrial based on its ID."""
    if trna_id is None:
        return False
    trna_id_upper = trna_id.upper()
    return "MITO-TRNA" in trna_id_upper or trna_id_upper.startswith("MITO-")


def get_label_offset_correction(trna_id: str) -> int:
    """
    Get the Sprinzl label offset correction for a tRNA.

    Some mito tRNAs have R2DT template labels that are shifted relative to
    standard Sprinzl numbering. This function returns the correction to apply
    to numeric labels to fix the alignment.

    Args:
        trna_id: The tRNA identifier (as returned by infer_trna_id_from_filename)

    Returns:
        Integer offset to add to numeric sprinzl_labels (0 if no correction needed)
    """
    if trna_id is None:
        return 0

    # Exact match only - keys in MITO_LABEL_OFFSET_CORRECTIONS must match
    # the trna_id exactly as returned by infer_trna_id_from_filename()
    return MITO_LABEL_OFFSET_CORRECTIONS.get(trna_id, 0)


def apply_label_offset(label: str, offset: int) -> str:
    """
    Apply an offset correction to a Sprinzl label.

    Only adjusts purely numeric labels (e.g., "34" -> "33" with offset -1).
    Non-numeric labels (e.g., "e5", "17a") are returned unchanged.

    Args:
        label: The original Sprinzl label
        offset: The offset to apply

    Returns:
        Corrected label string
    """
    if not label or offset == 0:
        return label

    # Only adjust purely numeric labels
    if label.isdigit():
        new_val = int(label) + offset
        if new_val >= 1:
            return str(new_val)

    return label


def should_exclude_trna(trna_id: str, include_mito: bool = False) -> bool:
    """
    Filter out structurally incompatible tRNAs that cannot be meaningfully aligned.

    Args:
        trna_id: The tRNA identifier
        include_mito: If True, include mitochondrial tRNAs (for mito-only coordinate generation)

    Exclusions for nuclear coordinates:
    1. SeC tRNAs: Structurally distinctive (~95 nt, extended variable arm)
    2. Mitochondrial tRNAs: Different architecture (60-75 nt, variable structure)
    3. Initiator Met (iMet/fMet): Special structures for ribosome binding
    4. Poorly annotated tRNAs: Missing critical sprinzl_labels from R2DT

    Exclusions for mitochondrial coordinates:
    1. Nuclear tRNAs: Different coordinate system
    2. SeC tRNAs: Not found in mitochondria
    3. Poorly annotated tRNAs: Missing critical labels
    """
    if trna_id is None:
        return False

    trna_id_upper = trna_id.upper()
    is_mito = is_mitochondrial_trna(trna_id)

    if include_mito:
        # Generating mitochondrial coordinates - only include mito tRNAs
        if not is_mito:
            return True  # Exclude non-mito tRNAs
        # Check poorly annotated exclusion list (includes some mito tRNAs)
        if trna_id in EXCLUDED_POORLY_ANNOTATED:
            return True
        return False
    else:
        # Generating nuclear coordinates - exclude mito tRNAs
        # Check poorly annotated exclusion list first
        if trna_id in EXCLUDED_POORLY_ANNOTATED:
            return True

        # Exclude selenocysteine - incompatible structure
        if "SEC" in trna_id_upper or "SELENOCYSTEINE" in trna_id_upper:
            return True

        # Exclude mitochondrial tRNAs - different structural architecture
        if is_mito:
            return True

        # Exclude initiator methionine tRNAs - different structural features
        if "IMET" in trna_id_upper or "INITIAT" in trna_id_upper or "FMET" in trna_id_upper:
            return True

        return False


def classify_trna_type(trna_id: str, include_mito: bool = False) -> str:
    """
    Classify tRNAs into structural types for dual coordinate system approach.

    Args:
        trna_id: The tRNA identifier
        include_mito: If True, classifying for mito coordinate generation

    Returns:
        'type1': Standard tRNAs with simple variable arm (most amino acids)
        'type2': Extended variable arm tRNAs (Leu, Ser, Tyr with e-positions)
        'exclude': Structurally incompatible tRNAs (SeC, mito, iMet)
    """
    # First check if this tRNA should be excluded entirely
    if should_exclude_trna(trna_id, include_mito=include_mito):
        return "exclude"

    if trna_id is None:
        return "exclude"

    trna_id_upper = trna_id.upper()

    # Type II: Extended variable arm tRNAs (Leu, Ser, Tyr)
    # These have e1-e24 positions in their extended variable arms
    if any(amino_acid in trna_id_upper for amino_acid in ["LEU", "SER", "TYR"]):
        return "type2"

    # Type I: All other nuclear elongator tRNAs
    # Standard 76nt structure with simple variable arm (positions 47-48)
    return "type1"


# --------------------- validation ---------------------

from typing import List, Optional


def validate_no_global_index_collisions(df: pd.DataFrame):
    """
    Detect and report global_index collisions that would break coordinate alignment.
    Uses the preferred label (the one actually used for coordinates) rather than
    raw sprinzl_label to avoid false positives from R2DT template differences.
    Exits with error if collisions are found.
    """
    # Build preferred labels using the same logic as coordinate generation
    pref_labels = build_pref_label(df)
    df_with_pref = df.copy()
    df_with_pref["pref_label"] = pref_labels

    # Group by global_index and find any with multiple distinct preferred labels
    collision_groups = []
    for global_idx, group in df_with_pref.groupby("global_index"):
        if pd.isna(global_idx):
            continue
        unique_pref_labels = group["pref_label"].unique()
        unique_pref_labels = [lbl for lbl in unique_pref_labels if pd.notna(lbl) and lbl != ""]
        if len(unique_pref_labels) > 1:
            # Found collision: multiple preferred labels mapping to same global_index
            collision_groups.append((global_idx, unique_pref_labels, group))

    if collision_groups:
        print("\n[ERROR] Global index collisions detected!")
        print(
            "Multiple structural positions share the same global_index, breaking coordinate alignment:\n"
        )
        for global_idx, labels, group in collision_groups:
            print(f"  global_index {global_idx}:")
            for label in labels:
                trna_examples = group[group["pref_label"] == label]["trna_id"].unique()[:3]
                print(f"    - position '{label}' (examples: {', '.join(trna_examples)})")
            print()

        print("This indicates that the coordinate system cannot properly handle the diversity")
        print("of tRNA structures present. Consider:")
        print("  1. Filtering out additional incompatible tRNA types")
        print("  2. Adjusting the sort_key() function for better position ordering")
        print("  3. Expanding reserved coordinate space for extended variable arms")
        sys.exit(1)

    print(
        f"[ok] Global index validation: No collisions detected among {df['global_index'].nunique()} unique positions"
    )


# --------------------- phase 1: JSON -> rows ----------------------


def collect_rows_from_json(fp: str, include_mito: bool = False):
    """
    Collect tRNA data rows from an R2DT enriched JSON file.

    Args:
        fp: Path to the enriched JSON file
        include_mito: If True, collecting for mito coordinates (include mito, exclude nuclear)
                      If False, collecting for nuclear coordinates (include nuclear, exclude mito)
    """
    with open(fp, "r") as f:
        J = json.load(f)
    mol = J["rnaComplexes"][0]["rnaMolecules"][0]
    seq = mol["sequence"]

    trna_id = infer_trna_id_from_filename(fp)

    # Filter based on tRNA type and mode
    if should_exclude_trna(trna_id, include_mito=include_mito):
        # Only print for significant exclusions, not mode-based filtering
        is_mito = is_mitochondrial_trna(trna_id)
        if include_mito and not is_mito:
            pass  # Don't log nuclear tRNAs excluded in mito mode
        elif not include_mito and is_mito:
            pass  # Don't log mito tRNAs excluded in nuclear mode
        else:
            print(f"Excluding incompatible tRNA: {trna_id}")
        return []

    # Check if this tRNA needs label offset correction
    label_offset = get_label_offset_correction(trna_id)
    if label_offset != 0:
        print(f"Applying label offset correction of {label_offset:+d} to {trna_id}")

    rows = []
    for s in seq:
        rname = s.get("residueName")
        if rname in ("5'", "3'"):
            continue
        ridx = int(s["residueIndex"])  # 1-based
        info = s.get("info", {}) or {}
        sprinzl_idx = info.get("templateResidueIndex", None)
        sprinzl_idx = int(sprinzl_idx) if isinstance(sprinzl_idx, int) else -1
        sprinzl_lbl = (info.get("templateNumberingLabel", "") or "").strip()

        # Apply label offset correction for mito tRNAs with shifted R2DT labels
        if label_offset != 0:
            sprinzl_lbl = apply_label_offset(sprinzl_lbl, label_offset)
            # Also adjust the sprinzl_index if it's valid
            if sprinzl_idx >= 1:
                sprinzl_idx = sprinzl_idx + label_offset
                if sprinzl_idx < 1:
                    sprinzl_idx = -1

        rows.append(
            {
                "trna_id": trna_id,
                "source_file": os.path.basename(fp),
                "seq_index": ridx,
                "sprinzl_index": sprinzl_idx,
                "sprinzl_label": sprinzl_lbl,
                "residue": rname,
            }
        )

    # Fill missing sprinzl_index by monotone inference along seq_index
    rows.sort(key=lambda r: r["seq_index"])
    n = len(rows)
    vals = [r["sprinzl_index"] for r in rows]

    # forward pass
    fwd = [None] * n
    last = None
    for i in range(n):
        v = vals[i]
        if isinstance(v, int) and v >= 1:
            fwd[i] = v
            last = v
        elif last is not None:
            fwd[i] = last + 1
            last += 1

    # backward pass
    bwd = [None] * n
    nxt = None
    for i in range(n - 1, -1, -1):
        v = vals[i]
        if isinstance(v, int) and v >= 1:
            bwd[i] = v
            nxt = v
        elif nxt is not None:
            bwd[i] = nxt - 1
            nxt -= 1

    for i, r in enumerate(rows):
        v = vals[i]
        if v >= 1:
            r["sprinzl_index"] = v
        elif fwd[i] is not None and bwd[i] is not None and fwd[i] == bwd[i] and 1 <= fwd[i] <= 76:
            r["sprinzl_index"] = fwd[i]
        elif fwd[i] is not None and 1 <= fwd[i] <= 76 and bwd[i] is None:
            r["sprinzl_index"] = fwd[i]
        elif bwd[i] is not None and 1 <= bwd[i] <= 76 and fwd[i] is None:
            r["sprinzl_index"] = bwd[i]
        else:
            r["sprinzl_index"] = -1

    # Apply label overrides for known R2DT labeling errors (manual fallback)
    if trna_id in LABEL_OVERRIDES:
        overrides = LABEL_OVERRIDES[trna_id]
        for row in rows:
            if row["seq_index"] in overrides:
                row["sprinzl_label"] = overrides[row["seq_index"]]

    return rows


# --------------- phase 2: label order & continuous ----------------


def normalize_label(lbl: str) -> str:
    """
    Normalize a Sprinzl label for consistent processing.

    - Strips whitespace
    - Converts float-formatted labels: "1.0" -> "1"
    - Handles None/empty/nan values

    Returns: Normalized label string, or empty string for invalid inputs.
    """
    if lbl is None or lbl == "nan":
        return ""
    s = str(lbl).strip()
    if s == "" or s == "nan":
        return ""

    # Normalize float-formatted labels: "1.0" -> "1"
    float_match = re.fullmatch(r"(\d+)\.0", s)
    if float_match:
        return float_match.group(1)

    return s


def sort_key(lbl: str):
    """
    Order labels for unified coordinate system.

    Sort order: 20 < 20A < 20B; dotted like 9.1; e-positions in biological hairpin order.

    Returns tuples where ALL third elements are strings to ensure type consistency
    in Python 3 comparisons. Uses zero-padding for numeric values.
    """
    s = normalize_label(lbl)
    if s == "":
        return (10**9, 2, "")

    # Type II extended variable arm positions (e1-e27) - biological hairpin ordering
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        e_label = s
        if e_label in E_POSITION_ORDER_MAP:
            bio_order = E_POSITION_ORDER_MAP[e_label]
            # Zero-pad to 3 digits for string comparison
            return (45, 2, f"{bio_order:03d}")
        else:
            return (45, 3, s)

    # Standard numeric positions with optional letter suffixes
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)

    # Dotted positions like 9.1 - convert to zero-padded string
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, f"{int(m.group(2)):03d}")

    return (10**9 - 1, 2, s)


def sort_key_type1(lbl: str):
    """
    Sort key for Type I tRNAs: Standard 76nt tRNAs with simple variable arm.
    No e-positions expected, simpler coordinate space.

    Returns tuples where ALL third elements are strings for type consistency.
    """
    s = normalize_label(lbl)
    if s == "":
        return (10**9, 2, "")

    # Standard numeric positions with optional letter suffixes
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)

    # Allow dotted positions like 9.1 - zero-pad for string comparison
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, f"{int(m.group(2)):03d}")

    # e-positions should not occur in Type I, but handle gracefully
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        print(f"Warning: e-position {s} found in Type I tRNA (unexpected)")
        return (10**8, 2, f"{int(m.group(1)):03d}")

    # Unknown labels at end
    return (10**9 - 1, 2, s)


def sort_key_type2(lbl: str):
    """
    Sort key for Type II tRNAs: Extended variable arm tRNAs (Leu, Ser, Tyr).
    Optimized for e-position handling with biological hairpin ordering.

    Returns tuples where ALL third elements are strings to ensure type consistency
    in Python 3 comparisons. Uses zero-padding for numeric values.
    """
    s = normalize_label(lbl)
    if s == "":
        return (10**9, 2, "")

    # Type II extended variable arm positions (e1-e27) - biological hairpin ordering
    # Uses E_POSITION_ORDER_MAP to sort in 5'→3' order along the RNA backbone
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        e_label = s  # e.g., "e1", "e12"
        if e_label in E_POSITION_ORDER_MAP:
            bio_order = E_POSITION_ORDER_MAP[e_label]
            # Place e-positions after position 45, in biological order
            # Zero-pad to 3 digits for string comparison
            return (45, 2, f"{bio_order:03d}")
        else:
            # Unknown e-position - place at end of e-region
            print(f"Warning: Unknown e-position {e_label} not in biological order map")
            return (45, 3, s)

    # Standard numeric positions with optional letter suffixes
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)

    # Allow dotted positions like 9.1 - convert to zero-padded string
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, f"{int(m.group(2)):03d}")

    # Unknown labels at end
    return (10**9 - 1, 2, s)


def build_pref_label(df: pd.DataFrame) -> pd.Series:
    """
    Use sprinzl_label only (templateNumberingLabel = canonical Sprinzl position).

    This is the correct field for functional alignment as it represents the
    canonical position (e.g., '18' is always position 18 regardless of insertions).

    Empty labels (insertions without canonical positions) are NOT filled with
    sprinzl_index values. Instead, they remain empty and get interpolated
    fractional coordinates in make_continuous_for_trna(). This preserves
    sequence order - insertions get coordinates between their labeled neighbors.

    Previous versions fell back to sprinzl_index for empty labels, but this
    caused ordering violations because inferred indices (e.g., "60") sorted
    incorrectly relative to actual Sprinzl labels (e.g., "52").
    """
    lbl = df["sprinzl_label"].astype("string").fillna("").str.strip()
    return lbl


def build_global_label_order(pref: pd.Series):
    uniq = sorted({p for p in pref if p not in ("", "nan")}, key=sort_key)
    to_ord = {u: i + 1 for i, u in enumerate(uniq)}  # 1..K
    return uniq, to_ord


def make_continuous_for_trna(sub: pd.DataFrame, ord_series: pd.Series) -> pd.Series:
    """
    For one tRNA (sorted by seq_index):
      - labeled sites -> integer ordinals
      - unlabeled internal runs -> fractions between adjacent ordinals
      - leading/trailing runs -> fractions near edge bins
    """
    sub = sub.sort_values("seq_index").copy()
    ords = ord_series.loc[sub.index].astype("Float64").to_numpy()
    n = len(sub)
    vals = [np.nan] * n
    i = 0
    while i < n:
        if not pd.isna(ords[i]):
            vals[i] = float(int(ords[i]))
            i += 1
            continue
        j = i
        while j < n and pd.isna(ords[j]):
            j += 1
        k = j - i
        left = int(ords[i - 1]) if i - 1 >= 0 and not pd.isna(ords[i - 1]) else None
        right = int(ords[j]) if j < n and not pd.isna(ords[j]) else None
        if left is not None and right is not None and right >= left + 1:
            for t in range(k):
                vals[i + t] = left + (t + 1) / (k + 1)
        elif left is None and right is not None:
            for t in range(k):
                vals[i + t] = right - (k - t) / (k + 1)
        elif left is not None and right is None:
            for t in range(k):
                vals[i + t] = left + (t + 1) / (k + 1)
        else:
            for t in range(k):
                vals[i + t] = np.nan
        i = j
    return pd.Series(vals, index=sub.index, dtype="float64")


# ------------------------ phase 3: regions -------------------------


def sprinzl_numeric_from_label(label: str):
    """Extract leading integer from Sprinzl label like '20a', '14:i1' -> 20, 14."""
    if label is None:
        return None
    m = re.match(r"(\d+)", str(label))
    return int(m.group(1)) if m else None


def assign_region_from_sprinzl(base_num: int) -> str:
    """
    Region buckets (Type I canonical; robust to insertions).
    See README text: acceptor-stem, D-stem/loop, anticodon-stem/loop,
    variable-region/arm, T-stem/loop, acceptor-tail.
    """
    if base_num is None:
        return "unknown"
    p = base_num
    if (1 <= p <= 7) or (66 <= p <= 72):
        return "acceptor-stem"
    if p == 73 or p > 73:
        return "acceptor-tail"  # include CCA if present
    if (10 <= p <= 13) or (22 <= p <= 25):
        return "D-stem"
    if 14 <= p <= 21:
        return "D-loop"
    if (27 <= p <= 31) or (39 <= p <= 43):
        return "anticodon-stem"
    if 32 <= p <= 38:
        return "anticodon-loop"
    if 44 <= p <= 46:
        return "variable-region"
    if 46 < p < 49:
        return "variable-arm"  # Type II long arm insertions
    if (49 <= p <= 53) or (61 <= p <= 65):
        return "T-stem"
    if 54 <= p <= 60:
        return "T-loop"
    return "unknown"


def compute_region_column(df: pd.DataFrame) -> pd.Series:
    # prefer label’s numeric part; fall back to sprinzl_index (1..76)
    base_from_label = df["sprinzl_label"].apply(sprinzl_numeric_from_label)
    idx_fallback = pd.to_numeric(df["sprinzl_index"], errors="coerce").where(
        lambda x: (x >= 1) & (x <= 76)
    )
    base_num = base_from_label.fillna(idx_fallback)
    return base_num.apply(assign_region_from_sprinzl)


# -------------------------------- main ---------------------------------


def build_global_label_order_type1(pref: pd.Series):
    """Build global label order for Type I tRNAs using type1-specific sort key."""
    uniq = sorted({p for p in pref if p not in ("", "nan")}, key=sort_key_type1)
    to_ord = {u: i + 1 for i, u in enumerate(uniq)}  # 1..K
    return uniq, to_ord


def build_global_label_order_type2(pref: pd.Series):
    """Build global label order for Type II tRNAs using type2-specific sort key."""
    uniq = sorted({p for p in pref if p not in ("", "nan")}, key=sort_key_type2)
    to_ord = {u: i + 1 for i, u in enumerate(uniq)}  # 1..K
    return uniq, to_ord


def generate_coordinates_for_type(all_rows, trna_type, output_file, allow_collisions=False):
    """
    Generate coordinates for a specific tRNA type (Type I or Type II).

    Args:
        all_rows: List of all tRNA data rows
        trna_type: 'type1' or 'type2'
        output_file: Path to output TSV file
        allow_collisions: Whether to allow coordinate collisions
    """
    # Filter rows to only include the specified type
    filtered_rows = []
    excluded_count = 0

    for row in all_rows:
        classification = classify_trna_type(row["trna_id"])
        if classification == trna_type:
            filtered_rows.append(row)
        else:
            excluded_count += 1

    if not filtered_rows:
        print(f"[error] No {trna_type} tRNAs found in dataset")
        return

    print(
        f"[info] Processing {len(filtered_rows)} {trna_type} tRNAs (excluded {excluded_count} other types)"
    )

    # Convert to DataFrame
    df = pd.DataFrame(filtered_rows).sort_values(["trna_id", "seq_index"]).reset_index(drop=True)

    # Build preferred labels
    pref = build_pref_label(df)

    # Use type-specific label ordering
    if trna_type == "type1":
        uniq_labels, to_ord = build_global_label_order_type1(pref)
    elif trna_type == "type2":
        uniq_labels, to_ord = build_global_label_order_type2(pref)
    else:
        raise ValueError(f"Unknown trna_type: {trna_type}")

    df["sprinzl_ordinal"] = pd.to_numeric(pref.map(to_ord), errors="coerce")

    # Continuous coordinate per tRNA
    ord_series = pref.map(to_ord)
    cont = []
    for _, sub in df.groupby("trna_id", sort=False):
        cont.append(make_continuous_for_trna(sub, ord_series))
    df["sprinzl_continuous"] = pd.concat(cont).sort_index().astype("float64")

    # Global equal-spaced index
    cont_round = df["sprinzl_continuous"].round(PRECISION)
    uniq_cont = sorted(cont_round.dropna().unique().tolist())
    cont_to_global = {v: i + 1 for i, v in enumerate(uniq_cont)}
    df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

    # Validate no collisions (unless allowing them)
    if allow_collisions:
        print(
            f"[info] {trna_type.upper()}: Collision validation bypassed due to --allow-collisions flag"
        )
    else:
        validate_no_global_index_collisions(df)

    # Region annotation
    df["region"] = compute_region_column(df)

    # Write output
    cols = [
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
    df.to_csv(output_file, sep="\t", index=False, columns=cols)

    # Stats
    print(f"[ok] {trna_type.upper()}: Wrote {output_file}")
    print(f"  Rows: {len(df)}  |  tRNAs: {df['trna_id'].nunique()}")
    print(f"  Unique labeled bins: {len(uniq_labels)}")
    print(f"  Unique global positions (K): {len(uniq_cont)}  |  rounding={PRECISION} d.p.")


def main():
    ap = argparse.ArgumentParser(
        description="Build global equal-spaced tRNA coordinates + regions from R2DT JSON."
    )
    ap.add_argument(
        "json_dir", help="Directory with R2DT *.enriched.json files (searched recursively)."
    )
    ap.add_argument("out_tsv", help="Output TSV path (or base name for dual system).")
    ap.add_argument(
        "--allow-collisions",
        action="store_true",
        help="Allow coordinate generation despite isoacceptor-level position collisions.",
    )
    ap.add_argument(
        "--dual-system",
        action="store_true",
        help="Generate separate coordinate files for Type I and Type II tRNAs.",
    )
    ap.add_argument(
        "--type",
        choices=["type1", "type2"],
        help="Generate coordinates for specific tRNA type only (overrides --dual-system).",
    )
    ap.add_argument(
        "--mito",
        action="store_true",
        help="Generate coordinates for mitochondrial tRNAs only (separate from nuclear).",
    )
    args = ap.parse_args()

    paths = glob(os.path.join(args.json_dir, "**", "*.enriched.json"), recursive=True)
    if not paths:
        print(f"[error] No *.enriched.json files found under: {args.json_dir}")
        sys.exit(2)

    # Collect tRNA data - pass include_mito to filter appropriately
    all_rows, skipped = [], 0
    for fp in sorted(paths):
        try:
            all_rows.extend(collect_rows_from_json(fp, include_mito=args.mito))
        except Exception as e:
            skipped += 1
            print(f"[warn] Skipping {fp} due to error: {e}")

    print(f"[info] JSON files parsed: {len(paths)}  |  skipped: {skipped}")

    # Determine which coordinate systems to generate
    if args.mito:
        # Generate coordinates for mitochondrial tRNAs only
        # (filtering already done during collection via include_mito=True)
        print("[info] Generating mitochondrial tRNA coordinate system")

        if not all_rows:
            print("[error] No mitochondrial tRNAs found in dataset")
            sys.exit(2)

        unique_trnas = len(set(r['trna_id'] for r in all_rows))
        print(f"[info] Found {unique_trnas} mitochondrial tRNAs")

        df = (
            pd.DataFrame(all_rows).sort_values(["trna_id", "seq_index"]).reset_index(drop=True)
        )

        # Build global label order
        pref = build_pref_label(df)
        uniq_labels, to_ord = build_global_label_order(pref)
        df["sprinzl_ordinal"] = pd.to_numeric(pref.map(to_ord), errors="coerce")

        # Generate continuous coordinates per-tRNA
        ord_series = pref.map(to_ord)
        cont = []
        for _, sub in df.groupby("trna_id", sort=False):
            cont.append(make_continuous_for_trna(sub, ord_series))
        df["sprinzl_continuous"] = pd.concat(cont).sort_index().astype("float64")

        # Map to integer global_index
        cont_round = df["sprinzl_continuous"].round(PRECISION)
        uniq_cont = sorted(cont_round.dropna().unique().tolist())
        cont_to_global = {v: i + 1 for i, v in enumerate(uniq_cont)}
        df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

        # Validate (collisions less likely in mito due to simpler structure)
        if args.allow_collisions:
            print("[info] Collision validation bypassed due to --allow-collisions flag")
        else:
            validate_no_global_index_collisions(df)

        df["region"] = compute_region_column(df)

        cols = [
            "trna_id", "source_file", "seq_index", "sprinzl_index", "sprinzl_label",
            "residue", "sprinzl_ordinal", "sprinzl_continuous", "global_index", "region",
        ]
        df.to_csv(args.out_tsv, sep="\t", index=False, columns=cols)

        print(f"[ok] Wrote {args.out_tsv}")
        print(f"  Rows: {len(df)}  |  tRNAs: {df['trna_id'].nunique()}")
        print(f"  Unique labeled bins: {len(uniq_labels)}")
        print(f"  Unique global positions (K): {len(uniq_cont)}  |  rounding={PRECISION} d.p.")

    elif args.type:
        # Generate coordinates for specific type only
        output_file = args.out_tsv
        generate_coordinates_for_type(all_rows, args.type, output_file, args.allow_collisions)

    elif args.dual_system:
        # Generate separate coordinate files for both types
        base_name = args.out_tsv
        if base_name.endswith(".tsv"):
            base_name = base_name[:-4]

        type1_file = f"{base_name}_type1.tsv"
        type2_file = f"{base_name}_type2.tsv"

        print("[info] Generating dual coordinate system:")
        print(f"  Type I (standard): {type1_file}")
        print(f"  Type II (extended): {type2_file}")

        generate_coordinates_for_type(all_rows, "type1", type1_file, args.allow_collisions)
        generate_coordinates_for_type(all_rows, "type2", type2_file, args.allow_collisions)

    else:
        # Unified coordinate system - single global_index for all tRNAs
        # This works because:
        # 1. Empty labels no longer get filled with wrong index values
        # 2. Sort keys return consistent types (all strings in third element)
        # 3. E-positions sort in biological (hairpin) order
        print("[info] Generating unified coordinate system")

        # Filter out excluded tRNAs (only process type1 and type2)
        filtered_rows = []
        for row in all_rows:
            classification = classify_trna_type(row["trna_id"])
            if classification in ["type1", "type2"]:
                filtered_rows.append(row)

        df = (
            pd.DataFrame(filtered_rows).sort_values(["trna_id", "seq_index"]).reset_index(drop=True)
        )

        # Build global label order using unified sort_key
        pref = build_pref_label(df)
        uniq_labels, to_ord = build_global_label_order(pref)
        df["sprinzl_ordinal"] = pd.to_numeric(pref.map(to_ord), errors="coerce")

        # Generate continuous coordinates per-tRNA (handles empty labels)
        ord_series = pref.map(to_ord)
        cont = []
        for _, sub in df.groupby("trna_id", sort=False):
            cont.append(make_continuous_for_trna(sub, ord_series))
        df["sprinzl_continuous"] = pd.concat(cont).sort_index().astype("float64")

        # Map continuous values to integer global_index
        cont_round = df["sprinzl_continuous"].round(PRECISION)
        uniq_cont = sorted(cont_round.dropna().unique().tolist())
        cont_to_global = {v: i + 1 for i, v in enumerate(uniq_cont)}
        df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

        # Validate no collisions (same global_index with different residues)
        if args.allow_collisions:
            print("[info] Collision validation bypassed due to --allow-collisions flag")
        else:
            validate_no_global_index_collisions(df)

        df["region"] = compute_region_column(df)

        cols = [
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
        df.to_csv(args.out_tsv, sep="\t", index=False, columns=cols)

        print(f"[ok] Wrote {args.out_tsv}")
        print(f"  Rows: {len(df)}  |  tRNAs: {df['trna_id'].nunique()}")
        print(f"  Unique labeled bins: {len(uniq_labels)}")
        print(f"  Unique global positions (K): {len(uniq_cont)}  |  rounding={PRECISION} d.p.")


if __name__ == "__main__":
    main()
