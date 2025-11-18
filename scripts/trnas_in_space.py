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

import os, sys, json, re, argparse
from glob import glob
import pandas as pd
import numpy as np

# ----------------------------- config -----------------------------
PRECISION = 6  # fixed rounding for sprinzl_continuous before uniquing

# ------------------------- helpers: files -------------------------


def infer_trna_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    name = base
    for suf in (".enriched.json", ".json"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    m = re.match(r"^(.*?)-[A-Z]_[A-Za-z0-9]+$", name)  # strip "-B_His" style suffixes if present
    return m.group(1) if m else name


def should_exclude_trna(trna_id: str) -> bool:
    """
    Filter out structurally incompatible tRNAs that cannot be meaningfully aligned
    with standard nuclear tRNA coordinate systems.

    Exclusions:
    1. SeC tRNAs: Structurally distinctive due to unique biological requirements:
       - Avoid normal stop codon recognition (would terminate translation)
       - Recruit specialized factors (SelA/SelB complex)
       - Coordinate with SECIS elements (distant mRNA signals)
       - Results in ~95 nucleotides vs. standard 76, with 17+ nucleotide extended
         variable arms that cannot be aligned with Type I/II tRNAs.

    2. Mitochondrial tRNAs: Fundamentally different architecture:
       - Often 60-75 nucleotides vs. 76 for nuclear tRNAs
       - Can lack certain structural features (e.g., D-loop)
       - Different Sprinzl numbering patterns
       - Designed for mitochondrial translation system, not cytoplasmic
       - Cannot be meaningfully aligned with nuclear tRNA coordinates

    Focus on nuclear tRNAs (nuc-tRNA-*) which have standardized Type I/II structures.
    """
    if trna_id is None:
        return False

    trna_id_upper = trna_id.upper()

    # Exclude selenocysteine - incompatible structure
    if 'SEC' in trna_id_upper or 'SELENOCYSTEINE' in trna_id_upper:
        return True

    # Exclude mitochondrial tRNAs - different structural architecture
    if 'MITO-TRNA' in trna_id_upper or trna_id_upper.startswith('MITO-'):
        return True

    # Exclude initiator methionine tRNAs - different structural features
    # Initiator tRNAs have modified structures for ribosome binding
    if 'IMET' in trna_id_upper or 'INITIAT' in trna_id_upper or 'FMET' in trna_id_upper:
        return True

    # Add other exclusions as needed
    return False


def classify_trna_type(trna_id: str) -> str:
    """
    Classify tRNAs into structural types for dual coordinate system approach.

    Returns:
        'type1': Standard tRNAs with simple variable arm (most amino acids)
        'type2': Extended variable arm tRNAs (Leu, Ser, Tyr with e-positions)
        'exclude': Structurally incompatible tRNAs (SeC, mito, iMet)
    """
    # First check if this tRNA should be excluded entirely
    if should_exclude_trna(trna_id):
        return 'exclude'

    if trna_id is None:
        return 'exclude'

    trna_id_upper = trna_id.upper()

    # Type II: Extended variable arm tRNAs (Leu, Ser, Tyr)
    # These have e1-e24 positions in their extended variable arms
    if any(amino_acid in trna_id_upper for amino_acid in ['LEU', 'SER', 'TYR']):
        return 'type2'

    # Type I: All other nuclear elongator tRNAs
    # Standard 76nt structure with simple variable arm (positions 47-48)
    return 'type1'


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
    df_with_pref['pref_label'] = pref_labels

    # Group by global_index and find any with multiple distinct preferred labels
    collision_groups = []
    for global_idx, group in df_with_pref.groupby('global_index'):
        if pd.isna(global_idx):
            continue
        unique_pref_labels = group['pref_label'].unique()
        unique_pref_labels = [lbl for lbl in unique_pref_labels if pd.notna(lbl) and lbl != '']
        if len(unique_pref_labels) > 1:
            # Found collision: multiple preferred labels mapping to same global_index
            collision_groups.append((global_idx, unique_pref_labels, group))

    if collision_groups:
        print("\n[ERROR] Global index collisions detected!")
        print("Multiple structural positions share the same global_index, breaking coordinate alignment:\n")
        for global_idx, labels, group in collision_groups:
            print(f"  global_index {global_idx}:")
            for label in labels:
                trna_examples = group[group['pref_label'] == label]['trna_id'].unique()[:3]
                print(f"    - position '{label}' (examples: {', '.join(trna_examples)})")
            print()

        print("This indicates that the coordinate system cannot properly handle the diversity")
        print("of tRNA structures present. Consider:")
        print("  1. Filtering out additional incompatible tRNA types")
        print("  2. Adjusting the sort_key() function for better position ordering")
        print("  3. Expanding reserved coordinate space for extended variable arms")
        sys.exit(1)

    print(f"[ok] Global index validation: No collisions detected among {df['global_index'].nunique()} unique positions")


# --------------------- phase 1: JSON -> rows ----------------------


def collect_rows_from_json(fp: str):
    with open(fp, "r") as f:
        J = json.load(f)
    mol = J["rnaComplexes"][0]["rnaMolecules"][0]
    seq = mol["sequence"]

    trna_id = infer_trna_id_from_filename(fp)

    # Filter out structurally incompatible tRNAs (e.g., SeC)
    if should_exclude_trna(trna_id):
        print(f"Excluding incompatible tRNA: {trna_id}")
        return []

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

    return rows


# --------------- phase 2: label order & continuous ----------------


def sort_key(lbl: str):
    """Order labels: 20 < 20A < 20B; allow dotted like 9.1; Type II extended arms (e1-e24) map to reserved space; unknowns at end."""
    if lbl is None or lbl == "" or lbl == "nan":
        return (10**9, 2, "")
    s = str(lbl)
    # Type II extended variable arm positions (e1-e24) map to reserved coordinate space
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        e_num = int(m.group(1))
        # Map to reserved space: after position 46, before regular 48
        return (46, 2, e_num)
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, int(m.group(2)))
    return (10**9 - 1, 2, s)


def sort_key_type1(lbl: str):
    """
    Sort key for Type I tRNAs: Standard 76nt tRNAs with simple variable arm.
    No e-positions expected, simpler coordinate space.
    """
    if lbl is None or lbl == "" or lbl == "nan":
        return (10**9, 2, "")
    s = str(lbl)

    # Standard numeric positions with optional letter suffixes
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)

    # Allow dotted positions like 9.1
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, int(m.group(2)))

    # e-positions should not occur in Type I, but handle gracefully
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        print(f"Warning: e-position {s} found in Type I tRNA (unexpected)")
        return (10**8, 2, int(m.group(1)))

    # Unknown labels at end
    return (10**9 - 1, 2, s)


def sort_key_type2(lbl: str):
    """
    Sort key for Type II tRNAs: Extended variable arm tRNAs (Leu, Ser, Tyr).
    Optimized for e-position handling without collision concerns.
    """
    if lbl is None or lbl == "" or lbl == "nan":
        return (10**9, 2, "")
    s = str(lbl)

    # Type II extended variable arm positions (e1-e24) - natural ordering
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        e_num = int(m.group(1))
        # Place e-positions after position 46, in natural order
        return (46, 2, e_num)

    # Standard numeric positions with optional letter suffixes
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)

    # Allow dotted positions like 9.1
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, int(m.group(2)))

    # Unknown labels at end
    return (10**9 - 1, 2, s)


def build_pref_label(df: pd.DataFrame) -> pd.Series:
    """
    Prefer sprinzl_index (numeric) for standard positions [1..76] to ensure consistency.
    Use sprinzl_label for structural insertions:
    - Type II extended arms ("e1", "e2", etc.) - crucial structural information
    - Other insertions (e.g., "20a", "47a") where index is invalid
    Returns strings with '' for unknown.
    """
    lbl = df["sprinzl_label"].astype("string").fillna("").str.strip()
    num = pd.to_numeric(df["sprinzl_index"], errors="coerce")
    num_ok = num.where((num >= 1) & (num <= 76))
    num_ok_str = num_ok.astype("Int64").astype(str).replace({"<NA>": ""})

    # Always preserve Type II extended arm positions (e1, e2, e3, etc.)
    is_extended_arm = lbl.str.match(r'^e\d+$', na=False)

    # Use numeric index when valid (1-76), otherwise fall back to label for insertions
    pref = num_ok_str.mask(num_ok_str.eq(""), other=lbl)

    # Override with extended arm labels to preserve structural information
    pref = pref.mask(is_extended_arm, other=lbl)

    return pref


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
        classification = classify_trna_type(row['trna_id'])
        if classification == trna_type:
            filtered_rows.append(row)
        else:
            excluded_count += 1

    if not filtered_rows:
        print(f"[error] No {trna_type} tRNAs found in dataset")
        return

    print(f"[info] Processing {len(filtered_rows)} {trna_type} tRNAs (excluded {excluded_count} other types)")

    # Convert to DataFrame
    df = pd.DataFrame(filtered_rows).sort_values(["trna_id", "seq_index"]).reset_index(drop=True)

    # Build preferred labels
    pref = build_pref_label(df)

    # Use type-specific label ordering
    if trna_type == 'type1':
        uniq_labels, to_ord = build_global_label_order_type1(pref)
    elif trna_type == 'type2':
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
        print(f"[info] {trna_type.upper()}: Collision validation bypassed due to --allow-collisions flag")
    else:
        validate_no_global_index_collisions(df)

    # Region annotation
    df["region"] = compute_region_column(df)

    # Write output
    cols = [
        "trna_id", "source_file", "seq_index", "sprinzl_index", "sprinzl_label",
        "residue", "sprinzl_ordinal", "sprinzl_continuous", "global_index", "region",
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
        "--allow-collisions", action="store_true",
        help="Allow coordinate generation despite isoacceptor-level position collisions."
    )
    ap.add_argument(
        "--dual-system", action="store_true",
        help="Generate separate coordinate files for Type I and Type II tRNAs."
    )
    ap.add_argument(
        "--type", choices=['type1', 'type2'],
        help="Generate coordinates for specific tRNA type only (overrides --dual-system)."
    )
    args = ap.parse_args()

    paths = glob(os.path.join(args.json_dir, "**", "*.enriched.json"), recursive=True)
    if not paths:
        print(f"[error] No *.enriched.json files found under: {args.json_dir}")
        sys.exit(2)

    # Collect all tRNA data
    all_rows, skipped = [], 0
    for fp in sorted(paths):
        try:
            all_rows.extend(collect_rows_from_json(fp))
        except Exception as e:
            skipped += 1
            print(f"[warn] Skipping {fp} due to error: {e}")

    print(f"[info] JSON files parsed: {len(paths)}  |  skipped: {skipped}")

    # Determine which coordinate systems to generate
    if args.type:
        # Generate coordinates for specific type only
        output_file = args.out_tsv
        generate_coordinates_for_type(all_rows, args.type, output_file, args.allow_collisions)

    elif args.dual_system:
        # Generate separate coordinate files for both types
        base_name = args.out_tsv
        if base_name.endswith('.tsv'):
            base_name = base_name[:-4]

        type1_file = f"{base_name}_type1.tsv"
        type2_file = f"{base_name}_type2.tsv"

        print(f"[info] Generating dual coordinate system:")
        print(f"  Type I (standard): {type1_file}")
        print(f"  Type II (extended): {type2_file}")

        generate_coordinates_for_type(all_rows, 'type1', type1_file, args.allow_collisions)
        generate_coordinates_for_type(all_rows, 'type2', type2_file, args.allow_collisions)

    else:
        # Legacy unified system (original behavior) - for backwards compatibility
        print("[warn] Using legacy unified coordinate system. Consider using --dual-system for better results.")

        # Filter out excluded tRNAs for unified system
        filtered_rows = []
        for row in all_rows:
            classification = classify_trna_type(row['trna_id'])
            if classification in ['type1', 'type2']:
                filtered_rows.append(row)

        df = pd.DataFrame(filtered_rows).sort_values(["trna_id", "seq_index"]).reset_index(drop=True)

        # Use original unified processing (this will likely have collisions)
        pref = build_pref_label(df)
        uniq_labels, to_ord = build_global_label_order(pref)
        df["sprinzl_ordinal"] = pd.to_numeric(pref.map(to_ord), errors="coerce")

        ord_series = pref.map(to_ord)
        cont = []
        for _, sub in df.groupby("trna_id", sort=False):
            cont.append(make_continuous_for_trna(sub, ord_series))
        df["sprinzl_continuous"] = pd.concat(cont).sort_index().astype("float64")

        cont_round = df["sprinzl_continuous"].round(PRECISION)
        uniq_cont = sorted(cont_round.dropna().unique().tolist())
        cont_to_global = {v: i + 1 for i, v in enumerate(uniq_cont)}
        df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

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


if __name__ == "__main__":
    main()
