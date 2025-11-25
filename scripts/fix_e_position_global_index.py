#!/usr/bin/env python3
"""
fix_e_position_global_index.py

This script fixes the global_index assignment for extended variable arm positions
("e" positions) in existing coordinate files. It recalculates global_index values
using the corrected sort_key function.

Usage:
    python scripts/fix_e_position_global_index.py <input.tsv> <output.tsv>
"""

import sys
import os
import pandas as pd

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(__file__))
import trnas_in_space  # noqa: E402

PRECISION = 6  # Must match PRECISION in trnas_in_space.py


def fix_global_index(input_tsv: str, output_tsv: str):
    """
    Reprocess a coordinate TSV file to fix global_index values.

    This reads the existing file, recalculates global_index using the fixed
    sort_key function, and writes the corrected output.
    """
    print(f"Reading {input_tsv}...")
    df = pd.read_csv(input_tsv, sep='\t')

    print(f"  Rows: {len(df)}")
    print(f"  tRNAs: {df['trna_id'].nunique()}")

    # Build preferred label (same logic as original script)
    lbl = df["sprinzl_label"].astype("string").fillna("").str.strip()
    num = pd.to_numeric(df["sprinzl_index"], errors="coerce")
    num_ok = num.where((num >= 1) & (num <= 76))
    num_ok_str = num_ok.astype("Int64").astype(str).replace({"<NA>": ""})
    pref = num_ok_str.mask(num_ok_str.eq(""), other=lbl)

    # Build global label order using FIXED sort_key
    uniq = sorted({p for p in pref if p not in ("", "nan")}, key=trnas_in_space.sort_key)
    to_ord = {u: i + 1 for i, u in enumerate(uniq)}

    print(f"  Unique labels: {len(uniq)}")
    print("  Sample labels around position 47-48:")
    for i, label in enumerate(uniq):
        if label in ["44", "45", "46", "47", "48", "49"] or label.startswith("e"):
            print(f"    {i+1:3d}: {label}")

    # Assign ordinals
    df["sprinzl_ordinal"] = pd.to_numeric(pref.map(to_ord), errors="coerce")

    # Recalculate continuous coordinates (same logic as original)
    ord_series = pref.map(to_ord)
    cont = []
    for _, sub in df.groupby("trna_id", sort=False):
        cont.append(trnas_in_space.make_continuous_for_trna(sub, ord_series))
    df["sprinzl_continuous"] = pd.concat(cont).sort_index().astype("float64")

    # Recalculate global_index
    cont_round = df["sprinzl_continuous"].round(PRECISION)
    uniq_cont = sorted(cont_round.dropna().unique().tolist())
    cont_to_global = {v: i + 1 for i, v in enumerate(uniq_cont)}
    df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

    print(f"  New global positions (K): {len(uniq_cont)}")

    # Check for "e" positions and report their new global_index values
    e_positions = df[df["sprinzl_label"].str.startswith("e", na=False)]
    if len(e_positions) > 0:
        print(f"\n  Found {len(e_positions)} 'e' position rows:")
        e_summary = e_positions.groupby("sprinzl_label")["global_index"].first().sort_values()
        for label, gidx in e_summary.items():
            print(f"    {label}: global_index {gidx}")

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

    df.to_csv(output_tsv, sep="\t", index=False, columns=cols)
    print(f"\n✓ Wrote {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    fix_global_index(input_file, output_file)
