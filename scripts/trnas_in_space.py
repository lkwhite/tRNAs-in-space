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


# --------------------- phase 1: JSON -> rows ----------------------


def collect_rows_from_json(fp: str):
    with open(fp, "r") as f:
        J = json.load(f)
    mol = J["rnaComplexes"][0]["rnaMolecules"][0]
    seq = mol["sequence"]

    trna_id = infer_trna_id_from_filename(fp)
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
    """Order labels: 20 < 20A < 20B; allow dotted like 9.1; e positions between 47 and 48; unknowns at end."""
    if lbl is None or lbl == "" or lbl == "nan":
        return (10**9, 2, "")
    s = str(lbl)
    # Handle "e" positions (extended variable arm): e1, e2, e12, etc.
    # Place them between position 47 and 48
    m = re.fullmatch(r"e(\d+)", s)
    if m:
        e_num = int(m.group(1))
        # Return (47, 1, e_num) to place between 47 (base=47, tier=0) and 48 (base=48, tier=0)
        return (47, 1, e_num)
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1))
        suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, int(m.group(2)))
    return (10**9 - 1, 2, s)


def build_pref_label(df: pd.DataFrame) -> pd.Series:
    """
    Prefer sprinzl_index (numeric) for standard positions [1..76] to ensure consistency.
    Use sprinzl_label only for insertions (e.g., "20a", "47a") where index is invalid.
    Returns strings with '' for unknown.
    """
    lbl = df["sprinzl_label"].astype("string").fillna("").str.strip()
    num = pd.to_numeric(df["sprinzl_index"], errors="coerce")
    num_ok = num.where((num >= 1) & (num <= 76))
    num_ok_str = num_ok.astype("Int64").astype(str).replace({"<NA>": ""})
    # Use numeric index when valid (1-76), otherwise fall back to label for insertions
    pref = num_ok_str.mask(num_ok_str.eq(""), other=lbl)
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


def main():
    ap = argparse.ArgumentParser(
        description="Build global equal-spaced tRNA coordinates + regions from R2DT JSON."
    )
    ap.add_argument(
        "json_dir", help="Directory with R2DT *.enriched.json files (searched recursively)."
    )
    ap.add_argument("out_tsv", help="Output TSV path.")
    args = ap.parse_args()

    paths = glob(os.path.join(args.json_dir, "**", "*.enriched.json"), recursive=True)
    if not paths:
        print(f"[error] No *.enriched.json files found under: {args.json_dir}")
        sys.exit(2)

    all_rows, skipped = [], 0
    for fp in sorted(paths):
        try:
            all_rows.extend(collect_rows_from_json(fp))
        except Exception as e:
            skipped += 1
            print(f"[warn] Skipping {fp} due to error: {e}")

    df = pd.DataFrame(all_rows).sort_values(["trna_id", "seq_index"]).reset_index(drop=True)

    # Preferred text label, else numeric 1..76
    pref = build_pref_label(df)

    # Global label order and ordinal
    uniq_labels, to_ord = build_global_label_order(pref)
    df["sprinzl_ordinal"] = pd.to_numeric(pref.map(to_ord), errors="coerce")

    # Continuous coordinate per tRNA
    ord_series = pref.map(to_ord)
    cont = []
    for _, sub in df.groupby("trna_id", sort=False):
        cont.append(make_continuous_for_trna(sub, ord_series))
    df["sprinzl_continuous"] = pd.concat(cont).sort_index().astype("float64")

    # Global equal-spaced index (rounded internally to PRECISION)
    cont_round = df["sprinzl_continuous"].round(PRECISION)
    uniq_cont = sorted(cont_round.dropna().unique().tolist())
    cont_to_global = {v: i + 1 for i, v in enumerate(uniq_cont)}
    df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

    # Region annotation from Sprinzl
    df["region"] = compute_region_column(df)

    # Write out
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

    # Stats
    print(f"[ok] Wrote {args.out_tsv}")
    print(f"  JSON files parsed: {len(paths)}  |  skipped: {skipped}")
    print(f"  Rows: {len(df)}  |  tRNAs: {df['trna_id'].nunique()}")
    print(f"  Unique labeled bins: {len(uniq_labels)}")
    print(f"  Unique global positions (K): {len(uniq_cont)}  |  rounding={PRECISION} d.p.")
    missing = (df["sprinzl_index"] < 1).sum()
    if missing:
        print(f"  Note: {missing} rows still have sprinzl_index == -1 (unresolvable from JSON).")


if __name__ == "__main__":
    main()
