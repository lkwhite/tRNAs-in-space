#!/usr/bin/env python3
"""
make_sprinzl_continuous.py

Usage:
  python make_sprinzl_continuous.py input.tsv [out.tsv] [precision]

Input TSV must have columns from r2dt_collect_sprinzl.py:
  trna_id  source_file  seq_index  sprinzl_index  sprinzl_label  residue

Outputs the same table plus:
  sprinzl_ordinal      # global integer bin based on label text (1<2<...<20<20A<20B<21<...<76)
  sprinzl_continuous   # per-tRNA continuous coordinate (fills unlabeled runs with fractions)
  global_index         # equal-spaced 1..K built from sorted unique sprinzl_continuous across all tRNAs

'precision' (default 12) controls rounding before uniquing to avoid FP jitter.
"""

import sys, re
import pandas as pd
import numpy as np

def sort_key(lbl: str):
    """Order labels: 20 < 20A < 20B; also handle dotted like 9.1."""
    if lbl is None or lbl == "" or lbl == "nan":
        return (10**9, 2, "")
    s = str(lbl)
    m = re.fullmatch(r"(\d+)([A-Za-z]+)?", s)
    if m:
        base = int(m.group(1)); suf = (m.group(2) or "").upper()
        return (base, 1 if suf else 0, suf)
    m = re.fullmatch(r"(\d+)\.(\d+)", s)
    if m:
        return (int(m.group(1)), 1, int(m.group(2)))
    return (10**9-1, 2, s)

def build_pref_label(df: pd.DataFrame) -> pd.Series:
    """
    Prefer sprinzl_label (text). If missing/empty, fall back to sprinzl_index
    when it is an integer in [1..76]. Returns a string Series with '' where unknown.
    """
    lbl = df["sprinzl_label"].astype("string").fillna("").str.strip()
    # numeric fallback in [1..76]
    num = pd.to_numeric(df["sprinzl_index"], errors="coerce")
    num_ok = num.where((num >= 1) & (num <= 76))
    # Convert to string, but keep '' for NaNs
    num_ok_str = num_ok.astype("Int64").astype(str).replace({"<NA>": ""})
    # If lbl is empty, take num_ok_str; otherwise keep lbl
    pref = lbl.mask(lbl.eq(""), other=num_ok_str)
    return pref

def build_global_label_order(pref: pd.Series):
    uniq = sorted({p for p in pref if p not in ("", "nan")}, key=sort_key)
    to_ord = {u: i + 1 for i, u in enumerate(uniq)}  # 1..K
    return uniq, to_ord

def make_continuous_for_trna(sub: pd.DataFrame, ord_series: pd.Series) -> pd.Series:
    """
    Build continuous coordinates for one tRNA:
      - exact labels -> integer ordinals
      - unlabeled runs -> fractions between neighboring integer ordinals
      - leading/trailing runs -> fractions just below/above edge bins
    """
    sub = sub.sort_values("seq_index").copy()
    ords = ord_series.loc[sub.index].astype("Float64").to_numpy()

    n = len(sub)
    vals = [np.nan]*n

    i = 0
    while i < n:
        if not pd.isna(ords[i]):
            vals[i] = float(int(ords[i]))
            i += 1
            continue
        # start unlabeled run
        j = i
        while j < n and pd.isna(ords[j]):
            j += 1
        k = j - i  # run length
        left = int(ords[i-1]) if i-1 >= 0 and not pd.isna(ords[i-1]) else None
        right = int(ords[j])   if j   <  n and not pd.isna(ords[j])   else None

        if left is not None and right is not None and right >= left + 1:
            # internal gap between bins [left, right]
            for t in range(k):
                vals[i+t] = left + (t+1)/(k+1)
        elif left is None and right is not None:
            # leading run before first labeled bin
            for t in range(k):
                vals[i+t] = right - (k - t)/(k+1)
        elif left is not None and right is None:
            # trailing run after last labeled bin
            for t in range(k):
                vals[i+t] = left + (t+1)/(k+1)
        else:
            # all unlabeled (rare)
            for t in range(k):
                vals[i+t] = np.nan
        i = j

    return pd.Series(vals, index=sub.index, dtype="float64")

def main():
    if len(sys.argv) < 2:
        print("Usage: python make_sprinzl_continuous.py input.tsv [out.tsv] [precision]")
        sys.exit(1)

    in_tsv = sys.argv[1]
    out_tsv = sys.argv[2] if len(sys.argv) >= 3 else in_tsv.replace(".tsv", "_with_continuous.tsv")
    precision = int(sys.argv[3]) if len(sys.argv) >= 4 else 12

    df = pd.read_csv(in_tsv, sep="\t")

    # 1) Preferred label per row (text first, else numeric 1..76)
    pref = build_pref_label(df)

    # 2) Global label order and numeric ordinal for labeled rows
    uniq_labels, to_ord = build_global_label_order(pref)
    sprinzl_ordinal = pref.map(to_ord)  # Ints for labeled rows, NaN otherwise
    df["sprinzl_ordinal"] = pd.to_numeric(sprinzl_ordinal, errors="coerce")

    # 3) sprinzl_continuous per tRNA using the global label order
    df = df.sort_values(["trna_id","seq_index"]).reset_index(drop=True)
    ord_series = pref.map(to_ord)  # numeric or NaN
    cont_chunks = []
    for tid, sub in df.groupby("trna_id", sort=False):
        cont_chunks.append(make_continuous_for_trna(sub, ord_series))
    df["sprinzl_continuous"] = pd.concat(cont_chunks).sort_index().astype("float64")

    # 4) global_index — equal-spaced, shared across all tRNAs
    cont_round = df["sprinzl_continuous"].round(precision)
    uniq_cont = sorted(cont_round.dropna().unique().tolist())
    cont_to_global = {v: i+1 for i, v in enumerate(uniq_cont)}  # 1..K
    df["global_index"] = cont_round.map(cont_to_global).astype("Int64")

    # Write
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[ok] wrote {out_tsv}")
    print(f"Unique global positions (K): {len(uniq_cont)}  |  rounding precision={precision}")

if __name__ == "__main__":
    main()

