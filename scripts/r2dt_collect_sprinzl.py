#!/usr/bin/env python3
"""
r2dt_collect_sprinzl.py

Usage:
  python r2dt_collect_sprinzl.py /path/to/r2dt_output_dir [out.tsv]

What it does:
  - Recursively finds all *.enriched.json files under the given directory.
  - For each file, loads the R2DT drawing JSON and extracts per-nucleotide rows:
        trna_id, source_file, seq_index, sprinzl_index, sprinzl_label, residue
  - Fills missing Sprinzl indices (-1) by monotonic interpolation along sequence:
        exact (from R2DT), interp (internal gap), lead (leading run), trail (trailing run).
    Unresolvable cases remain -1.
  - Writes one combined TSV (default: r2dt_sprinzl_all.tsv).
"""

import os, sys, json, csv, re
from glob import glob

def infer_trna_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    name = base
    for suf in (".enriched.json", ".json"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    m = re.match(r"^(.*?)-[A-Z]_[A-Za-z0-9]+$", name)  # strip -B_His etc.
    return m.group(1) if m else name

def process_file(fp: str):
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
        ridx = s["residueIndex"]                          # sequence index as in JSON
        info = s.get("info", {}) or {}
        sprinzl_idx = info.get("templateResidueIndex", None)
        sprinzl_idx = sprinzl_idx if isinstance(sprinzl_idx, int) else -1
        sprinzl_lbl = (info.get("templateNumberingLabel", "") or "").strip()

        rows.append({
            "trna_id": trna_id,
            "source_file": os.path.basename(fp),
            "seq_index": ridx,
            "sprinzl_index": sprinzl_idx,                 # will be overwritten after fill
            "sprinzl_label": sprinzl_lbl,
            "residue": rname,
        })

    # ---- fill missing Sprinzl indices by monotonic inference along sequence ----
    rows.sort(key=lambda r: r["seq_index"])
    n = len(rows)
    vals = [r["sprinzl_index"] for r in rows]

    # forward pass (carry upward)
    fwd = [None]*n
    last = None
    for i in range(n):
        v = vals[i]
        if isinstance(v, int) and v >= 1:
            fwd[i] = v
            last = v
        elif last is not None:
            fwd[i] = last + 1
            last += 1

    # backward pass (carry downward)
    bwd = [None]*n
    nxt = None
    for i in range(n-1, -1, -1):
        v = vals[i]
        if isinstance(v, int) and v >= 1:
            bwd[i] = v
            nxt = v
        elif nxt is not None:
            bwd[i] = nxt - 1
            nxt -= 1

    # decide final sprinzl_index (overwrite)
    for i, r in enumerate(rows):
        v = vals[i]
        if v >= 1:
            r["sprinzl_index"] = v
            continue
        # internal gap bounded and agreeing
        if fwd[i] is not None and bwd[i] is not None and fwd[i] == bwd[i] and 1 <= fwd[i] <= 76:
            r["sprinzl_index"] = fwd[i]
        # only leading side known
        elif fwd[i] is not None and 1 <= fwd[i] <= 76 and bwd[i] is None:
            r["sprinzl_index"] = fwd[i]
        # only trailing side known
        elif bwd[i] is not None and 1 <= bwd[i] <= 76 and fwd[i] is None:
            r["sprinzl_index"] = bwd[i]
        else:
            r["sprinzl_index"] = -1  # leave as missing

    return rows

def main():
    if len(sys.argv) < 2:
        print("Usage: python r2dt_collect_sprinzl.py /path/to/r2dt_output_dir [out.tsv]")
        sys.exit(1)

    in_dir = sys.argv[1]
    out_tsv = sys.argv[2] if len(sys.argv) >= 3 else "r2dt_sprinzl_all.tsv"

    paths = glob(os.path.join(in_dir, "**", "*.enriched.json"), recursive=True)
    if not paths:
        print(f"[error] No *.enriched.json files found under: {in_dir}")
        sys.exit(2)

    all_rows = []
    for fp in sorted(paths):
        try:
            all_rows.extend(process_file(fp))
        except Exception as e:
            print(f"[warn] Skipping {fp} due to error: {e}")

    cols = ["trna_id", "source_file", "seq_index", "sprinzl_index", "sprinzl_label", "residue"]
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)

    print(f"[ok] Wrote {len(all_rows)} rows from {len(paths)} files to {out_tsv}")
    print("[tip] Filter canonical positions with: 1 <= sprinzl_index <= 76")

if __name__ == "__main__":
    main()

