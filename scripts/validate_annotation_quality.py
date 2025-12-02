#!/usr/bin/env python3
"""
Validate R2DT annotation quality and detect poorly-annotated tRNAs.

This script analyzes R2DT enriched JSON files to identify tRNAs with annotation
issues that may need to be excluded from coordinate generation. Run this after
regenerating R2DT annotations to check for:

1. High percentage of empty sprinzl_labels (>20% threshold)
2. Missing anticodon labels (positions 34-36)
3. Large gaps in annotation (>5 consecutive unlabeled positions)

Usage:
    python scripts/validate_annotation_quality.py outputs/hg38_jsons
    python scripts/validate_annotation_quality.py outputs/sacCer_jsons --threshold 15
    python scripts/validate_annotation_quality.py outputs/hg38_jsons --mito-only

Output:
    Prints a report of tRNAs with potential annotation issues.
    Use this to update EXCLUDED_POORLY_ANNOTATED in trnas_in_space.py.

When to run:
    - After regenerating R2DT annotations (new R2DT models)
    - When adding a new organism
    - When investigating alignment issues
"""

import argparse
import json
import os
import re
from glob import glob
from pathlib import Path
from typing import Dict, List, Tuple


def infer_trna_id_from_filename(path: str) -> str:
    """Extract tRNA ID from filename."""
    base = os.path.basename(path)
    name = base
    for suf in (".enriched.json", ".json"):
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    m = re.match(r"^(.*?)-[A-Z]_[A-Za-z0-9]+$", name)
    trna_id = m.group(1) if m else name
    # Convert T to U in anticodon
    def replace_anticodon_t_with_u(match):
        return match.group(1) + match.group(2).replace("T", "U") + match.group(3)
    trna_id = re.sub(r"(tRNA-[A-Za-z0-9]+-)([ACGTU]{3})(-)", replace_anticodon_t_with_u, trna_id, flags=re.IGNORECASE)
    return trna_id


def is_mitochondrial_trna(trna_id: str) -> bool:
    """Check if a tRNA is mitochondrial."""
    if trna_id is None:
        return False
    trna_id_upper = trna_id.upper()
    return "MITO-TRNA" in trna_id_upper or trna_id_upper.startswith("MITO-")


def analyze_json(filepath: str) -> Dict:
    """
    Analyze a single R2DT enriched JSON file for annotation quality.

    Returns dict with:
        - trna_id: tRNA identifier
        - total_positions: total nucleotide positions
        - empty_labels: count of positions with empty templateNumberingLabel
        - empty_pct: percentage of empty labels
        - max_gap: longest consecutive run of empty labels
        - gap_ranges: list of (start, end) tuples for gaps > 3 positions
        - has_anticodon_labels: whether positions 34-36 have labels
        - anticodon_labels: dict of {position: label} for 34-36
    """
    with open(filepath) as f:
        data = json.load(f)

    mol = data["rnaComplexes"][0]["rnaMolecules"][0]
    seq = mol["sequence"]

    trna_id = infer_trna_id_from_filename(filepath)

    positions = []
    for item in seq:
        if item.get("residueName") in ("5'", "3'"):
            continue
        info = item.get("info", {}) or {}
        label = (info.get("templateNumberingLabel", "") or "").strip()
        idx = info.get("templateResidueIndex", -1)
        positions.append({
            "label": label,
            "index": idx,
            "residue": item.get("residueName", "?"),
        })

    total = len(positions)
    empty_count = sum(1 for p in positions if p["label"] == "")
    empty_pct = 100 * empty_count / total if total > 0 else 0

    # Find gaps (consecutive empty labels)
    gaps = []
    current_gap_start = None
    for i, p in enumerate(positions):
        if p["label"] == "":
            if current_gap_start is None:
                current_gap_start = i
        else:
            if current_gap_start is not None:
                gap_len = i - current_gap_start
                if gap_len > 3:
                    gaps.append((current_gap_start + 1, i))  # 1-indexed
                current_gap_start = None
    # Handle gap at end
    if current_gap_start is not None:
        gap_len = len(positions) - current_gap_start
        if gap_len > 3:
            gaps.append((current_gap_start + 1, len(positions)))

    max_gap = 0
    for start, end in gaps:
        max_gap = max(max_gap, end - start)

    # Check anticodon labels (positions with templateResidueIndex 34, 35, 36)
    anticodon_labels = {}
    for p in positions:
        if p["index"] in (34, 35, 36):
            anticodon_labels[p["index"]] = p["label"]

    has_anticodon = all(
        anticodon_labels.get(pos, "") != ""
        for pos in (34, 35, 36)
    )

    return {
        "trna_id": trna_id,
        "filepath": filepath,
        "total_positions": total,
        "empty_labels": empty_count,
        "empty_pct": empty_pct,
        "max_gap": max_gap,
        "gap_ranges": gaps,
        "has_anticodon_labels": has_anticodon,
        "anticodon_labels": anticodon_labels,
        "is_mito": is_mitochondrial_trna(trna_id),
    }


def validate_directory(json_dir: str, threshold: float = 20.0,
                       mito_only: bool = False, nuclear_only: bool = False) -> List[Dict]:
    """
    Validate all R2DT JSON files in a directory.

    Args:
        json_dir: Directory containing *.enriched.json files
        threshold: Percentage threshold for flagging empty labels (default 20%)
        mito_only: Only check mitochondrial tRNAs
        nuclear_only: Only check nuclear tRNAs

    Returns:
        List of dicts for tRNAs with potential issues
    """
    paths = glob(os.path.join(json_dir, "**", "*.enriched.json"), recursive=True)

    if not paths:
        print(f"No *.enriched.json files found in {json_dir}")
        return []

    issues = []

    for fp in sorted(paths):
        result = analyze_json(fp)

        # Filter by mito/nuclear
        if mito_only and not result["is_mito"]:
            continue
        if nuclear_only and result["is_mito"]:
            continue

        # Check for issues
        has_issue = False
        issue_reasons = []

        if result["empty_pct"] > threshold:
            has_issue = True
            issue_reasons.append(f"{result['empty_pct']:.1f}% empty labels")

        if not result["has_anticodon_labels"]:
            has_issue = True
            issue_reasons.append("missing anticodon labels (34-36)")

        if result["max_gap"] > 5:
            has_issue = True
            issue_reasons.append(f"large gap ({result['max_gap']} positions)")

        if has_issue:
            result["issue_reasons"] = issue_reasons
            issues.append(result)

    return issues


def print_report(issues: List[Dict], json_dir: str):
    """Print a formatted report of annotation issues."""
    print()
    print("=" * 80)
    print("R2DT ANNOTATION QUALITY REPORT")
    print("=" * 80)
    print(f"Directory: {json_dir}")
    print()

    if not issues:
        print("No annotation issues detected!")
        print()
        print("All tRNAs have:")
        print("  - Less than 20% empty labels")
        print("  - Anticodon positions (34-36) labeled")
        print("  - No gaps longer than 5 positions")
        return

    # Separate mito and nuclear
    mito_issues = [i for i in issues if i["is_mito"]]
    nuclear_issues = [i for i in issues if not i["is_mito"]]

    if nuclear_issues:
        print(f"NUCLEAR tRNAs WITH ISSUES ({len(nuclear_issues)}):")
        print("-" * 80)
        for issue in sorted(nuclear_issues, key=lambda x: -x["empty_pct"]):
            print(f"\n  {issue['trna_id']}")
            print(f"    Empty labels: {issue['empty_labels']}/{issue['total_positions']} ({issue['empty_pct']:.1f}%)")
            if issue["gap_ranges"]:
                print(f"    Gaps: {issue['gap_ranges']}")
            if not issue["has_anticodon_labels"]:
                print(f"    Anticodon labels: {issue['anticodon_labels']} (INCOMPLETE)")
            print(f"    Issues: {', '.join(issue['issue_reasons'])}")
        print()

    if mito_issues:
        print(f"MITOCHONDRIAL tRNAs WITH ISSUES ({len(mito_issues)}):")
        print("-" * 80)
        for issue in sorted(mito_issues, key=lambda x: -x["empty_pct"]):
            print(f"\n  {issue['trna_id']}")
            print(f"    Empty labels: {issue['empty_labels']}/{issue['total_positions']} ({issue['empty_pct']:.1f}%)")
            if issue["gap_ranges"]:
                print(f"    Gaps: {issue['gap_ranges']}")
            if not issue["has_anticodon_labels"]:
                print(f"    Anticodon labels: {issue['anticodon_labels']} (INCOMPLETE)")
            print(f"    Issues: {', '.join(issue['issue_reasons'])}")
        print()

    # Print exclusion list format
    print("=" * 80)
    print("SUGGESTED ADDITIONS TO EXCLUDED_POORLY_ANNOTATED:")
    print("=" * 80)
    print()
    print("Add to scripts/trnas_in_space.py EXCLUDED_POORLY_ANNOTATED set:")
    print()
    for issue in issues:
        reasons = ", ".join(issue["issue_reasons"])
        print(f'    "{issue["trna_id"]}",  # {reasons}')
    print()
    print("Review each tRNA before adding to exclusion list.")
    print("See docs/EXCLUDED_TRNAS.md for documentation requirements.")


def main():
    parser = argparse.ArgumentParser(
        description="Validate R2DT annotation quality and detect poorly-annotated tRNAs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Check all tRNAs in a directory
    python scripts/validate_annotation_quality.py outputs/hg38_jsons

    # Check only mitochondrial tRNAs
    python scripts/validate_annotation_quality.py outputs/sacCer_jsons --mito-only

    # Use a stricter threshold (15% instead of 20%)
    python scripts/validate_annotation_quality.py outputs/hg38_jsons --threshold 15

When to run this script:
    - After regenerating R2DT annotations (new R2DT version/models)
    - When adding a new organism to the pipeline
    - When investigating coordinate alignment issues
    - Before releasing updated coordinate files
        """
    )
    parser.add_argument("json_dir", help="Directory containing R2DT *.enriched.json files")
    parser.add_argument("--threshold", type=float, default=20.0,
                        help="Percentage threshold for flagging empty labels (default: 20)")
    parser.add_argument("--mito-only", action="store_true",
                        help="Only check mitochondrial tRNAs")
    parser.add_argument("--nuclear-only", action="store_true",
                        help="Only check nuclear tRNAs")

    args = parser.parse_args()

    if args.mito_only and args.nuclear_only:
        parser.error("Cannot specify both --mito-only and --nuclear-only")

    issues = validate_directory(
        args.json_dir,
        threshold=args.threshold,
        mito_only=args.mito_only,
        nuclear_only=args.nuclear_only
    )

    print_report(issues, args.json_dir)


if __name__ == "__main__":
    main()
