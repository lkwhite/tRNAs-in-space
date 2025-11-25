#!/usr/bin/env python3
"""
Bulk processing script for generating tRNA global coordinates for multiple organisms.

This script automates the workflow:
1. Read organism configuration from config/organisms.yaml
2. For each organism with status='pending':
   - Verify FASTA file exists
   - Run R2DT annotation (Docker)
   - Generate global coordinates
   - Validate output
   - Update status

Usage:
    python scripts/process_organisms.py [--organisms ORGANISM_IDS] [--skip-r2dt] [--dry-run]

Examples:
    # Process all pending organisms
    python scripts/process_organisms.py

    # Process specific organisms
    python scripts/process_organisms.py --organisms mm39,dm6,ce11

    # Skip R2DT and only generate coordinates (if JSONs already exist)
    python scripts/process_organisms.py --skip-r2dt

    # Dry run to see what would be processed
    python scripts/process_organisms.py --dry-run
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

# Ensure we're in the project root
PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)


def load_config():
    """Load organism configuration from YAML file."""
    config_path = PROJECT_ROOT / "config" / "organisms.yaml"
    if not config_path.exists():
        print(f"Error: Configuration file not found at {config_path}")
        sys.exit(1)

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    return config


def get_organisms_to_process(config, organism_ids=None):
    """
    Get list of organisms to process.

    Args:
        config: Loaded YAML configuration
        organism_ids: Optional list of specific organism IDs to process

    Returns:
        List of organism dictionaries to process
    """
    all_organisms = config["organisms"]

    if organism_ids:
        # Process specific organisms
        organisms = [org for org in all_organisms if org["organism_id"] in organism_ids]
        if not organisms:
            print(f"Error: No organisms found with IDs: {organism_ids}")
            sys.exit(1)
    else:
        # Process all pending organisms (priority 1)
        organisms = [
            org
            for org in all_organisms
            if org.get("status") == "pending" and org.get("priority") == 1
        ]

    return organisms


def check_fasta_exists(organism):
    """Check if FASTA file exists for organism."""
    gtrnadb_id = organism["gtrnadb_id"]

    # Try common naming patterns
    possible_names = [
        f"fastas/{gtrnadb_id}-tRNAs.fa",
        f"fastas/{gtrnadb_id}-mito-and-nuclear-tRNAs.fa",
        f"fastas/{organism['organism_id']}-tRNAs.fa",
    ]

    for fasta_path in possible_names:
        if Path(fasta_path).exists():
            return fasta_path

    return None


def run_r2dt(organism, fasta_path):
    """
    Run R2DT annotation on FASTA file.

    Args:
        organism: Organism dictionary from config
        fasta_path: Path to input FASTA file

    Returns:
        Path to output JSON directory
    """
    organism_id = organism["organism_id"]
    json_output_dir = PROJECT_ROOT / "outputs" / "jsons" / organism_id
    json_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print(f"Running R2DT annotation for {organism['name']}")
    print(f"{'='*80}")
    print(f"Input FASTA: {fasta_path}")
    print(f"Output directory: {json_output_dir}")

    # Build Docker command
    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{PROJECT_ROOT}:/data",
        "rnacentral/r2dt",
        "r2dt.py",
        "gtrnadb",
        "draw",
        f"/data/{fasta_path}",
        f"/data/{json_output_dir.relative_to(PROJECT_ROOT)}",
    ]

    print(f"\nCommand: {' '.join(cmd)}")
    print("\nThis may take several minutes...\n")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("R2DT completed successfully!")
        return json_output_dir
    except subprocess.CalledProcessError as e:
        print(f"Error running R2DT: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return None


def generate_coordinates(organism, json_dir):
    """
    Generate global coordinates from R2DT JSON output.

    Args:
        organism: Organism dictionary from config
        json_dir: Path to R2DT JSON output directory

    Returns:
        Path to output TSV file
    """
    organism_id = organism["organism_id"]
    output_tsv = PROJECT_ROOT / "outputs" / f"{organism_id}_global_coords.tsv"

    print(f"\n{'='*80}")
    print(f"Generating global coordinates for {organism['name']}")
    print(f"{'='*80}")
    print(f"Input JSON directory: {json_dir}")
    print(f"Output TSV: {output_tsv}")

    # Build command
    cmd = ["python", "scripts/trnas_in_space.py", str(json_dir), str(output_tsv)]

    print(f"\nCommand: {' '.join(cmd)}\n")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        print("Coordinate generation completed successfully!")
        return output_tsv
    except subprocess.CalledProcessError as e:
        print(f"Error generating coordinates: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return None


def validate_output(output_tsv):
    """
    Perform basic validation on output TSV.

    Args:
        output_tsv: Path to output TSV file

    Returns:
        Dictionary with validation results
    """
    print(f"\n{'='*80}")
    print(f"Validating output: {output_tsv}")
    print(f"{'='*80}")

    if not output_tsv.exists():
        print("Error: Output file does not exist!")
        return {"success": False, "error": "File not found"}

    # Read and validate
    import pandas as pd

    try:
        df = pd.read_csv(output_tsv, sep="\t")

        # Expected columns
        expected_cols = [
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

        # Check columns
        missing_cols = set(expected_cols) - set(df.columns)
        if missing_cols:
            return {"success": False, "error": f"Missing columns: {missing_cols}"}

        # Get statistics
        n_trnas = df["trna_id"].nunique()
        n_nucleotides = len(df)
        n_global_positions = df["global_index"].max()

        stats = {
            "success": True,
            "n_trnas": n_trnas,
            "n_nucleotides": n_nucleotides,
            "n_global_positions": int(n_global_positions),
            "file_size_kb": round(output_tsv.stat().st_size / 1024, 2),
        }

        print("✓ Valid TSV file")
        print(f"✓ tRNAs: {n_trnas}")
        print(f"✓ Nucleotides: {n_nucleotides}")
        print(f"✓ Global positions: {n_global_positions}")
        print(f"✓ File size: {stats['file_size_kb']} KB")

        return stats

    except Exception as e:
        print(f"Error validating output: {e}")
        return {"success": False, "error": str(e)}


def process_organism(organism, skip_r2dt=False, dry_run=False):
    """
    Process a single organism through the full pipeline.

    Args:
        organism: Organism dictionary from config
        skip_r2dt: Skip R2DT and use existing JSONs
        dry_run: Don't actually run commands

    Returns:
        Dictionary with processing results
    """
    organism_id = organism["organism_id"]
    organism_name = organism["name"]

    print(f"\n\n{'#'*80}")
    print(f"# Processing: {organism_name} ({organism_id})")
    print(f"{'#'*80}")

    result = {
        "organism_id": organism_id,
        "organism_name": organism_name,
        "success": False,
        "steps_completed": [],
    }

    # Check if FASTA exists
    fasta_path = check_fasta_exists(organism)
    if not fasta_path:
        print(f"\n❌ ERROR: FASTA file not found for {organism_name}")
        print(f"   Expected in fastas/ with name pattern: {organism['gtrnadb_id']}-tRNAs.fa")
        print("   Please download from GtRNAdb: http://gtrnadb.ucsc.edu/")
        result["error"] = "FASTA file not found"
        return result

    print(f"✓ FASTA file found: {fasta_path}")
    result["fasta_path"] = fasta_path

    if dry_run:
        print("\n[DRY RUN] Would process this organism")
        result["success"] = True
        return result

    # Run R2DT
    json_dir = None
    if not skip_r2dt:
        json_dir = run_r2dt(organism, fasta_path)
        if not json_dir:
            result["error"] = "R2DT annotation failed"
            return result
        result["steps_completed"].append("r2dt")
    else:
        # Use existing JSON directory
        json_dir = PROJECT_ROOT / "outputs" / "jsons" / organism_id
        if not json_dir.exists():
            print(f"❌ ERROR: JSON directory not found: {json_dir}")
            result["error"] = "JSON directory not found (use without --skip-r2dt)"
            return result
        print(f"✓ Using existing JSON directory: {json_dir}")

    # Generate coordinates
    output_tsv = generate_coordinates(organism, json_dir)
    if not output_tsv:
        result["error"] = "Coordinate generation failed"
        return result
    result["steps_completed"].append("coordinates")
    result["output_tsv"] = str(output_tsv)

    # Validate output
    validation = validate_output(output_tsv)
    if not validation["success"]:
        result["error"] = f"Validation failed: {validation.get('error')}"
        return result
    result["steps_completed"].append("validation")
    result["validation"] = validation

    # Success!
    result["success"] = True
    print(f"\n✅ Successfully processed {organism_name}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Bulk process organisms for tRNA global coordinates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--organisms",
        help="Comma-separated list of organism IDs to process (default: all pending)",
        default=None,
    )
    parser.add_argument(
        "--skip-r2dt", action="store_true", help="Skip R2DT annotation and use existing JSON files"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be processed without actually running",
    )

    args = parser.parse_args()

    # Load configuration
    config = load_config()

    # Get organisms to process
    organism_ids = args.organisms.split(",") if args.organisms else None
    organisms = get_organisms_to_process(config, organism_ids)

    if not organisms:
        print("No organisms to process!")
        sys.exit(0)

    print(f"\n{'='*80}")
    print("tRNAs in Space - Bulk Organism Processing")
    print(f"{'='*80}")
    print(f"Mode: {'DRY RUN' if args.dry_run else 'LIVE'}")
    print(f"Skip R2DT: {args.skip_r2dt}")
    print(f"Organisms to process: {len(organisms)}")
    for org in organisms:
        print(f"  - {org['name']} ({org['organism_id']})")
    print(f"{'='*80}\n")

    if args.dry_run:
        print("This is a DRY RUN - no actual processing will occur\n")

    # Process each organism
    results = []
    for organism in organisms:
        result = process_organism(organism, skip_r2dt=args.skip_r2dt, dry_run=args.dry_run)
        results.append(result)

    # Summary
    print(f"\n\n{'='*80}")
    print("PROCESSING SUMMARY")
    print(f"{'='*80}")

    successful = [r for r in results if r["success"]]
    failed = [r for r in results if not r["success"]]

    print(f"\nTotal organisms: {len(results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")

    if successful:
        print("\n✅ Successfully processed:")
        for r in successful:
            print(f"   - {r['organism_name']}")
            if not args.dry_run and "validation" in r:
                v = r["validation"]
                print(
                    f"     {v['n_trnas']} tRNAs, {v['n_nucleotides']} nucleotides, {v['file_size_kb']} KB"
                )

    if failed:
        print("\n❌ Failed:")
        for r in failed:
            print(f"   - {r['organism_name']}: {r.get('error', 'Unknown error')}")

    # Save results to JSON
    if not args.dry_run:
        results_file = PROJECT_ROOT / "outputs" / "processing_results.json"
        with open(results_file, "w") as f:
            json.dump({"timestamp": datetime.now().isoformat(), "results": results}, f, indent=2)
        print(f"\nResults saved to: {results_file}")

    print(f"\n{'='*80}\n")

    # Exit with error code if any failed
    sys.exit(0 if len(failed) == 0 else 1)


if __name__ == "__main__":
    main()
