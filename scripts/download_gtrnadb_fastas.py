#!/usr/bin/env python3
"""
Download FASTA files from GtRNAdb for configured organisms.

This script attempts to download tRNA FASTA files from GtRNAdb based on
the organism configuration in config/organisms.yaml.

Usage:
    python scripts/download_gtrnadb_fastas.py [--organisms ORGANISM_IDS] [--manual]

Examples:
    # Download all pending organisms
    python scripts/download_gtrnadb_fastas.py

    # Download specific organisms
    python scripts/download_gtrnadb_fastas.py --organisms mm39,dm6

    # Just print download URLs for manual download
    python scripts/download_gtrnadb_fastas.py --manual
"""

import argparse
import os
import sys
import urllib.request
import urllib.error
from pathlib import Path
import yaml
import time

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


def construct_gtrnadb_url(organism):
    """
    Construct GtRNAdb download URL for an organism.

    GtRNAdb URLs typically follow patterns like:
    - http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa
    - http://gtrnadb.ucsc.edu/genomes/bacteria/Ecoli_K12_MG1655/ecoliK12MG1655-tRNAs.fa

    Args:
        organism: Organism dictionary from config

    Returns:
        Tuple of (url, filename) or None if URL cannot be constructed
    """
    kingdom = organism["kingdom"].lower()
    gtrnadb_id = organism["gtrnadb_id"]
    organism_id = organism["organism_id"]

    # Determine kingdom path
    kingdom_map = {
        "eukaryota": "eukaryota",
        "bacteria": "bacteria",
        "archaea": "archaea",
        "plantae": "eukaryota",  # Plants are in eukaryota
    }
    kingdom_path = kingdom_map.get(kingdom, kingdom)

    # Try different URL patterns
    base_url = "http://gtrnadb.ucsc.edu/genomes"

    # Pattern 1: {gtrnadb_id}-mature-tRNAs.fa (common for eukaryotes)
    # Pattern 2: {gtrnadb_id}-tRNAs.fa (common for bacteria)
    # Pattern 3: {organism_id}-tRNAs.fa

    # Construct organism name for URL path
    # This is tricky as GtRNAdb uses different naming conventions
    org_name_patterns = []

    if kingdom == "Eukaryota" or kingdom == "Plantae":
        # Common eukaryote patterns
        name_parts = organism["name"].split()
        if len(name_parts) >= 2:
            # e.g., "Homo sapiens" -> "Hsapi38" or "Hsapi"
            genus_initial = name_parts[0][0].upper()
            species_abbrev = name_parts[1][:3].lower()  # or [:4]
            org_name_patterns.append(
                f"{genus_initial}{species_abbrev}{organism_id.replace(organism_id[:len(genus_initial + species_abbrev)], '')}"
            )
            org_name_patterns.append(f"{genus_initial}{species_abbrev}")

    # Add the gtrnadb_id as a pattern
    org_name_patterns.append(gtrnadb_id)

    # Filename patterns
    filename_patterns = [
        f"{gtrnadb_id}-mature-tRNAs.fa",
        f"{gtrnadb_id}-tRNAs.fa",
        f"{organism_id}-mature-tRNAs.fa",
        f"{organism_id}-tRNAs.fa",
    ]

    # Add mito/nuclear variant for eukaryotes
    if "mitochondria" in organism.get("compartments", []):
        filename_patterns.extend(
            [
                f"{gtrnadb_id}-mito-and-nuclear-tRNAs.fa",
                f"{organism_id}-mito-and-nuclear-tRNAs.fa",
            ]
        )

    # Generate all possible URLs
    urls = []
    for org_pattern in org_name_patterns:
        for filename in filename_patterns:
            url = f"{base_url}/{kingdom_path}/{org_pattern}/{filename}"
            urls.append((url, filename))

    return urls


def download_file(url, output_path, timeout=30):
    """
    Download a file from URL to output path.

    Args:
        url: URL to download from
        output_path: Path to save file
        timeout: Download timeout in seconds

    Returns:
        True if successful, False otherwise
    """
    try:
        print(f"  Trying: {url}")
        with urllib.request.urlopen(url, timeout=timeout) as response:
            content = response.read()

            # Verify it's a FASTA file
            if not content.startswith(b">"):
                print("  ❌ Not a FASTA file (doesn't start with '>')")
                return False

            # Write to file
            with open(output_path, "wb") as f:
                f.write(content)

            # Get file size
            file_size = len(content)
            print(f"  ✓ Downloaded {file_size} bytes")
            return True

    except urllib.error.HTTPError as e:
        print(f"  ❌ HTTP Error {e.code}: {e.reason}")
        return False
    except urllib.error.URLError as e:
        print(f"  ❌ URL Error: {e.reason}")
        return False
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return False


def download_organism_fasta(organism, manual_only=False):
    """
    Download FASTA file for an organism.

    Args:
        organism: Organism dictionary from config
        manual_only: If True, only print URLs without downloading

    Returns:
        Dictionary with download results
    """
    organism_id = organism["organism_id"]
    organism_name = organism["name"]

    print(f"\n{'='*80}")
    print(f"Organism: {organism_name} ({organism_id})")
    print(f"{'='*80}")

    # Check if file already exists
    fastas_dir = PROJECT_ROOT / "fastas"
    fastas_dir.mkdir(exist_ok=True)

    existing_files = list(fastas_dir.glob(f"*{organism_id}*.fa")) + list(
        fastas_dir.glob(f"*{organism['gtrnadb_id']}*.fa")
    )

    if existing_files:
        print(f"✓ FASTA file already exists: {existing_files[0]}")
        return {"success": True, "file": str(existing_files[0]), "already_exists": True}

    # Construct possible URLs
    urls = construct_gtrnadb_url(organism)

    if manual_only:
        print("\nManual download URLs to try:")
        for url, filename in urls:
            print(f"  {url}")
        print(f"\nSave as: fastas/{organism['gtrnadb_id']}-tRNAs.fa")
        return {"success": False, "manual": True}

    # Try each URL
    print("\nAttempting automated download...")
    for url, filename in urls:
        output_path = fastas_dir / f"{organism['gtrnadb_id']}-tRNAs.fa"

        if download_file(url, output_path):
            print(f"\n✅ Successfully downloaded: {output_path}")
            return {"success": True, "file": str(output_path), "url": url}

        time.sleep(0.5)  # Be nice to the server

    # All URLs failed
    print(f"\n❌ Automated download failed for {organism_name}")
    print("\nPlease manually download from GtRNAdb:")
    print("  1. Visit: http://gtrnadb.ucsc.edu/")
    print(f"  2. Navigate to: {organism['kingdom']} → {organism_name}")
    print("  3. Download mature tRNA sequences (FASTA)")
    print(f"  4. Save as: fastas/{organism['gtrnadb_id']}-tRNAs.fa")

    return {"success": False, "error": "All download URLs failed"}


def main():
    parser = argparse.ArgumentParser(
        description="Download tRNA FASTA files from GtRNAdb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--organisms",
        help="Comma-separated list of organism IDs to download (default: all pending)",
        default=None,
    )
    parser.add_argument(
        "--manual",
        action="store_true",
        help="Only print URLs for manual download, do not attempt automated download",
    )

    args = parser.parse_args()

    # Load configuration
    config = load_config()

    # Get organisms to download
    all_organisms = config["organisms"]

    if args.organisms:
        organism_ids = args.organisms.split(",")
        organisms = [org for org in all_organisms if org["organism_id"] in organism_ids]
    else:
        # Download for all pending organisms
        organisms = [
            org
            for org in all_organisms
            if org.get("status") == "pending" and org.get("priority") == 1
        ]

    if not organisms:
        print("No organisms to download!")
        sys.exit(0)

    print(f"\n{'='*80}")
    print("GtRNAdb FASTA Download Tool")
    print(f"{'='*80}")
    print(f"Mode: {'MANUAL URLS ONLY' if args.manual else 'AUTOMATED DOWNLOAD'}")
    print(f"Organisms: {len(organisms)}")
    for org in organisms:
        print(f"  - {org['name']} ({org['organism_id']})")
    print(f"{'='*80}\n")

    # Download each organism
    results = []
    for organism in organisms:
        result = download_organism_fasta(organism, manual_only=args.manual)
        result["organism_id"] = organism["organism_id"]
        result["organism_name"] = organism["name"]
        results.append(result)

    # Summary
    print(f"\n\n{'='*80}")
    print("DOWNLOAD SUMMARY")
    print(f"{'='*80}")

    successful = [r for r in results if r["success"]]
    failed = [r for r in results if not r["success"] and not r.get("manual")]
    # Note: 'manual' results are those where --manual flag was used

    print(f"\nTotal organisms: {len(results)}")
    print(f"Already exist: {len([r for r in successful if r.get('already_exists')])}")
    print(f"Downloaded: {len([r for r in successful if not r.get('already_exists')])}")
    print(f"Failed (manual download required): {len(failed)}")

    if failed:
        print("\n⚠️  Manual download required for:")
        for r in failed:
            print(f"   - {r['organism_name']}")
        print("\nVisit http://gtrnadb.ucsc.edu/ to download these organisms manually.")

    print(f"\n{'='*80}\n")

    sys.exit(0 if len(failed) == 0 else 1)


if __name__ == "__main__":
    main()
