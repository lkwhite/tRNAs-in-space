#!/usr/bin/env python3
"""
Phase 2: Sequence alignment and Sprinzl position mapping.

This module aligns Modomics tRNA sequences to gtRNAdb sequences
and maps modification positions to Sprinzl coordinates.
"""

import csv
import json
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from Bio import pairwise2
except ImportError:
    print("ERROR: BioPython not installed. Run: pip install biopython")
    import sys

    sys.exit(1)

logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Result of aligning a Modomics tRNA to a gtRNAdb tRNA."""

    modomics_id: int
    modomics_name: str
    gtRNAdb_trna_id: str
    alignment_score: float
    alignment_identity: float
    aligned_modomics: str
    aligned_gtRNAdb: str
    position_mapping: Dict[int, int]  # modomics_pos -> gtRNAdb_pos


@dataclass
class SprinzlMapping:
    """Complete mapping from Modomics modification to Sprinzl position."""

    modomics_id: int
    modomics_name: str
    species: str
    trna_type: str
    anticodon: str
    gtRNAdb_trna_id: str
    alignment_score: float
    alignment_identity: float
    modification_position_modomics: int
    modification_char: str
    modification_name: str
    modification_short_name: str
    unmodified_char: str
    position_gtRNAdb: int
    sprinzl_index: int
    sprinzl_label: str
    region: str


class SpeciesNameMapper:
    """Maps species names between Modomics and gtRNAdb formats."""

    # Manual mapping for key species
    SPECIES_MAP = {
        "Escherichia coli": "ecoliK12",
        "Escherichia coli K-12": "ecoliK12",
        "Escherichia coli str. K-12": "ecoliK12",
        "Saccharomyces cerevisiae": "sacCer",
        "Homo sapiens": "hg38",
    }

    @classmethod
    def modomics_to_gtRNAdb(cls, modomics_species: str) -> Optional[str]:
        """Convert Modomics species name to gtRNAdb organism ID."""
        return cls.SPECIES_MAP.get(modomics_species)

    @classmethod
    def get_supported_species(cls) -> List[str]:
        """Get list of supported species names (Modomics format)."""
        return list(cls.SPECIES_MAP.keys())


class AnticodonNormalizer:
    """Normalizes anticodons to handle modified bases."""

    # Map modified anticodon codes to canonical bases
    MODIFIED_BASES = {
        "I": "G",  # Inosine
        "V": "U",  # Various uridine modifications
        "Y": "U",  # Pseudouridine
        "D": "U",  # Dihydrouridine
        "T": "U",  # Ribothymidine
        "m1G": "G",
        "m2G": "G",
        "m7G": "G",
    }

    @classmethod
    def normalize(cls, anticodon: str) -> str:
        """
        Normalize anticodon by replacing modified bases with canonical bases.

        Args:
            anticodon: Original anticodon (may contain modification codes)

        Returns:
            Normalized anticodon with only A, C, G, U
        """
        normalized = anticodon
        for mod, canonical in cls.MODIFIED_BASES.items():
            normalized = normalized.replace(mod, canonical)

        # Remove any remaining non-standard characters
        valid = set("ACGU")
        normalized = "".join(c if c in valid else "N" for c in normalized)

        return normalized


class GTRNAdbLoader:
    """Loads gtRNAdb FASTA and global coordinate files."""

    def __init__(self, fasta_dir: str, output_dir: str):
        """
        Initialize loader.

        Args:
            fasta_dir: Directory containing gtRNAdb FASTA files
            output_dir: Directory containing global_coords.tsv files
        """
        self.fasta_dir = Path(fasta_dir)
        self.output_dir = Path(output_dir)

    def load_fasta(self, organism_id: str) -> Dict[str, str]:
        """
        Load gtRNAdb FASTA file for an organism.

        Args:
            organism_id: Organism ID (e.g., 'ecoliK12')

        Returns:
            Dictionary mapping tRNA ID to sequence (RNA format with U)
        """
        # Try different possible FASTA filenames
        possible_names = [
            f"{organism_id}-tRNAs.fa",
            f"{organism_id}-mature-tRNAs.fa",
            f"{organism_id}-mito-and-nuclear-tRNAs.fa",  # Yeast/human format
            f"{organism_id}MG1655-tRNAs.fa",  # E. coli special case
        ]

        fasta_file = None
        for name in possible_names:
            path = self.fasta_dir / name
            if path.exists():
                fasta_file = path
                break

        if not fasta_file:
            raise FileNotFoundError(
                f"Could not find FASTA file for {organism_id} in {self.fasta_dir}"
            )

        logger.info(f"Loading FASTA: {fasta_file}")

        sequences = {}
        current_id = None
        current_seq: list[str] = []

        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # Save previous sequence
                    if current_id:
                        # Convert DNA (T) to RNA (U)
                        seq = "".join(current_seq).replace("T", "U")
                        sequences[current_id] = seq

                    # Start new sequence
                    current_id = line[1:].split()[0]  # Take first word after >
                    current_seq = []
                else:
                    current_seq.append(line)

            # Save last sequence
            if current_id:
                seq = "".join(current_seq).replace("T", "U")
                sequences[current_id] = seq

        logger.info(f"Loaded {len(sequences)} tRNA sequences")
        return sequences

    def load_global_coords(self, organism_id: str) -> Dict[str, List[dict]]:
        """
        Load global coordinates TSV file for an organism.

        Args:
            organism_id: Organism ID (e.g., 'ecoliK12')

        Returns:
            Dictionary mapping tRNA ID to list of position records
        """
        coords_file = self.output_dir / f"{organism_id}_global_coords.tsv"

        if not coords_file.exists():
            raise FileNotFoundError(f"Global coords file not found: {coords_file}")

        logger.info(f"Loading global coords: {coords_file}")

        coords = defaultdict(list)

        with open(coords_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                trna_id = row["trna_id"]
                coords[trna_id].append(
                    {
                        "seq_index": int(row["seq_index"]),
                        "sprinzl_index": int(row["sprinzl_index"]),
                        "sprinzl_label": row["sprinzl_label"],
                        "residue": row["residue"],
                        "region": row["region"],
                    }
                )

        logger.info(f"Loaded coordinates for {len(coords)} tRNAs")
        return dict(coords)


class ModomicsAligner:
    """Aligns Modomics tRNA sequences to gtRNAdb sequences."""

    def __init__(
        self, match: int = 2, mismatch: int = -1, gap_open: float = -2, gap_extend: float = -0.5
    ):
        """
        Initialize aligner with scoring parameters.

        Args:
            match: Score for matching bases
            mismatch: Penalty for mismatches
            gap_open: Penalty for opening a gap
            gap_extend: Penalty for extending a gap
        """
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend

    def align_sequences(self, seq1: str, seq2: str) -> Tuple[str, str, float]:
        """
        Perform global alignment of two sequences.

        Args:
            seq1: First sequence
            seq2: Second sequence

        Returns:
            Tuple of (aligned_seq1, aligned_seq2, score)
        """
        alignments = pairwise2.align.globalms(
            seq1,
            seq2,
            self.match,
            self.mismatch,
            self.gap_open,
            self.gap_extend,
            one_alignment_only=True,
        )

        if not alignments:
            return seq1, seq2, 0.0

        alignment = alignments[0]
        return alignment.seqA, alignment.seqB, alignment.score

    def calculate_identity(self, aligned1: str, aligned2: str) -> float:
        """
        Calculate sequence identity percentage.

        Args:
            aligned1: First aligned sequence (may contain gaps)
            aligned2: Second aligned sequence (may contain gaps)

        Returns:
            Identity as percentage (0-100)
        """
        matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != "-")
        total = len([c for c in aligned1 if c != "-"])

        if total == 0:
            return 0.0

        return (matches / total) * 100.0

    def create_position_mapping(
        self, aligned_modomics: str, aligned_gtRNAdb: str
    ) -> Dict[int, int]:
        """
        Create mapping from Modomics positions to gtRNAdb positions.

        Args:
            aligned_modomics: Aligned Modomics sequence (with gaps)
            aligned_gtRNAdb: Aligned gtRNAdb sequence (with gaps)

        Returns:
            Dictionary mapping 1-based Modomics position to 1-based gtRNAdb position
        """
        mapping = {}
        modomics_pos = 0
        gtRNAdb_pos = 0

        for mod_char, gtR_char in zip(aligned_modomics, aligned_gtRNAdb):
            if mod_char != "-":
                modomics_pos += 1
            if gtR_char != "-":
                gtRNAdb_pos += 1

            # Only map if both positions are not gaps
            if mod_char != "-" and gtR_char != "-":
                mapping[modomics_pos] = gtRNAdb_pos

        return mapping

    def find_best_match(
        self,
        modomics_seq: str,
        modomics_anticodon: str,
        modomics_subtype: str,
        gtRNAdb_sequences: Dict[str, str],
        min_identity: float = 80.0,
    ) -> Optional[AlignmentResult]:
        """
        Find best matching gtRNAdb sequence for a Modomics tRNA.

        Args:
            modomics_seq: Modomics unmodified sequence
            modomics_anticodon: Modomics anticodon
            modomics_subtype: Modomics tRNA subtype (amino acid)
            gtRNAdb_sequences: Dictionary of gtRNAdb sequences
            min_identity: Minimum sequence identity threshold (%)

        Returns:
            AlignmentResult or None if no good match found
        """
        # Normalize anticodon
        norm_anticodon = AnticodonNormalizer.normalize(modomics_anticodon)

        best_result = None
        best_score = -float("inf")

        for trna_id, gtRNAdb_seq in gtRNAdb_sequences.items():
            # Quick filter: check if anticodon and subtype match
            # tRNA ID formats:
            #   tRNA-{Subtype}-{Anticodon}-{Copy}-{Variant}  (E. coli)
            #   nuc-tRNA-{Subtype}-{Anticodon}-{Copy}-{Variant}  (yeast/human nuclear)
            #   mito-tRNA-{Subtype}-{Anticodon}-{Copy}-{Variant}  (yeast/human mito)
            parts = trna_id.split("-")

            # Handle prefix (nuc-, mito-, or no prefix)
            if parts[0] in ["nuc", "mito"]:
                # Format: nuc-tRNA-Ala-AGC-1-1
                if len(parts) >= 5:
                    gtR_subtype = parts[2]
                    gtR_anticodon = parts[3]
                else:
                    continue
            else:
                # Format: tRNA-Ala-AGC-1-1
                if len(parts) >= 3:
                    gtR_subtype = parts[1]
                    gtR_anticodon = parts[2]
                else:
                    continue

            # Skip if amino acid doesn't match
            if gtR_subtype.lower() != modomics_subtype.lower():
                continue

            # Skip if anticodon doesn't match (allowing for normalized anticodon)
            if gtR_anticodon.upper() != norm_anticodon.upper():
                continue

            # Perform alignment
            aligned_mod, aligned_gtR, score = self.align_sequences(modomics_seq, gtRNAdb_seq)
            identity = self.calculate_identity(aligned_mod, aligned_gtR)

            # Keep best match
            if score > best_score and identity >= min_identity:
                best_score = score
                position_mapping = self.create_position_mapping(aligned_mod, aligned_gtR)

                best_result = AlignmentResult(
                    modomics_id=0,  # Will be filled by caller
                    modomics_name="",  # Will be filled by caller
                    gtRNAdb_trna_id=trna_id,
                    alignment_score=score,
                    alignment_identity=identity,
                    aligned_modomics=aligned_mod,
                    aligned_gtRNAdb=aligned_gtR,
                    position_mapping=position_mapping,
                )

        return best_result


class SprinzlMapper:
    """Maps Modomics modifications to Sprinzl positions."""

    def __init__(self, modomics_json: str, fasta_dir: str, output_dir: str):
        """
        Initialize mapper.

        Args:
            modomics_json: Path to modomics_modifications.json
            fasta_dir: Directory containing gtRNAdb FASTA files
            output_dir: Directory containing global_coords.tsv files
        """
        self.modomics_json = Path(modomics_json)
        self.loader = GTRNAdbLoader(fasta_dir, output_dir)
        self.aligner = ModomicsAligner()

        # Load Modomics data
        with open(self.modomics_json, "r") as f:
            self.modomics_data = json.load(f)

        logger.info(f"Loaded {len(self.modomics_data)} Modomics tRNAs")

    def process_species(
        self, modomics_species: str, min_identity: float = 80.0
    ) -> List[SprinzlMapping]:
        """
        Process all tRNAs for a species and create Sprinzl mappings.

        Args:
            modomics_species: Species name in Modomics format
            min_identity: Minimum alignment identity threshold (%)

        Returns:
            List of SprinzlMapping objects
        """
        # Map species name
        organism_id = SpeciesNameMapper.modomics_to_gtRNAdb(modomics_species)
        if not organism_id:
            logger.warning(f"No gtRNAdb organism ID for species: {modomics_species}")
            return []

        logger.info(f"\n=== Processing {modomics_species} → {organism_id} ===")

        # Load gtRNAdb data
        try:
            gtRNAdb_sequences = self.loader.load_fasta(organism_id)
            global_coords = self.loader.load_global_coords(organism_id)
        except FileNotFoundError as e:
            logger.error(f"Could not load gtRNAdb data: {e}")
            return []

        # Filter Modomics tRNAs for this species
        species_trnas = {
            trna_id: trna_data
            for trna_id, trna_data in self.modomics_data.items()
            if trna_data["species"] == modomics_species
        }

        logger.info(f"Found {len(species_trnas)} Modomics tRNAs for {modomics_species}")

        # Process each tRNA
        mappings = []
        aligned_count = 0

        for trna_id, trna_data in species_trnas.items():
            # Skip if no unmodified sequence
            if not trna_data.get("unmodified_sequence"):
                logger.warning(f"Skipping {trna_data['name']}: no unmodified sequence")
                continue

            # Skip if no modifications
            if not trna_data.get("modifications"):
                continue

            # Find best alignment
            alignment = self.aligner.find_best_match(
                modomics_seq=trna_data["unmodified_sequence"],
                modomics_anticodon=trna_data["anticodon"],
                modomics_subtype=trna_data["subtype"],
                gtRNAdb_sequences=gtRNAdb_sequences,
                min_identity=min_identity,
            )

            if not alignment:
                logger.debug(
                    f"No alignment for {trna_data['name']} ({trna_data['subtype']}-{trna_data['anticodon']})"
                )
                continue

            aligned_count += 1

            # Get Sprinzl coordinates for the matched gtRNAdb tRNA
            if alignment.gtRNAdb_trna_id not in global_coords:
                logger.warning(f"No global coords for {alignment.gtRNAdb_trna_id}")
                continue

            coords = global_coords[alignment.gtRNAdb_trna_id]

            # Map each modification
            for mod in trna_data["modifications"]:
                modomics_pos = mod["position"]

                # Get gtRNAdb position from alignment
                if modomics_pos not in alignment.position_mapping:
                    logger.debug(f"Position {modomics_pos} not in alignment mapping")
                    continue

                gtRNAdb_pos = alignment.position_mapping[modomics_pos]

                # Find Sprinzl position
                sprinzl_info = None
                for coord in coords:
                    if coord["seq_index"] == gtRNAdb_pos:
                        sprinzl_info = coord
                        break

                if not sprinzl_info:
                    logger.debug(f"No Sprinzl position for gtRNAdb pos {gtRNAdb_pos}")
                    continue

                # Create mapping
                mapping = SprinzlMapping(
                    modomics_id=int(trna_id),
                    modomics_name=trna_data["name"],
                    species=modomics_species,
                    trna_type=trna_data["subtype"],
                    anticodon=trna_data["anticodon"],
                    gtRNAdb_trna_id=alignment.gtRNAdb_trna_id,
                    alignment_score=alignment.alignment_score,
                    alignment_identity=alignment.alignment_identity,
                    modification_position_modomics=modomics_pos,
                    modification_char=mod["modified_char"],
                    modification_name=mod.get("modification_name", "Unknown"),
                    modification_short_name=mod.get("short_name", ""),
                    unmodified_char=mod["unmodified_char"],
                    position_gtRNAdb=gtRNAdb_pos,
                    sprinzl_index=sprinzl_info["sprinzl_index"],
                    sprinzl_label=sprinzl_info["sprinzl_label"],
                    region=sprinzl_info["region"],
                )

                mappings.append(mapping)

        logger.info(f"Successfully aligned {aligned_count}/{len(species_trnas)} tRNAs")
        logger.info(f"Created {len(mappings)} modification→Sprinzl mappings")

        return mappings

    def export_to_tsv(self, mappings: List[SprinzlMapping], output_path: str):
        """
        Export mappings to TSV file.

        Args:
            mappings: List of SprinzlMapping objects
            output_path: Output TSV file path
        """
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        fieldnames = [
            "modomics_id",
            "modomics_name",
            "species",
            "trna_type",
            "anticodon",
            "gtRNAdb_trna_id",
            "alignment_score",
            "alignment_identity",
            "modification_position_modomics",
            "modification_char",
            "modification_name",
            "modification_short_name",
            "unmodified_char",
            "position_gtRNAdb",
            "sprinzl_index",
            "sprinzl_label",
            "region",
        ]

        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()

            for mapping in mappings:
                writer.writerow(
                    {
                        "modomics_id": mapping.modomics_id,
                        "modomics_name": mapping.modomics_name,
                        "species": mapping.species,
                        "trna_type": mapping.trna_type,
                        "anticodon": mapping.anticodon,
                        "gtRNAdb_trna_id": mapping.gtRNAdb_trna_id,
                        "alignment_score": f"{mapping.alignment_score:.1f}",
                        "alignment_identity": f"{mapping.alignment_identity:.1f}",
                        "modification_position_modomics": mapping.modification_position_modomics,
                        "modification_char": mapping.modification_char,
                        "modification_name": mapping.modification_name,
                        "modification_short_name": mapping.modification_short_name,
                        "unmodified_char": mapping.unmodified_char,
                        "position_gtRNAdb": mapping.position_gtRNAdb,
                        "sprinzl_index": mapping.sprinzl_index,
                        "sprinzl_label": mapping.sprinzl_label,
                        "region": mapping.region,
                    }
                )

        logger.info(f"Exported {len(mappings)} mappings to {out_path}")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Align Modomics tRNAs to gtRNAdb and map to Sprinzl positions"
    )
    parser.add_argument(
        "--modomics-json", required=True, help="Path to modomics_modifications.json"
    )
    parser.add_argument(
        "--fasta-dir", default="fastas", help="Directory containing gtRNAdb FASTA files"
    )
    parser.add_argument(
        "--output-dir", default="outputs", help="Directory containing global_coords.tsv files"
    )
    parser.add_argument(
        "--species", help='Species to process (Modomics format, e.g., "Escherichia coli")'
    )
    parser.add_argument("--all-species", action="store_true", help="Process all supported species")
    parser.add_argument("--output", help="Output TSV file path")
    parser.add_argument(
        "--min-identity",
        type=float,
        default=80.0,
        help="Minimum alignment identity threshold (default: 80.0)",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Initialize mapper
    mapper = SprinzlMapper(
        modomics_json=args.modomics_json, fasta_dir=args.fasta_dir, output_dir=args.output_dir
    )

    # Determine which species to process
    if args.all_species:
        species_list = SpeciesNameMapper.get_supported_species()
    elif args.species:
        species_list = [args.species]
    else:
        print("ERROR: Must specify --species or --all-species")
        return 1

    # Process species
    all_mappings = []
    for species in species_list:
        mappings = mapper.process_species(species, min_identity=args.min_identity)
        all_mappings.extend(mappings)

    # Export results
    if all_mappings:
        output_path = args.output or "outputs/modomics/modomics_to_sprinzl_mapping.tsv"
        mapper.export_to_tsv(all_mappings, output_path)

        print(f"\n✓ Successfully created {len(all_mappings)} modification→Sprinzl mappings")
        print(f"✓ Output saved to: {output_path}")
    else:
        print("\n✗ No mappings created")
        return 1

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
