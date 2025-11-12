#!/usr/bin/env python3
"""
Parser for Modomics tRNA FASTA files.

This module provides functions to:
- Parse Modomics FASTA files (modified and unmodified versions)
- Extract metadata from headers
- Detect modification positions by comparing sequences
- Generate structured output for downstream analysis
"""

import re
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
import logging

from .modification_codes import ModificationCodec

logger = logging.getLogger(__name__)


@dataclass
class ModomicsTRNA:
    """Data structure for a Modomics tRNA entry."""
    modomics_id: int
    name: str
    so_term: str
    trna_type: str
    subtype: str
    anticodon: str
    cellular_localization: str
    species: str
    modified_sequence: str
    unmodified_sequence: Optional[str] = None
    modifications: Optional[List[dict]] = None

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return asdict(self)


class ModomicsParser:
    """Parser for Modomics FASTA files."""

    def __init__(self, codec: Optional[ModificationCodec] = None):
        """
        Initialize the parser.

        Args:
            codec: ModificationCodec instance for decoding modifications
        """
        self.codec = codec
        self.trnas: Dict[int, ModomicsTRNA] = {}

    def parse_header(self, header: str) -> dict:
        """
        Parse a Modomics FASTA header.

        Expected format:
        >id:1|Name:tdbR00000010|SOterm:SO:0000254|Type:tRNA|Subtype:Ala|Feature:VGC|Cellular_Localization:|Species:Escherichia coli

        Args:
            header: FASTA header line (with or without leading >)

        Returns:
            Dictionary with parsed metadata
        """
        header = header.lstrip('>')
        fields = {}

        # Split by | and parse key:value pairs
        for part in header.split('|'):
            if ':' in part:
                key, value = part.split(':', 1)
                fields[key.strip()] = value.strip()

        return {
            'modomics_id': int(fields.get('id', -1)),
            'name': fields.get('Name', ''),
            'so_term': fields.get('SOterm', ''),
            'trna_type': fields.get('Type', ''),
            'subtype': fields.get('Subtype', ''),
            'anticodon': fields.get('Feature', ''),
            'cellular_localization': fields.get('Cellular_Localization', ''),
            'species': fields.get('Species', ''),
        }

    def parse_fasta(self, fasta_path: str, sequence_type: str = 'modified') -> Dict[int, ModomicsTRNA]:
        """
        Parse a Modomics FASTA file.

        Args:
            fasta_path: Path to FASTA file
            sequence_type: Either 'modified' or 'unmodified'

        Returns:
            Dictionary mapping modomics_id to ModomicsTRNA objects
        """
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        trnas = {}
        current_header = None
        current_sequence = []

        with open(fasta_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Save previous entry
                    if current_header is not None:
                        metadata = self.parse_header(current_header)
                        seq = ''.join(current_sequence)

                        if sequence_type == 'modified':
                            trna = ModomicsTRNA(
                                modified_sequence=seq,
                                **metadata
                            )
                        else:
                            # For unmodified, we'll merge with existing entries later
                            trna = ModomicsTRNA(
                                modified_sequence='',
                                unmodified_sequence=seq,
                                **metadata
                            )

                        trnas[metadata['modomics_id']] = trna

                    # Start new entry
                    current_header = line
                    current_sequence = []
                else:
                    current_sequence.append(line)

            # Save last entry
            if current_header is not None:
                metadata = self.parse_header(current_header)
                seq = ''.join(current_sequence)

                if sequence_type == 'modified':
                    trna = ModomicsTRNA(
                        modified_sequence=seq,
                        **metadata
                    )
                else:
                    trna = ModomicsTRNA(
                        modified_sequence='',
                        unmodified_sequence=seq,
                        **metadata
                    )

                trnas[metadata['modomics_id']] = trna

        logger.info(f"Parsed {len(trnas)} tRNA sequences from {fasta_path}")
        return trnas

    def merge_sequences(self,
                       modified_trnas: Dict[int, ModomicsTRNA],
                       unmodified_trnas: Dict[int, ModomicsTRNA]) -> Dict[int, ModomicsTRNA]:
        """
        Merge modified and unmodified sequence data.

        Args:
            modified_trnas: Dictionary of tRNAs with modified sequences
            unmodified_trnas: Dictionary of tRNAs with unmodified sequences

        Returns:
            Dictionary of tRNAs with both sequences populated
        """
        merged = {}

        for mod_id, mod_trna in modified_trnas.items():
            if mod_id in unmodified_trnas:
                unmod_trna = unmodified_trnas[mod_id]
                mod_trna.unmodified_sequence = unmod_trna.unmodified_sequence
            else:
                logger.warning(f"No unmodified sequence found for Modomics ID {mod_id}")

            merged[mod_id] = mod_trna

        return merged

    def detect_modifications(self, trna: ModomicsTRNA) -> List[dict]:
        """
        Detect modification positions by comparing modified and unmodified sequences.

        Args:
            trna: ModomicsTRNA object with both sequences populated

        Returns:
            List of modification dictionaries with position and code information
        """
        if not trna.unmodified_sequence:
            logger.warning(f"Cannot detect modifications for {trna.name}: missing unmodified sequence")
            return []

        modifications = []
        mod_seq = trna.modified_sequence
        unmod_seq = trna.unmodified_sequence

        if len(mod_seq) != len(unmod_seq):
            logger.warning(
                f"Sequence length mismatch for {trna.name}: "
                f"modified={len(mod_seq)}, unmodified={len(unmod_seq)}"
            )
            # Still proceed with comparison up to shorter length
            min_len = min(len(mod_seq), len(unmod_seq))
        else:
            min_len = len(mod_seq)

        for pos in range(min_len):
            mod_char = mod_seq[pos]
            unmod_char = unmod_seq[pos]

            if mod_char != unmod_char:
                mod_info = {
                    'position': pos + 1,  # 1-based position
                    'modified_char': mod_char,
                    'unmodified_char': unmod_char,
                }

                # Add decoded information if codec available
                if self.codec:
                    decoded = self.codec.decode(mod_char)
                    if decoded:
                        mod_info['modification_name'] = decoded['name']
                        mod_info['short_name'] = decoded['short_name']
                        mod_info['reference_base'] = decoded['reference_base']
                        mod_info['modomics_db_id'] = decoded['modomics_db_id']
                    else:
                        mod_info['modification_name'] = 'Unknown'
                        logger.warning(f"Unknown modification code '{mod_char}' at position {pos+1} in {trna.name}")

                modifications.append(mod_info)

        return modifications

    def process_all_modifications(self, trnas: Dict[int, ModomicsTRNA]):
        """
        Detect modifications for all tRNAs in place.

        Args:
            trnas: Dictionary of ModomicsTRNA objects (modified in place)
        """
        for trna in trnas.values():
            trna.modifications = self.detect_modifications(trna)

        total_mods = sum(len(t.modifications or []) for t in trnas.values())
        logger.info(f"Detected {total_mods} modifications across {len(trnas)} tRNAs")

    def export_to_json(self, trnas: Dict[int, ModomicsTRNA], output_path: str):
        """
        Export parsed tRNA data to JSON.

        Args:
            trnas: Dictionary of ModomicsTRNA objects
            output_path: Path for output JSON file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to serializable format
        data = {
            str(trna_id): trna.to_dict()
            for trna_id, trna in trnas.items()
        }

        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        logger.info(f"Exported {len(trnas)} tRNAs to {output_path}")

    def get_species_list(self, trnas: Dict[int, ModomicsTRNA]) -> List[str]:
        """Get sorted list of unique species in the dataset."""
        species = set(trna.species for trna in trnas.values())
        return sorted(species)

    def get_statistics(self, trnas: Dict[int, ModomicsTRNA]) -> dict:
        """
        Generate statistics about the parsed tRNA dataset.

        Args:
            trnas: Dictionary of ModomicsTRNA objects

        Returns:
            Dictionary with statistics
        """
        species_counts = {}
        subtype_counts = {}
        total_modifications = 0

        for trna in trnas.values():
            # Count by species
            species_counts[trna.species] = species_counts.get(trna.species, 0) + 1

            # Count by subtype
            subtype_counts[trna.subtype] = subtype_counts.get(trna.subtype, 0) + 1

            # Count modifications
            if trna.modifications:
                total_modifications += len(trna.modifications)

        return {
            'total_trnas': len(trnas),
            'total_species': len(species_counts),
            'total_modifications': total_modifications,
            'avg_modifications_per_trna': total_modifications / len(trnas) if trnas else 0,
            'species_counts': species_counts,
            'subtype_counts': subtype_counts,
        }


def parse_modomics_files(modified_fasta: str,
                         unmodified_fasta: str,
                         codes_csv: str,
                         output_json: Optional[str] = None) -> Dict[int, ModomicsTRNA]:
    """
    Convenience function to parse Modomics files and detect modifications.

    Args:
        modified_fasta: Path to modified sequences FASTA
        unmodified_fasta: Path to unmodified sequences FASTA
        codes_csv: Path to modification codes CSV
        output_json: Optional path to save output JSON

    Returns:
        Dictionary of parsed ModomicsTRNA objects with modifications detected
    """
    # Load codec
    codec = ModificationCodec(codes_csv)

    # Initialize parser
    parser = ModomicsParser(codec=codec)

    # Parse both FASTA files
    logger.info("Parsing modified sequences...")
    modified_trnas = parser.parse_fasta(modified_fasta, sequence_type='modified')

    logger.info("Parsing unmodified sequences...")
    unmodified_trnas = parser.parse_fasta(unmodified_fasta, sequence_type='unmodified')

    # Merge sequences
    logger.info("Merging sequences...")
    trnas = parser.merge_sequences(modified_trnas, unmodified_trnas)

    # Detect modifications
    logger.info("Detecting modifications...")
    parser.process_all_modifications(trnas)

    # Export if requested
    if output_json:
        parser.export_to_json(trnas, output_json)

    # Print statistics
    stats = parser.get_statistics(trnas)
    logger.info(f"\n=== Statistics ===")
    logger.info(f"Total tRNAs: {stats['total_trnas']}")
    logger.info(f"Total species: {stats['total_species']}")
    logger.info(f"Total modifications: {stats['total_modifications']}")
    logger.info(f"Avg modifications per tRNA: {stats['avg_modifications_per_trna']:.2f}")

    return trnas


if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Parse Modomics tRNA FASTA files')
    parser.add_argument('--modified', required=True, help='Path to modified sequences FASTA')
    parser.add_argument('--unmodified', required=True, help='Path to unmodified sequences FASTA')
    parser.add_argument('--codes', required=True, help='Path to modification codes CSV')
    parser.add_argument('--output', help='Output JSON path (default: outputs/modomics/modomics_modifications.json)')

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    output_path = args.output or 'outputs/modomics/modomics_modifications.json'

    trnas = parse_modomics_files(
        modified_fasta=args.modified,
        unmodified_fasta=args.unmodified,
        codes_csv=args.codes,
        output_json=output_path
    )

    print(f"\n✓ Parsed {len(trnas)} tRNAs with modifications")
    print(f"✓ Output saved to: {output_path}")
