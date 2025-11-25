#!/usr/bin/env python3
"""
Modification code decoder for Modomics single-character codes.

This module provides functions to:
- Parse the modomicscodes.csv file
- Map single-character codes to modification names and metadata
- Handle bidirectional lookups (code->name and name->code)
"""

import csv
from pathlib import Path
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)


class ModificationCodec:
    """Decoder for Modomics modification codes."""

    def __init__(self, codes_csv_path: str):
        """
        Initialize the codec from a modomicscodes.csv file.

        Args:
            codes_csv_path: Path to modomicscodes.csv
        """
        self.codes_csv_path = Path(codes_csv_path)
        self.code_to_info: Dict[str, dict] = {}
        self.name_to_code: Dict[str, str] = {}
        self.short_name_to_code: Dict[str, str] = {}

        self._load_codes()

    def _load_codes(self):
        """Load and parse the modification codes CSV file."""
        if not self.codes_csv_path.exists():
            raise FileNotFoundError(f"Modification codes file not found: {self.codes_csv_path}")

        with open(self.codes_csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                code = row.get("RNAMods code (2023)", "").strip()
                if not code:
                    logger.warning(f"Skipping row with empty code: {row}")
                    continue

                info = {
                    "name": row.get("Name", "").strip(),
                    "short_name": row.get("Short Name", "").strip(),
                    "modomics_code_new": row.get("MODOMICS code new (2023)", "").strip(),
                    "moiety_type": row.get("Moiety type", "").strip(),
                    "reference_base": row.get("Reference NucleoBase", "").strip(),
                    "modomics_db_id": row.get("MODOMICS Database ID", "").strip(),
                }

                self.code_to_info[code] = info

                # Build reverse lookups
                if info["name"]:
                    self.name_to_code[info["name"]] = code
                if info["short_name"]:
                    self.short_name_to_code[info["short_name"]] = code

        logger.info(
            f"Loaded {len(self.code_to_info)} modification codes from {self.codes_csv_path}"
        )

    def decode(self, code: str) -> Optional[dict]:
        """
        Get modification information from a single-character code.

        Args:
            code: Single-character modification code

        Returns:
            Dictionary with modification information, or None if not found
        """
        return self.code_to_info.get(code)

    def get_modification_name(self, code: str) -> Optional[str]:
        """
        Get the full modification name from a code.

        Args:
            code: Single-character modification code

        Returns:
            Full modification name, or None if not found
        """
        info = self.decode(code)
        return info["name"] if info else None

    def get_short_name(self, code: str) -> Optional[str]:
        """
        Get the short modification name from a code.

        Args:
            code: Single-character modification code

        Returns:
            Short modification name (e.g., "m1A", "Y"), or None if not found
        """
        info = self.decode(code)
        return info["short_name"] if info else None

    def get_reference_base(self, code: str) -> Optional[str]:
        """
        Get the reference nucleobase for a modification.

        Args:
            code: Single-character modification code

        Returns:
            Reference base (A, C, G, or U), or None if not found
        """
        info = self.decode(code)
        return info["reference_base"] if info else None

    def encode(self, name: str) -> Optional[str]:
        """
        Get the single-character code from a modification name.

        Args:
            name: Full or short modification name

        Returns:
            Single-character code, or None if not found
        """
        # Try full name first
        code = self.name_to_code.get(name)
        if code:
            return code

        # Try short name
        return self.short_name_to_code.get(name)

    def is_modified_base(self, char: str) -> bool:
        """
        Check if a character represents a modified base.

        Args:
            char: Single character to check

        Returns:
            True if the character is a known modification code
        """
        return char in self.code_to_info

    def is_standard_base(self, char: str) -> bool:
        """
        Check if a character is a standard RNA base.

        Args:
            char: Single character to check

        Returns:
            True if the character is A, C, G, or U
        """
        return char.upper() in ["A", "C", "G", "U"]

    def get_all_codes(self) -> list:
        """Get list of all known modification codes."""
        return list(self.code_to_info.keys())

    def get_statistics(self) -> dict:
        """Get statistics about loaded modifications."""
        base_counts = {"A": 0, "C": 0, "G": 0, "U": 0, "Unknown": 0}

        for info in self.code_to_info.values():
            ref_base = info.get("reference_base", "Unknown")
            if ref_base in base_counts:
                base_counts[ref_base] += 1
            else:
                base_counts["Unknown"] += 1

        return {
            "total_modifications": len(self.code_to_info),
            "by_reference_base": base_counts,
        }


def load_modification_codec(codes_csv_path: str) -> ModificationCodec:
    """
    Convenience function to load modification codec.

    Args:
        codes_csv_path: Path to modomicscodes.csv

    Returns:
        Initialized ModificationCodec instance
    """
    return ModificationCodec(codes_csv_path)


if __name__ == "__main__":
    # Test the codec
    import sys

    if len(sys.argv) > 1:
        codec_path = sys.argv[1]
    else:
        codec_path = "docs/development/modomicscodes.csv"

    logging.basicConfig(level=logging.INFO)
    codec = load_modification_codec(codec_path)

    print("\n=== Modification Codec Statistics ===")
    stats = codec.get_statistics()
    print(f"Total modifications: {stats['total_modifications']}")
    print("\nBy reference base:")
    for base, count in stats["by_reference_base"].items():
        print(f"  {base}: {count}")

    print("\n=== Example Lookups ===")
    # Test with some example codes that should be in the file
    test_codes = ["D", "T", "Y", "K", "Ѣ"]
    for code in test_codes:
        info = codec.decode(code)
        if info:
            print(f"{code}: {info['short_name']} - {info['name'][:50]}...")
        else:
            print(f"{code}: Not found")
