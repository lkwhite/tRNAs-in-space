#!/usr/bin/env python3
"""
Scraper for Modomics mitochondrial tRNA data.

Modomics bulk FASTA downloads exclude mitochondrial tRNAs, but individual
sequence pages contain the data. This script fetches those pages and
extracts the modification data.

Usage:
    python scripts/modomics/scrape_mito_trnas.py --output outputs/modomics/sacCer_mito_trnas.json
"""

import argparse
import json
import logging
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import requests
    from bs4 import BeautifulSoup
except ImportError:
    print("Required packages not installed. Run:")
    print("  pip install requests beautifulsoup4")
    exit(1)

from .modification_codes import ModificationCodec

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Yeast mitochondrial tRNAs from Modomics website
# Format: (modomics_id, tdbR_name, amino_acid, anticodon)
# Anticodons normalized to standard ACGU (Modomics uses !, $, N for modified bases)
YEAST_MITO_TRNAS = [
    (146, "tdbR00000164", "Ile", "GAU"),
    (147, "tdbR00000237", "Leu", "UAA"),   # !AA in Modomics notation
    (148, "tdbR00000183", "Lys", "UUU"),   # $UU in Modomics notation
    (149, "tdbR00000280", "Met", "CAU"),
    (153, "tdbR00000395", "Ser", "GCU"),
    (157, "tdbR00000551", "Tyr", "GUA"),
    (159, "tdbR00000518", "Ini", "CAU"),   # Initiator Met
    (727, "tdbR00000075", "Phe", "GAA"),
    (734, "tdbR00000123", "Gly", "UCC"),   # NCC in Modomics
    (740, "tdbR00000142", "His", "GUG"),
    (792, "tdbR00000321", "Pro", "UGG"),
    (798, "tdbR00000361", "Arg", "ACG"),
    (799, "tdbR00000362", "Arg", "UCU"),   # NCU in Modomics
    (806, "tdbR00000396", "Ser", "UGA"),
    (807, "tdbR00000397", "Ser", "UGA"),
    (818, "tdbR00000441", "Thr", "UAG"),
    (831, "tdbR00000490", "Trp", "UCA"),   # !CA in Modomics
]

BASE_URL = "https://genesilico.pl/modomics/sequences"

# Modification code to unmodified base mapping
# These are the common single-character codes used in Modomics sequences
MOD_CODE_TO_BASE = {
    "D": "U",  # dihydrouridine
    "P": "U",  # pseudouridine (pY)
    "T": "U",  # ribothymidine (m5U)
    "6": "A",  # t6A (threonylcarbamoyladenosine)
    "7": "G",  # m7G
    "K": "G",  # m1G
    "Y": "U",  # pseudouridine (alternate)
    "R": "G",  # m2,2G (dimethylguanosine)
    "L": "G",  # m2G
    "N": "G",  # ms2t6A or similar
    "Ч": "A",  # i6A (isopentenyladenosine)
    "Ѣ": "A",  # m1A
    "ʍ": "U",  # mcm5U
    "ǻ": "A",  # t6A variant
    "B": "C",  # Cm (2'-O-methylcytidine)
    "#": "G",  # Gm (2'-O-methylguanosine)
    "?": "C",  # m5C
    "1": "A",  # m1A
    "8": "A",  # ms2t6A
    "O": "A",  # m1I (1-methylinosine from A)
    "I": "A",  # inosine
    "H": "U",  # s2U or similar
    "!": "U",  # cmnm5U (5-carboxymethylaminomethyluridine)
    "$": "U",  # mcm5s2U or cmnm5s2U
    "2": "U",  # s2U (2-thiouridine)
}


@dataclass
class ModomicsTRNA:
    """Data structure matching existing modomics_modifications.json format."""
    modomics_id: int
    name: str
    so_term: str
    trna_type: str
    subtype: str
    anticodon: str
    cellular_localization: str
    species: str
    modified_sequence: str
    unmodified_sequence: str
    modifications: List[dict]

    def to_dict(self) -> dict:
        return asdict(self)


class ModomicsScraper:
    """Scraper for individual Modomics sequence pages."""

    def __init__(self, codec: Optional[ModificationCodec] = None, delay: float = 1.5):
        """
        Initialize scraper.

        Args:
            codec: ModificationCodec for looking up modification details
            delay: Seconds to wait between requests (be polite to server)
        """
        self.codec = codec
        self.delay = delay
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "tRNAs-in-space research project (academic use)"
        })

    def fetch_page(self, modomics_id: int) -> Optional[str]:
        """Fetch HTML content for a Modomics sequence page."""
        url = f"{BASE_URL}/{modomics_id}"
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            logger.error(f"Failed to fetch {url}: {e}")
            return None

    def parse_sequence_page(self, html: str, modomics_id: int,
                           amino_acid: str, anticodon: str) -> Optional[ModomicsTRNA]:
        """
        Parse a Modomics sequence page.

        Args:
            html: Raw HTML content
            modomics_id: Modomics database ID
            amino_acid: Amino acid type (Ile, Leu, etc.)
            anticodon: Anticodon sequence

        Returns:
            ModomicsTRNA object with extracted data
        """
        soup = BeautifulSoup(html, "html.parser")

        # Extract tdbR name from "Original Source" link or table
        tdbr_name = ""
        source_link = soup.find("a", href=re.compile(r"tpsic\.igcz\.poznan\.pl"))
        if source_link:
            tdbr_name = source_link.text.strip()

        # Extract SO term
        so_term = ""
        so_link = soup.find("a", href=re.compile(r"sequenceontology\.org"))
        if so_link:
            so_term = so_link.text.strip()

        # Find modified sequence from full text
        # The UNICODE encoded sequence contains single-character modification codes
        full_text = soup.get_text()

        # Look for the modified sequence - it's a string of mostly ACGU with some mod codes
        # Pattern: starts with A/G/C/U, contains mix of bases and mod codes, ends with CCA
        # Special chars: ! ($=cmnm5s2U), $ (mcm5s2U), etc.
        modified_seq = ""
        seq_match = re.search(
            r'([ACGU][ACGUDPTYKRLNБЧѢʍǻBH#?1678OI!\$2]{50,100}CCA)',
            full_text
        )
        if seq_match:
            modified_seq = seq_match.group(1)

        # Derive unmodified sequence by replacing modification codes with bases
        unmodified_seq = self._derive_unmodified_sequence(modified_seq)

        # Extract modifications from the modification position table
        modifications = []
        tables = soup.find_all("table", class_=["table", "table-bordered"])

        for table in tables:
            rows = table.find_all("tr")
            if len(rows) >= 3:
                # Check if this is a modification positions table
                header_text = rows[0].get_text()
                if "Position" in header_text:
                    # Parse the table - positions are in row 1, modifications in row 2
                    pos_cells = rows[1].find_all(["td", "th"])
                    mod_cells = rows[2].find_all(["td", "th"])

                    if len(pos_cells) > 1 and len(mod_cells) > 1:
                        # First cell is label, rest are data
                        positions = []
                        for cell in pos_cells[1:]:
                            try:
                                positions.append(int(cell.get_text().strip()))
                            except ValueError:
                                continue

                        mod_names = []
                        for cell in mod_cells[1:]:
                            mod_names.append(cell.get_text().strip())

                        # Match positions to modifications
                        for i, (pos, mod_name) in enumerate(zip(positions, mod_names)):
                            if not mod_name or mod_name == "":
                                continue

                            # Look up modification by short name first (table shows short names)
                            mod_info = None
                            mod_char = ""

                            if self.codec:
                                # First try looking up by short name (e.g., "Y", "t6A", "m5U")
                                code = self.codec.short_name_to_code.get(mod_name)
                                if code:
                                    mod_info = self.codec.decode(code)
                                    mod_char = code
                                else:
                                    # If short name lookup failed, try as single-char code
                                    if len(mod_name) == 1:
                                        mod_info = self.codec.decode(mod_name)
                                        if mod_info:
                                            mod_char = mod_name

                            # Get unmodified base from codec or fallback mapping
                            unmod_base = "N"
                            if mod_info:
                                unmod_base = mod_info.get("reference_base", "N")
                            elif mod_char in MOD_CODE_TO_BASE:
                                unmod_base = MOD_CODE_TO_BASE[mod_char]

                            modifications.append({
                                "position": pos,
                                "modified_char": mod_char,
                                "unmodified_char": unmod_base,
                                "modification_name": mod_info["name"] if mod_info else mod_name,
                                "short_name": mod_info["short_name"] if mod_info else mod_name,
                                "reference_base": unmod_base,
                                "modomics_db_id": mod_info["modomics_db_id"] if mod_info else "",
                            })

                        # Only use first matching table
                        if modifications:
                            break

        # If no table found, extract from sequence comparison
        if not modifications and modified_seq:
            modifications = self._extract_mods_from_sequences(modified_seq, unmodified_seq)

        return ModomicsTRNA(
            modomics_id=modomics_id,
            name=tdbr_name,
            so_term=so_term,
            trna_type="tRNA",
            subtype=amino_acid,
            anticodon=anticodon,
            cellular_localization="mitochondrion",
            species="Saccharomyces cerevisiae",
            modified_sequence=modified_seq,
            unmodified_sequence=unmodified_seq,
            modifications=modifications,
        )

    def _derive_unmodified_sequence(self, modified_seq: str) -> str:
        """Derive unmodified sequence by replacing modification codes with bases."""
        unmodified = ""
        for char in modified_seq:
            if char in MOD_CODE_TO_BASE:
                unmodified += MOD_CODE_TO_BASE[char]
            else:
                unmodified += char
        return unmodified

    def _extract_mods_from_sequences(self, modified: str, unmodified: str) -> List[dict]:
        """Extract modifications by comparing modified and unmodified sequences."""
        mods = []
        # Sequences should be same length after removing modification markers
        for i, (m, u) in enumerate(zip(modified, unmodified)):
            if m != u and m not in "ACGU":
                mod_info = self.codec.decode(m) if self.codec else None
                mods.append({
                    "position": i + 1,  # 1-indexed
                    "modified_char": m,
                    "unmodified_char": u,
                    "modification_name": mod_info["name"] if mod_info else f"unknown ({m})",
                    "short_name": mod_info["short_name"] if mod_info else m,
                    "reference_base": mod_info["reference_base"] if mod_info else u,
                    "modomics_db_id": mod_info["modomics_db_id"] if mod_info else "",
                })
        return mods

    def scrape_by_tdbr(self, tdbr_name: str, amino_acid: str, anticodon: str) -> Optional[ModomicsTRNA]:
        """
        Search for a tRNA by its tdbR name and scrape it.

        The tdbR name is used to find the Modomics ID via search.
        """
        # Search URL
        search_url = f"{BASE_URL}/?search={tdbr_name}"
        try:
            response = self.session.get(search_url, timeout=30)
            response.raise_for_status()

            soup = BeautifulSoup(response.text, "html.parser")

            # Find link to the specific sequence
            link = soup.find("a", href=re.compile(rf"/sequences/\d+"))
            if link:
                match = re.search(r"/sequences/(\d+)", link["href"])
                if match:
                    modomics_id = int(match.group(1))
                    time.sleep(self.delay)
                    return self.scrape_by_id(modomics_id, amino_acid, anticodon)

            logger.warning(f"Could not find Modomics ID for {tdbr_name}")
            return None

        except requests.RequestException as e:
            logger.error(f"Search failed for {tdbr_name}: {e}")
            return None

    def scrape_by_id(self, modomics_id: int, amino_acid: str, anticodon: str) -> Optional[ModomicsTRNA]:
        """Scrape a tRNA by its Modomics ID."""
        html = self.fetch_page(modomics_id)
        if not html:
            return None
        return self.parse_sequence_page(html, modomics_id, amino_acid, anticodon)

    def scrape_all_yeast_mito(self) -> List[ModomicsTRNA]:
        """Scrape all yeast mitochondrial tRNAs."""
        results = []

        for entry in YEAST_MITO_TRNAS:
            modomics_id, tdbr_name, amino_acid, anticodon = entry

            logger.info(f"Scraping {amino_acid}-{anticodon} ({tdbr_name})...")

            if modomics_id:
                trna = self.scrape_by_id(modomics_id, amino_acid, anticodon)
            else:
                trna = self.scrape_by_tdbr(tdbr_name, amino_acid, anticodon)

            if trna:
                results.append(trna)
                logger.info(f"  Found: {len(trna.modifications)} modifications")
            else:
                logger.warning(f"  Failed to scrape")

            time.sleep(self.delay)

        return results


def main():
    parser = argparse.ArgumentParser(description="Scrape Modomics mitochondrial tRNA data")
    parser.add_argument(
        "--output",
        default="outputs/modomics/sacCer_mito_trnas.json",
        help="Output JSON file"
    )
    parser.add_argument(
        "--codes-csv",
        default="docs/archive/modomics-integration/modomicscodes.csv",
        help="Path to modomicscodes.csv"
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=1.5,
        help="Delay between requests in seconds"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Don't save output, just print what would be scraped"
    )
    args = parser.parse_args()

    # Load modification codec
    codec = None
    if Path(args.codes_csv).exists():
        codec = ModificationCodec(args.codes_csv)
        logger.info(f"Loaded {len(codec.code_to_info)} modification codes")
    else:
        logger.warning(f"Modification codes file not found: {args.codes_csv}")

    # Create scraper
    scraper = ModomicsScraper(codec=codec, delay=args.delay)

    if args.dry_run:
        print("Would scrape the following yeast mitochondrial tRNAs:")
        for entry in YEAST_MITO_TRNAS:
            modomics_id, tdbr_name, amino_acid, anticodon = entry
            print(f"  - {amino_acid}-{anticodon}: {tdbr_name} (ID: {modomics_id or 'search'})")
        return

    # Scrape all
    results = scraper.scrape_all_yeast_mito()

    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    output_data = {str(t.modomics_id): t.to_dict() for t in results}

    with open(output_path, "w") as f:
        json.dump(output_data, f, indent=2)

    logger.info(f"Saved {len(results)} tRNAs to {output_path}")

    # Summary
    total_mods = sum(len(t.modifications) for t in results)
    print(f"\nSummary:")
    print(f"  tRNAs scraped: {len(results)}/{len(YEAST_MITO_TRNAS)}")
    print(f"  Total modifications: {total_mods}")


if __name__ == "__main__":
    main()
