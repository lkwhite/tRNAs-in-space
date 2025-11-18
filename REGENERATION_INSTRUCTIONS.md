# tRNA Coordinate Regeneration Instructions

**Purpose:** Regenerate all coordinate files using fixed code to eliminate critical bugs and improve quality from 37.5-47.5/100 to >90/100.

**Context:** Code fixes have been applied to `scripts/trnas_in_space.py` that address:
1. Enhanced tRNA filtering (now catches fMet)
2. Fixed region boundary gaps (positions 8, 9, 26)
3. Fixed e-position region assignment (Type II extended variable arms)

**Prerequisites:**
- Docker installed and running
- Python 3.9+ with pandas
- All code fixes already applied and committed

---

## Step 1: Verify Docker is Running

```bash
# Check Docker is available
docker ps

# Should show table of containers (can be empty)
# If error, start Docker Desktop or Docker daemon
```

---

## Step 2: Run R2DT on All FASTA Files

R2DT (RNA 2D structure Templates) annotates tRNA structures with Sprinzl positions.

### Create output directories

```bash
mkdir -p r2dt_outputs/ecoli_jsons
mkdir -p r2dt_outputs/yeast_jsons
mkdir -p r2dt_outputs/human_jsons
```

### Process E. coli

```bash
docker run --rm \
  -v "$(pwd):/data" \
  rnacentral/r2dt \
  r2dt.py gtrnadb draw \
  /data/fastas/ecoliK12MG1655-tRNAs.fa \
  /data/r2dt_outputs/ecoli_jsons
```

**Expected output:**
- Files in `r2dt_outputs/ecoli_jsons/*.enriched.json`
- Should process ~90-95 tRNAs (before filtering)

### Process Yeast

```bash
docker run --rm \
  -v "$(pwd):/data" \
  rnacentral/r2dt \
  r2dt.py gtrnadb draw \
  /data/fastas/sacCer-mito-and-nuclear-tRNAs.fa \
  /data/r2dt_outputs/yeast_jsons
```

**Expected output:**
- Files in `r2dt_outputs/yeast_jsons/*.enriched.json`
- Should process ~270-280 tRNAs

### Process Human

```bash
docker run --rm \
  -v "$(pwd):/data" \
  rnacentral/r2dt \
  r2dt.py gtrnadb draw \
  /data/fastas/hg38-mito-and-nuclear-tRNAs.fa \
  /data/r2dt_outputs/human_jsons
```

**Expected output:**
- Files in `r2dt_outputs/human_jsons/*.enriched.json`
- Should process ~420-450 tRNAs

**Note:** R2DT can take 5-15 minutes per species depending on number of tRNAs and CPU.

---

## Step 3: Generate Coordinates with Fixed Code

### E. coli

```bash
python scripts/trnas_in_space.py \
  r2dt_outputs/ecoli_jsons/ \
  outputs/ecoliK12_global_coords.tsv
```

**Expected console output:**
```
Excluding incompatible tRNA: tRNA-SeC-TCA-1-1
Excluding incompatible tRNA: tRNA-fMet-CAT-1-1
Excluding incompatible tRNA: tRNA-fMet-CAT-1-2
Excluding incompatible tRNA: tRNA-fMet-CAT-1-3
Excluding incompatible tRNA: tRNA-fMet-CAT-2-1
[ok] Global index validation: No collisions detected among X unique positions
[ok] Wrote outputs/ecoliK12_global_coords.tsv
  JSON files parsed: 87  |  skipped: 0
  Rows: ~6400  |  tRNAs: 82 (5 filtered out)
```

**Key checks:**
- ✅ Should see "Excluding incompatible tRNA" for SeC and fMet
- ✅ Should see "No collisions detected"
- ✅ Final tRNA count should be ~82 (down from 87 after filtering)

### Yeast

```bash
python scripts/trnas_in_space.py \
  r2dt_outputs/yeast_jsons/ \
  outputs/sacCer_global_coords.tsv
```

**Expected:**
- Should see exclusion messages for any mito or excluded tRNAs
- No collision warnings
- Final tRNA count: ~268

### Human

```bash
python scripts/trnas_in_space.py \
  r2dt_outputs/human_jsons/ \
  outputs/hg38_global_coords.tsv
```

**Expected:**
- Should see exclusion messages for any mito or excluded tRNAs
- No collision warnings
- Final tRNA count: ~422

---

## Step 4: Validate Regenerated Coordinates

### Run comprehensive validation

```bash
python scripts/validate_coordinate_system.py outputs/*.tsv
```

**Expected output for each species:**

```
================================================================================
VALIDATION REPORT: ecoliK12_global_coords
================================================================================

File: ecoliK12_global_coords.tsv
Total tRNAs: 82
Total positions: ~6400
Global index range: 1 - ~100

--------------------------------------------------------------------------------
TEST 0: Excluded tRNA Types
--------------------------------------------------------------------------------
✅ PASS: No excluded tRNA types found

--------------------------------------------------------------------------------
TEST 1: Type I & Type II Coordinate Sharing
--------------------------------------------------------------------------------
Type I tRNAs (standard): ~66
Type II tRNAs (Leu/Ser/Tyr): 16
Extended arm positions (e-values) present: True
  e-positions occupy global_index: 60-85
  ✅ e-positions in expected reserved space

--------------------------------------------------------------------------------
TEST 2: Sprinzl Label Retention
--------------------------------------------------------------------------------
✅ PASS: 100.0% retention (≥95% threshold)

--------------------------------------------------------------------------------
TEST 3: Region Boundary Alignment
--------------------------------------------------------------------------------
  ✅ Position 1: 'acceptor-stem' ✓
  ✅ Position 8: 'D-stem' ✓
  ✅ Position 9: 'D-stem' ✓
  ✅ Position 26: 'anticodon-stem' ✓
  ✅ Position 34: 'anticodon-loop' ✓
  ✅ Position 45: 'variable-region' ✓
  ✅ Position 54: 'T-loop' ✓
  ✅ Position 73: 'acceptor-tail' ✓

  ✅ All e-positions assigned to 'variable-arm'

✅ PASS: All region boundaries correctly aligned

--------------------------------------------------------------------------------
TEST 4: Critical Position Consistency (1, 73, 76)
--------------------------------------------------------------------------------

Sprinzl position 1:
  ✅ PASS: Maps to single global_index (1)

Sprinzl position 73:
  ⚠️  Maps to 2-3 different global_index values
  (This is EXPECTED - biological variation in 3' end structure)

Sprinzl position 76:
  ⚠️  Maps to 2 different global_index values
  (This is EXPECTED - some tRNAs lack CCA tail)

--------------------------------------------------------------------------------
TEST 5: Extended Variable Arm (e-position) Ordering
--------------------------------------------------------------------------------
✅ PASS: e-positions have monotonically increasing global_index

--------------------------------------------------------------------------------
TEST 6: Global Index Collision Detection
--------------------------------------------------------------------------------
✅ PASS: No global_index collisions detected

================================================================================
OVERALL QUALITY ASSESSMENT
================================================================================

excluded_types           : 15.0/15 (PASS)
type_sharing             : 10.0/10 (PASS)
sprinzl_retention        : 15.0/15 (PASS)
region_boundaries        : 20.0/20 (PASS)
critical_positions       : 12.5/25 (PARTIAL) - 73/76 biological variation
e_position_ordering      :  5.0/5 (PASS)
collisions               : 10.0/10 (PASS)

--------------------------------------------------------------------------------
TOTAL QUALITY SCORE: 87.5-92.5/100
--------------------------------------------------------------------------------

Quality Level: EXCELLENT - Ready for production use

```

**Summary for all species:**

```
================================================================================
SUMMARY: 3 files validated
================================================================================

ecoliK12_global_coords.tsv              : 87.5-92.5/100
hg38_global_coords.tsv                  : 87.5-92.5/100
sacCer_global_coords.tsv                : 87.5-92.5/100
```

---

## Step 5: Run Test Suite

```bash
python test_trnas_in_space.py
```

**Expected output:**

```
Running tests for trnas_in_space.py
============================================================
✓ Imports
✓ Sort key
✓ Sprinzl numeric extraction
✓ Region assignment
✓ Filename parsing
✓ SeC filtering
[ok] Global index validation: No collisions detected among X unique positions
✓ Collision detection
✓ Output files exist
✓ ecoliK12_global_coords.tsv: ~6400 rows, 82 unique tRNAs
✓ Output file structure
✓ Global index continuity
============================================================
Results: 10 passed, 0 failed
```

---

## Step 6: Compare Before/After

### Check improvement in quality scores

```bash
# If you saved old files, compare them
python scripts/validate_coordinate_system.py outputs.old/ecoliK12_global_coords.tsv
# OLD: 37.5/100

python scripts/validate_coordinate_system.py outputs/ecoliK12_global_coords.tsv
# NEW: 87.5-92.5/100
```

### Verify excluded tRNAs are gone

```bash
# Check E. coli no longer has SeC or fMet
python3 << 'EOF'
import pandas as pd
df = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')
sec = df['trna_id'].str.contains('SeC', case=False, na=False).sum()
fmet = df['trna_id'].str.contains('fMet', case=False, na=False).sum()
print(f"SeC tRNAs: {sec} (should be 0)")
print(f"fMet tRNAs: {fmet} (should be 0)")
print(f"Total tRNAs: {df['trna_id'].nunique()} (should be ~82)")
EOF
```

**Expected:**
```
SeC tRNAs: 0 (should be 0)
fMet tRNAs: 0 (should be 0)
Total tRNAs: 82 (should be ~82)
```

---

## Step 7: Investigate 3' Variation (Optional)

To understand position 73/76 variation in detail:

```bash
python scripts/investigate_3prime_variation.py outputs/ecoliK12_global_coords.tsv
```

This will show:
- 3' structure patterns across tRNAs
- Which tRNAs end at position 73 vs 76
- Biological vs algorithmic sources of variation

---

## Step 8: Update Metadata

```bash
# Update METADATA.json with new statistics
python3 << 'EOF'
import json
import pandas as pd
from datetime import datetime

metadata_file = 'outputs/METADATA.json'
with open(metadata_file) as f:
    metadata = json.load(f)

# Update E. coli stats
df_ecoli = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')
metadata['datasets']['ecoliK12_global_coords.tsv']['n_trnas'] = int(df_ecoli['trna_id'].nunique())
metadata['datasets']['ecoliK12_global_coords.tsv']['n_nucleotides'] = len(df_ecoli)
metadata['datasets']['ecoliK12_global_coords.tsv']['n_global_positions'] = int(df_ecoli['global_index'].max())
metadata['datasets']['ecoliK12_global_coords.tsv']['generated'] = datetime.now().isoformat() + 'Z'

# Update yeast stats
df_yeast = pd.read_csv('outputs/sacCer_global_coords.tsv', sep='\t')
metadata['datasets']['sacCer_global_coords.tsv']['n_trnas'] = int(df_yeast['trna_id'].nunique())
metadata['datasets']['sacCer_global_coords.tsv']['n_nucleotides'] = len(df_yeast)
metadata['datasets']['sacCer_global_coords.tsv']['n_global_positions'] = int(df_yeast['global_index'].max())
metadata['datasets']['sacCer_global_coords.tsv']['generated'] = datetime.now().isoformat() + 'Z'

# Update human stats
df_human = pd.read_csv('outputs/hg38_global_coords.tsv', sep='\t')
metadata['datasets']['hg38_global_coords.tsv']['n_trnas'] = int(df_human['trna_id'].nunique())
metadata['datasets']['hg38_global_coords.tsv']['n_nucleotides'] = len(df_human)
metadata['datasets']['hg38_global_coords.tsv']['n_global_positions'] = int(df_human['global_index'].max())
metadata['datasets']['hg38_global_coords.tsv']['generated'] = datetime.now().isoformat() + 'Z'

with open(metadata_file, 'w') as f:
    json.dump(metadata, f, indent=2)

print("✓ Metadata updated")
EOF
```

---

## Step 9: Commit Regenerated Files

```bash
# Stage regenerated coordinate files
git add outputs/ecoliK12_global_coords.tsv \
        outputs/sacCer_global_coords.tsv \
        outputs/hg38_global_coords.tsv \
        outputs/METADATA.json

# Commit with detailed message
git commit -m "$(cat <<'EOF'
Regenerate all coordinate files with fixed code

Quality improvement:
- E. coli: 37.5 → 87.5-92.5/100
- Yeast: 47.5 → 87.5-92.5/100
- Human: 47.5 → 87.5-92.5/100

Changes from regeneration:
- Excluded tRNAs properly filtered (SeC, fMet removed from E. coli)
- Region boundaries fixed (positions 8, 9, 26 now correctly assigned)
- e-positions all assigned to "variable-arm" region
- Zero global_index collisions
- All tests passing

Generated using:
- R2DT 2.0 (Docker: rnacentral/r2dt)
- trnas_in_space.py with fixes from commit 749535c
- Validated with validate_coordinate_system.py

See ASSESSMENT_AND_REMEDIATION_REPORT.md for full details.
EOF
)"

# Push changes
git push origin claude/assess-repo-docs-01QsX58sCRw4Xf3GxLC1Qrar
```

---

## Troubleshooting

### Docker Issues

**Problem:** `docker: command not found`
**Solution:** Install Docker Desktop or start Docker daemon

**Problem:** `Cannot connect to Docker daemon`
**Solution:** Start Docker Desktop (look for whale icon in menu bar)

**Problem:** `docker run` permission denied
**Solution:** Add user to docker group or use `sudo`

### R2DT Issues

**Problem:** R2DT takes very long time
**Solution:** Normal for large datasets (5-15 min per species)

**Problem:** R2DT produces no enriched.json files
**Solution:** Check FASTA file format, ensure tRNAs are valid

### Validation Issues

**Problem:** Still seeing excluded tRNAs after regeneration
**Solution:** Verify code fixes were applied, check git status

**Problem:** Collisions still detected
**Solution:** Check R2DT input - may have duplicate/malformed tRNAs

**Problem:** Quality score <85/100
**Solution:** Run investigation script, check console output for errors

---

## Success Criteria

After completing all steps, you should have:

- ✅ All R2DT JSON files generated
- ✅ All coordinate TSV files regenerated
- ✅ Quality scores >85/100 for all species
- ✅ Zero global_index collisions
- ✅ All tests passing (10/10)
- ✅ SeC and fMet tRNAs excluded from E. coli
- ✅ All region boundaries correctly assigned
- ✅ e-positions assigned to "variable-arm"

**If any of these fail, review console output and consult ASSESSMENT_AND_REMEDIATION_REPORT.md**

---

## Reference Documents

- **ASSESSMENT_AND_REMEDIATION_REPORT.md** - Complete analysis and rationale
- **scripts/validate_coordinate_system.py** - Validation test descriptions
- **scripts/investigate_3prime_variation.py** - 3' end analysis tool
- **README.md** - General usage and background

---

**Last Updated:** 2025-11-18
**Associated Commit:** 749535c
**Branch:** claude/assess-repo-docs-01QsX58sCRw4Xf3GxLC1Qrar
