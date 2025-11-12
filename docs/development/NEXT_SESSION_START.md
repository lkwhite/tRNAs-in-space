# Quick Start for Next Session

## 🎯 Goal
Continue Modomics integration - map tRNA modifications to Sprinzl positions

## ✅ What's Done
- Phase 1 parser infrastructure is built and tested
- Directory structure created
- Configuration files in place
- Documentation written

## 📋 What You Need
**Before starting next session, add this file:**
```
docs/development/unmodified_tRNA_all_all_rna_sequences.fasta
```
- Same format as modified version (same headers)
- Sequences contain only A, C, G, U (no modification codes)

## 🚀 Quick Start Commands

### Option 1: Test Phase 1 Parser
```
Test the Modomics Phase 1 parser. The unmodified FASTA is now at:
docs/development/unmodified_tRNA_all_all_rna_sequences.fasta

Run the parser and show me statistics:

python scripts/modomics/parse_modomics.py \
    --modified docs/development/modified_tRNA_all_all_rna_sequences.fasta \
    --unmodified docs/development/unmodified_tRNA_all_all_rna_sequences.fasta \
    --codes docs/development/modomicscodes.csv \
    --output outputs/modomics/modomics_modifications.json
```

### Option 2: Continue to Phase 2
```
Continue Modomics integration from docs/development/MODOMICS_INTEGRATION.md

The unmodified FASTA is ready. Please:
1. Test Phase 1 parser
2. Build Phase 2: sequence alignment and Sprinzl position mapping
3. Start with E. coli as test case

Focus on creating scripts/modomics/align_to_sprinzl.py
```

### Option 3: Just Build Phase 2
```
Build Phase 2 of Modomics integration (alignment to Sprinzl positions).

See plan in: docs/development/MODOMICS_INTEGRATION.md

Create scripts/modomics/align_to_sprinzl.py that:
1. Loads modomics_modifications.json
2. Normalizes species names (use config/modomics.yaml)
3. Aligns to gtRNAdb sequences
4. Maps modification positions to Sprinzl indices
5. Outputs modomics_to_sprinzl_mapping.tsv

Test with E. coli first.
```

## 📚 Key Documentation
- **Main guide:** `docs/development/MODOMICS_INTEGRATION.md`
- **Script docs:** `scripts/modomics/README.md`
- **Config:** `config/modomics.yaml`

## 🔍 Check Setup
```bash
# Verify files exist
ls -la docs/development/*modified*.fasta
ls -la scripts/modomics/
ls -la config/modomics.yaml

# Test decoder (should work now)
python scripts/modomics/modification_codes.py docs/development/modomicscodes.csv
```

## 🌳 Git Branch
```
claude/modomics-species-mapping-011CV4PiWXABfn89mspa1222
```

## 📦 Dependencies Needed for Phase 2
```bash
pip install biopython
```

## 💡 Context for Claude
When starting next session, Claude should:
1. Read `docs/development/MODOMICS_INTEGRATION.md` for full context
2. Check if unmodified FASTA exists
3. Test Phase 1 parser first
4. Proceed to Phase 2 implementation

---
**Last session:** Phase 1 setup complete
**Next step:** Test parser → Build alignment pipeline
