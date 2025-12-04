# Analysis Guidelines for Global tRNA Coordinates

This guide provides practical recommendations for researchers using the tRNAs-in-space global coordinate system. Different analysis types have varying levels of confidence and validation requirements.

## Choosing Coordinate Files

Each organism has a **single unified coordinate file** containing all tRNA types:

| Organism | Coordinate File | tRNAs | Positions |
|----------|-----------------|-------|-----------|
| E. coli K12 | `ecoliK12_global_coords.tsv` | 82 | 122 |
| S. cerevisiae | `sacCer_global_coords.tsv` | 267 | 115 |
| H. sapiens | `hg38_global_coords.tsv` | 416 | 122 |

Both Type I (standard) and Type II (Leu, Ser, Tyr with extended variable arms) tRNAs are included in the same file. Mitochondrial tRNAs are in separate files generated with `--mito` flag.

See [README.md](README.md#coordinate-file-organization) for details on file organization.

## Analysis Confidence Levels

### ✅ High-Confidence Analyses

These analyses are well-supported by the coordinate system design and have minimal limitations.

#### Cross-Amino-Acid Comparisons
**What**: Compare modification patterns, conservation, or experimental signals across different amino acid families.

**Examples**:
- Pseudouridine frequency: Alanine vs. Leucine vs. Phenylalanine tRNAs
- Structural conservation: Acceptor stem patterns across all amino acids
- Cross-species modification patterns by amino acid family

**Why High-Confidence**:
- Major structural differences averaged out across amino acid families
- Large sample sizes reduce impact of individual tRNA variations
- Biologically meaningful patterns emerge at this level

**Code Example**:
```python
# Cross-amino-acid modification analysis
aa_comparison = df.groupby(['amino_acid', 'global_index']).agg({
    'modification_frequency': 'mean',
    'conservation_score': 'mean'
})

# Statistical comparison
from scipy import stats
ala_mods = df[df['amino_acid'] == 'Ala']['modification_frequency']
leu_mods = df[df['amino_acid'] == 'Leu']['modification_frequency']
t_stat, p_value = stats.ttest_ind(ala_mods, leu_mods)
```

#### Structural Domain Analysis
**What**: Analyze patterns within specific tRNA structural regions.

**Supported Domains**:
- **Acceptor stem** (positions 1-7, 66-72): Amino acid attachment
- **D-stem and loop** (positions 8-25): Tertiary interactions
- **Anticodon region** (positions 26-44): Codon recognition
- **T-stem and loop** (positions 49-65): Ribosome binding

**Examples**:
- Modification hotspots in anticodon loops
- Conservation patterns in acceptor stems
- Structural variation in D-loops

**Code Example**:
```python
# Domain-specific analysis
domain_patterns = df.groupby(['region', 'global_index']).agg({
    'modification_rate': 'mean',
    'sequence_conservation': 'mean'
}).reset_index()

# Compare domains
acceptor_data = domain_patterns[domain_patterns['region'] == 'acceptor-stem']
anticodon_data = domain_patterns[domain_patterns['region'] == 'anticodon-loop']
```

#### Type I vs Type II Comparative Studies
**What**: Compare standard tRNAs with extended variable arm tRNAs.

**Applications**:
- Extended arm impact on surrounding regions
- Structural compensation mechanisms
- Evolution of tRNA complexity

**File Selection**: Type I and Type II tRNAs are in the same unified coordinate file. Filter by amino acid to separate them (Leu, Ser, Tyr are Type II; all others are Type I).

**Code Example**:
```python
# Load unified coordinate file
import pandas as pd

df = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')

# Extract amino acid from trna_id (e.g., "tRNA-Ala-AGC-1-1" -> "Ala")
df['amino_acid'] = df['trna_id'].str.extract(r'tRNA-([A-Za-z]+)-')

# Separate Type I and Type II
type2_aas = ['Leu', 'Ser', 'Tyr']
type1_df = df[~df['amino_acid'].isin(type2_aas)]
type2_df = df[df['amino_acid'].isin(type2_aas)]

# Compare T-loop regions
t_loop_type1 = type1_df[type1_df['region'] == 'T-loop'].groupby('sprinzl_label')['residue'].value_counts()
t_loop_type2 = type2_df[type2_df['region'] == 'T-loop'].groupby('sprinzl_label')['residue'].value_counts()
```

### ⚠️ Moderate-Confidence Analyses

These analyses are supported but require additional validation for robust conclusions.

#### Position-Specific Modification Studies
**What**: Focus on specific structural positions across multiple tRNAs.

**Requirements**:
- Validate critical positions against literature
- Use multiple validation approaches
- Consider biological significance vs. technical precision

**When to Use Caution**:
- Single-position analyses without broader context
- Claims about specific position numbers without structural validation
- Fine-grained spatial relationships between adjacent positions

**Validation Approach**:
```python
# Position-specific analysis with validation
position_data = df[df['global_index'] == 47]  # Example position

# Validate against known biology
known_modifications = {
    'Phe': 'pseudouridine',  # Literature reference
    'Leu': 'methylation',    # Literature reference
}

# Check consistency
observed_mods = position_data.groupby('amino_acid')['modification'].first()
validation_results = {}
for aa, expected in known_modifications.items():
    observed = observed_mods.get(aa, None)
    validation_results[aa] = observed == expected
```

#### Inter-Species Comparative Analysis
**What**: Compare tRNA features between species using global coordinates.

**Considerations**:
- Species-specific annotation variations
- Different R2DT template matches
- Evolutionary distance effects on structural conservation

**Best Practices**:
```python
# Cross-species analysis with validation
species_comparison = pd.merge(
    human_data.groupby(['amino_acid', 'global_index'])['feature'].mean(),
    ecoli_data.groupby(['amino_acid', 'global_index'])['feature'].mean(),
    left_index=True, right_index=True, suffixes=('_human', '_ecoli')
)

# Focus on high-confidence regions
reliable_positions = species_comparison[
    (species_comparison['global_index'] < 65) |  # Before extended arm
    (species_comparison['global_index'] > 85)    # After extended arm
]
```

### 🔴 Use-With-Caution Scenarios

These analyses may be affected by coordinate system limitations and require careful validation.

#### Fine-Grained Isoacceptor Comparisons
**What**: Detailed comparisons between tRNAs of the same amino acid family.

**Limitations**:
- Minor coordinate offsets between isoacceptors (e.g., Ala-AGC vs. Ala-TGC)
- Template matching variations affecting position assignments
- Small sample sizes amplify technical artifacts

**Alternative Approaches**:
```python
# Individual tRNA family analysis (alternative to global coordinates)
ala_trnas = df[df['amino_acid'] == 'Ala']

# Use original Sprinzl positions for within-family analysis
ala_analysis = ala_trnas.pivot_table(
    index='trna_id',
    columns='sprinzl_label',  # Use original labels, not global_index
    values='signal'
)
```

#### Extended Variable Arm Detailed Analysis
**What**: Precise analysis of Type II extended arm internal structure.

**Considerations**:
- Extended arm positions (e1-e24) are mapped to reserved space
- Internal structure of extended arms varies between amino acids
- Template matching may not capture all extended arm variations

**Validation Requirements**:
```python
# Extended arm analysis with validation
extended_arm_data = df[df['sprinzl_label'].str.match(r'^e\d+$', na=False)]

# Validate against known Type II structures
expected_type2_aas = ['Leu', 'Ser', 'Tyr']
observed_type2_aas = extended_arm_data['amino_acid'].unique()

print(f"Expected Type II: {expected_type2_aas}")
print(f"Observed Type II: {observed_type2_aas}")
print(f"Validation: {set(observed_type2_aas).issubset(set(expected_type2_aas))}")
```

## Quality Assessment Guidelines

### Data Validation Checklist

Before conducting analysis, verify:

1. **tRNA Type Distribution**:
```python
# Check tRNA type representation
type_counts = df.groupby(['amino_acid', 'trna_type'])['trna_id'].nunique()
print("tRNA type distribution:")
print(type_counts)
```

2. **Coordinate Range Coverage**:
```python
# Verify coordinate coverage
coord_coverage = df.groupby('global_index')['trna_id'].nunique()
low_coverage = coord_coverage[coord_coverage < 5]  # Positions with <5 tRNAs
print(f"Positions with low coverage: {len(low_coverage)}")
```

3. **Excluded tRNA Verification**:
```python
# Confirm excluded types are not present
excluded_check = df['trna_id'].str.contains('SeC|Sec|mito|iMet', case=False).any()
print(f"Excluded tRNAs present: {excluded_check}")  # Should be False
```

### Statistical Power Considerations

**Minimum Sample Sizes**:
- **Cross-amino-acid comparisons**: ≥3 amino acid families
- **Position-specific analysis**: ≥5 tRNAs per position
- **Type I vs Type II**: ≥3 representative amino acids per type

**Multiple Testing Correction**:
```python
from statsmodels.stats.multitest import multipletests

# Correct for multiple position testing
p_values = [...]  # Your position-wise p-values
rejected, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')
```

## Common Analysis Patterns

### Pattern 1: Modification Frequency Heatmap
```python
import seaborn as sns
import matplotlib.pyplot as plt

# Create modification frequency matrix
mod_matrix = df.pivot_table(
    index='amino_acid',
    columns='global_index',
    values='modification_frequency',
    aggfunc='mean'
)

# Generate heatmap
plt.figure(figsize=(20, 8))
sns.heatmap(mod_matrix, cmap='viridis', cbar_kws={'label': 'Modification Frequency'})
plt.title('tRNA Modification Patterns Across Global Coordinates')
plt.xlabel('Global Coordinate Position')
plt.ylabel('Amino Acid')
plt.tight_layout()
```

### Pattern 2: Domain Conservation Analysis
```python
# Calculate per-domain conservation
domain_stats = df.groupby(['region', 'amino_acid']).agg({
    'conservation_score': ['mean', 'std', 'count'],
    'global_index': ['min', 'max']
}).round(3)

# Identify most/least conserved domains
domain_conservation = df.groupby('region')['conservation_score'].mean().sort_values(ascending=False)
print("Domain conservation ranking:")
print(domain_conservation)
```

### Pattern 3: Type Comparison with Statistics
```python
import scipy.stats as stats

# Statistical comparison between Type I and II
type1_positions = df[df['tRNA_type'] == 'Type_I']['global_index'].values
type2_positions = df[df['tRNA_type'] == 'Type_II']['global_index'].values

# Kolmogorov-Smirnov test for distribution differences
ks_stat, ks_p = stats.ks_2samp(type1_positions, type2_positions)
print(f"Type I vs Type II distribution difference: KS={ks_stat:.3f}, p={ks_p:.3e}")
```

## Troubleshooting Common Issues

### Issue: Unexpected Position Gaps
**Symptoms**: Missing data at expected coordinate positions
**Diagnosis**: Check if positions fall in extended arm region (65-85)
**Solution**: Verify tRNA types in your dataset

### Issue: Low Statistical Power
**Symptoms**: High p-values, inconsistent results
**Diagnosis**: Check sample sizes per group
**Solution**: Combine amino acids or use broader categories

### Issue: Inconsistent Cross-Species Results
**Symptoms**: Different patterns between species
**Diagnosis**: Species-specific annotation differences
**Solution**: Focus on well-annotated positions, validate with literature

## Reporting Best Practices

### Methods Section Template
```
"tRNA comparative analysis was performed using the tRNAs-in-space global coordinate
system (White et al., 2025), which provides standardized position mapping for nuclear
elongator tRNAs. Analysis focused on [Type I/Type II/both] tRNAs from [species],
with [N] tRNAs representing [X] amino acid families. Position-specific analyses
were validated against known structural features, and multiple testing correction
was applied using the Benjamini-Hochberg method."
```

### Results Interpretation Guidelines
- Report confidence level (high/moderate/caution) for each analysis type
- Validate key findings against literature when possible
- Acknowledge coordinate system limitations for fine-grained analyses
- Provide biological context for statistical significance

---

For detailed system scope and limitations, see [docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SCOPE.md](docs/archive/coordinate-fixes/COORDINATE_SYSTEM_SCOPE.md).

For validation methods and quality assessment, see [docs/archive/coordinate-fixes/COORDINATE_SYSTEM_VALIDATION.md](docs/archive/coordinate-fixes/COORDINATE_SYSTEM_VALIDATION.md).