# Coordinate System Validation: Research Quality Assessment

This guide provides methods for researchers to validate the quality and appropriateness of global tRNA coordinates for their specific research questions.

## Pre-Analysis Validation

Before conducting research with global tRNA coordinates, validate that your dataset and analysis approach are appropriate for the coordinate system's scope.

### Dataset Composition Check

#### 1. Verify tRNA Type Distribution
```python
import pandas as pd

def validate_trna_composition(df):
    """Check tRNA type distribution in your dataset."""

    # Check amino acid representation
    aa_counts = df.groupby('amino_acid')['trna_id'].nunique().sort_values(ascending=False)

    # Identify Type I vs Type II amino acids
    type1_aas = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Lys', 'Met', 'Phe', 'Pro', 'Thr', 'Trp', 'Val']
    type2_aas = ['Leu', 'Ser', 'Tyr']

    type1_present = df['amino_acid'].isin(type1_aas).any()
    type2_present = df['amino_acid'].isin(type2_aas).any()

    print("=== Dataset Composition Validation ===")
    print(f"Total amino acid families: {len(aa_counts)}")
    print(f"Total unique tRNAs: {df['trna_id'].nunique()}")
    print(f"Type I amino acids present: {type1_present}")
    print(f"Type II amino acids present: {type2_present}")

    if type1_present and type2_present:
        print("✅ Dataset suitable for Type I/II comparative analysis")
    elif type1_present:
        print("⚠️  Dataset contains only Type I tRNAs - limited to standard structure analysis")
    elif type2_present:
        print("⚠️  Dataset contains only Type II tRNAs - extended arm analysis only")

    return {
        'aa_counts': aa_counts,
        'type1_present': type1_present,
        'type2_present': type2_present
    }

# Example usage
df = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')
composition = validate_trna_composition(df)
```

#### 2. Check for Excluded tRNA Types
```python
def check_excluded_types(df):
    """Verify that excluded tRNA types are not present."""

    excluded_patterns = ['SeC', 'Sec', 'mito', 'iMet', 'initiator']
    excluded_found = {}

    for pattern in excluded_patterns:
        found = df['trna_id'].str.contains(pattern, case=False, na=False).sum()
        excluded_found[pattern] = found
        if found > 0:
            print(f"⚠️  Found {found} {pattern} tRNAs - these should be excluded")
            excluded_examples = df[df['trna_id'].str.contains(pattern, case=False, na=False)]['trna_id'].unique()[:3]
            print(f"   Examples: {list(excluded_examples)}")

    total_excluded = sum(excluded_found.values())
    if total_excluded == 0:
        print("✅ No excluded tRNA types detected")
    else:
        print(f"❌ Total excluded types found: {total_excluded}")
        print("   Consider reprocessing data with updated filtering")

    return excluded_found

excluded_check = check_excluded_types(df)
```

#### 3. Coordinate Coverage Assessment
```python
def assess_coordinate_coverage(df):
    """Evaluate coverage across global coordinate positions."""

    # Count tRNAs per position
    coverage = df.groupby('global_index')['trna_id'].nunique().sort_index()

    # Identify regions
    standard_region = coverage[(coverage.index < 65) | (coverage.index > 85)]
    extended_arm_region = coverage[(coverage.index >= 65) & (coverage.index <= 85)]

    # Coverage statistics
    min_coverage = coverage.min()
    max_coverage = coverage.max()
    mean_coverage = coverage.mean()

    print("=== Coordinate Coverage Assessment ===")
    print(f"Position range: {coverage.index.min()} - {coverage.index.max()}")
    print(f"Coverage range: {min_coverage} - {max_coverage} tRNAs per position")
    print(f"Mean coverage: {mean_coverage:.1f} tRNAs per position")

    # Low coverage positions
    low_coverage_threshold = 3
    low_coverage_positions = coverage[coverage < low_coverage_threshold]
    if len(low_coverage_positions) > 0:
        print(f"⚠️  {len(low_coverage_positions)} positions with <{low_coverage_threshold} tRNAs:")
        print(f"   Positions: {list(low_coverage_positions.index)}")

    # Extended arm coverage
    if len(extended_arm_region) > 0:
        print(f"✅ Extended arm region (65-85): {extended_arm_region.sum()} total tRNA-positions")
    else:
        print("⚠️  No extended arm coverage - dataset may lack Type II tRNAs")

    return coverage

coverage = assess_coordinate_coverage(df)
```

## Analysis-Specific Validation

### Position-Specific Analysis Validation

```python
def validate_position_analysis(df, positions_of_interest):
    """Validate specific positions for detailed analysis."""

    validation_results = {}

    for pos in positions_of_interest:
        pos_data = df[df['global_index'] == pos]

        if len(pos_data) == 0:
            validation_results[pos] = {'status': 'missing', 'message': 'Position not found in dataset'}
            continue

        # Check amino acid representation
        aa_counts = pos_data['amino_acid'].value_counts()
        unique_aas = len(aa_counts)

        # Check structural region
        regions = pos_data['region'].value_counts()
        primary_region = regions.index[0] if len(regions) > 0 else 'unknown'

        # Check for extended arm positions
        is_extended_arm = (pos >= 65) and (pos <= 85)

        validation_results[pos] = {
            'status': 'valid',
            'amino_acids': unique_aas,
            'total_trnas': len(pos_data),
            'primary_region': primary_region,
            'is_extended_arm': is_extended_arm,
            'aa_distribution': dict(aa_counts)
        }

        print(f"Position {pos}:")
        print(f"  Region: {primary_region}")
        print(f"  tRNAs: {len(pos_data)} across {unique_aas} amino acids")
        if is_extended_arm:
            print(f"  ⚠️  Extended arm position - Type II tRNAs only")

    return validation_results

# Example: Validate positions of interest
important_positions = [32, 34, 36, 47, 55]  # anticodon, variable, T-loop
position_validation = validate_position_analysis(df, important_positions)
```

### Cross-Species Validation

```python
def validate_cross_species_analysis(df1, df2, species1_name, species2_name):
    """Validate datasets for cross-species comparative analysis."""

    print(f"=== Cross-Species Validation: {species1_name} vs {species2_name} ===")

    # Compare amino acid coverage
    aa1 = set(df1['amino_acid'].unique())
    aa2 = set(df2['amino_acid'].unique())

    common_aas = aa1.intersection(aa2)
    species1_only = aa1 - aa2
    species2_only = aa2 - aa1

    print(f"Common amino acids: {len(common_aas)} ({list(sorted(common_aas))})")
    if species1_only:
        print(f"{species1_name} only: {list(sorted(species1_only))}")
    if species2_only:
        print(f"{species2_name} only: {list(sorted(species2_only))}")

    # Compare coordinate coverage
    coords1 = set(df1['global_index'].unique())
    coords2 = set(df2['global_index'].unique())

    common_coords = coords1.intersection(coords2)
    coverage_similarity = len(common_coords) / len(coords1.union(coords2))

    print(f"Coordinate coverage similarity: {coverage_similarity:.2f}")

    if coverage_similarity > 0.8:
        print("✅ Good coordinate coverage overlap for cross-species analysis")
    else:
        print("⚠️  Limited coordinate overlap - validate species-specific differences")

    return {
        'common_amino_acids': common_aas,
        'coordinate_overlap': coverage_similarity,
        'suitable_for_comparison': coverage_similarity > 0.8
    }

# Example cross-species validation
# df_human = pd.read_csv('outputs/hg38_global_coords.tsv', sep='\t')
# df_ecoli = pd.read_csv('outputs/ecoliK12_global_coords.tsv', sep='\t')
# cross_species_val = validate_cross_species_analysis(df_human, df_ecoli, "Human", "E. coli")
```

## Quality Metrics and Diagnostics

### Structural Consistency Checks

```python
def check_structural_consistency(df):
    """Validate structural region assignments and coordinate consistency."""

    print("=== Structural Consistency Validation ===")

    # Check region boundaries
    region_boundaries = df.groupby('region')['global_index'].agg(['min', 'max']).sort_values('min')

    print("Region coordinate ranges:")
    for region, bounds in region_boundaries.iterrows():
        print(f"  {region}: {bounds['min']} - {bounds['max']}")

    # Check for expected regions
    expected_regions = ['acceptor-stem', 'D-stem', 'D-loop', 'anticodon-stem',
                       'anticodon-loop', 'variable-region', 'T-stem', 'T-loop']

    present_regions = set(df['region'].unique())
    missing_regions = set(expected_regions) - present_regions
    unexpected_regions = present_regions - set(expected_regions) - {'unknown'}

    if missing_regions:
        print(f"⚠️  Missing expected regions: {missing_regions}")
    if unexpected_regions:
        print(f"ℹ️  Additional regions found: {unexpected_regions}")

    # Check extended arm consistency
    extended_arm_data = df[df['sprinzl_label'].str.match(r'^e\d+$', na=False)]
    if len(extended_arm_data) > 0:
        extended_arm_coords = extended_arm_data['global_index'].unique()
        expected_range = set(range(65, 86))  # 65-85
        actual_range = set(extended_arm_coords)

        if actual_range.issubset(expected_range):
            print(f"✅ Extended arm positions properly mapped to reserved space (65-85)")
        else:
            unexpected_extended = actual_range - expected_range
            print(f"⚠️  Extended arm positions outside reserved space: {unexpected_extended}")

    return region_boundaries

structural_check = check_structural_consistency(df)
```

### Data Quality Scoring

```python
def calculate_quality_score(df):
    """Calculate overall data quality score for the coordinate dataset."""

    score_components = {}

    # 1. tRNA type diversity (0-25 points)
    aa_count = df['amino_acid'].nunique()
    type_diversity_score = min(25, aa_count * 2)  # Max 25 for 12+ amino acids
    score_components['type_diversity'] = type_diversity_score

    # 2. Coordinate coverage (0-25 points)
    coverage = df.groupby('global_index')['trna_id'].nunique()
    coverage_consistency = 25 * (1 - coverage.std() / coverage.mean())  # Lower variation = higher score
    coverage_consistency = max(0, min(25, coverage_consistency))
    score_components['coverage_consistency'] = coverage_consistency

    # 3. Structural annotation completeness (0-25 points)
    unknown_fraction = (df['region'] == 'unknown').sum() / len(df)
    annotation_score = 25 * (1 - unknown_fraction)
    score_components['annotation_completeness'] = annotation_score

    # 4. Extended arm handling (0-25 points)
    has_extended_arms = df['sprinzl_label'].str.match(r'^e\d+$', na=False).any()
    if has_extended_arms:
        extended_in_range = df[(df['global_index'] >= 65) & (df['global_index'] <= 85) &
                              df['sprinzl_label'].str.match(r'^e\d+$', na=False)].shape[0]
        total_extended = df[df['sprinzl_label'].str.match(r'^e\d+$', na=False)].shape[0]
        extended_score = 25 * (extended_in_range / total_extended) if total_extended > 0 else 25
    else:
        extended_score = 25  # Full score if no extended arms expected
    score_components['extended_arm_handling'] = extended_score

    # Total score
    total_score = sum(score_components.values())

    print("=== Data Quality Assessment ===")
    for component, score in score_components.items():
        print(f"{component}: {score:.1f}/25")
    print(f"Total Quality Score: {total_score:.1f}/100")

    if total_score >= 80:
        print("✅ High quality dataset suitable for most analyses")
    elif total_score >= 60:
        print("⚠️  Moderate quality dataset - validate specific analyses")
    else:
        print("❌ Low quality dataset - consider data reprocessing")

    return score_components, total_score

quality_assessment = calculate_quality_score(df)
```

## Literature Validation

### Known Position Validation

```python
def validate_against_known_positions(df):
    """Validate coordinate assignments against known tRNA biology."""

    known_positions = {
        'anticodon_positions': [32, 33, 34, 35, 36],  # Anticodon loop
        'variable_region': list(range(45, 49)),        # Type I variable
        'extended_variable': list(range(65, 86)),      # Type II extended (e1-e24)
        'T_loop': list(range(54, 61)),                 # T-loop region
        'acceptor_stem_5prime': list(range(1, 8)),     # 5' acceptor
        'acceptor_stem_3prime': list(range(66, 73))    # 3' acceptor
    }

    validation_results = {}

    for region_name, expected_positions in known_positions.items():
        # Check if positions exist in dataset
        positions_present = df[df['global_index'].isin(expected_positions)]['global_index'].unique()
        coverage = len(positions_present) / len(expected_positions)

        validation_results[region_name] = {
            'expected_positions': len(expected_positions),
            'positions_present': len(positions_present),
            'coverage': coverage
        }

        print(f"{region_name}:")
        print(f"  Expected: {len(expected_positions)} positions")
        print(f"  Present: {len(positions_present)} positions")
        print(f"  Coverage: {coverage:.2f}")

        if coverage >= 0.8:
            print(f"  ✅ Good coverage")
        else:
            print(f"  ⚠️  Limited coverage")

    return validation_results

literature_validation = validate_against_known_positions(df)
```

## Validation Reporting

### Generate Validation Report

```python
def generate_validation_report(df, output_file='coordinate_validation_report.txt'):
    """Generate comprehensive validation report for dataset."""

    with open(output_file, 'w') as f:
        f.write("tRNAs-in-Space Global Coordinate System Validation Report\n")
        f.write("=" * 60 + "\n\n")

        # Dataset summary
        f.write(f"Dataset: {df['trna_id'].nunique()} unique tRNAs\n")
        f.write(f"Amino acids: {df['amino_acid'].nunique()} families\n")
        f.write(f"Global coordinate range: {df['global_index'].min()} - {df['global_index'].max()}\n")
        f.write(f"Total positions: {len(df)}\n\n")

        # Run all validations and capture results
        print("Running validation checks...")

        # Composition check
        f.write("Composition Validation:\n")
        composition = validate_trna_composition(df)
        f.write(f"- Type I amino acids: {composition['type1_present']}\n")
        f.write(f"- Type II amino acids: {composition['type2_present']}\n\n")

        # Quality score
        f.write("Quality Assessment:\n")
        score_components, total_score = calculate_quality_score(df)
        f.write(f"- Overall quality score: {total_score:.1f}/100\n")
        for component, score in score_components.items():
            f.write(f"- {component}: {score:.1f}/25\n")
        f.write("\n")

        # Recommendations
        f.write("Analysis Recommendations:\n")
        if total_score >= 80:
            f.write("- Dataset suitable for all supported analysis types\n")
        elif total_score >= 60:
            f.write("- Validate specific positions for detailed analyses\n")
            f.write("- Consider limitations for fine-grained comparisons\n")
        else:
            f.write("- Consider data reprocessing with updated filtering\n")
            f.write("- Focus on high-coverage positions only\n")

    print(f"Validation report written to: {output_file}")

# Generate comprehensive report
generate_validation_report(df)
```

## Quick Validation Checklist

For rapid quality assessment, use this checklist:

```python
def quick_validation_checklist(df):
    """Quick go/no-go validation for coordinate analysis."""

    checks = []

    # 1. No excluded types
    excluded = df['trna_id'].str.contains('SeC|Sec|mito|iMet', case=False, na=False).any()
    checks.append(("No excluded tRNA types", not excluded))

    # 2. Reasonable tRNA count
    trna_count = df['trna_id'].nunique()
    checks.append(("Sufficient tRNA count (>20)", trna_count > 20))

    # 3. Amino acid diversity
    aa_count = df['amino_acid'].nunique()
    checks.append(("Amino acid diversity (>5)", aa_count > 5))

    # 4. Coordinate range
    coord_range = df['global_index'].max() - df['global_index'].min()
    checks.append(("Reasonable coordinate range (>60)", coord_range > 60))

    # 5. Extended arms present (if applicable)
    has_extended = df['sprinzl_label'].str.match(r'^e\d+$', na=False).any()
    has_type2 = df['amino_acid'].isin(['Leu', 'Ser', 'Tyr']).any()
    extended_check = (not has_type2) or has_extended
    checks.append(("Extended arms properly detected", extended_check))

    print("Quick Validation Checklist:")
    all_passed = True
    for check_name, passed in checks:
        status = "✅" if passed else "❌"
        print(f"{status} {check_name}")
        if not passed:
            all_passed = False

    if all_passed:
        print("\n✅ Dataset ready for analysis")
    else:
        print("\n⚠️  Dataset issues detected - see detailed validation")

    return all_passed

# Quick check
ready_for_analysis = quick_validation_checklist(df)
```

---

This validation framework ensures that researchers can confidently assess whether their datasets and analysis approaches are appropriate for the global tRNA coordinate system's capabilities and limitations.

For analysis guidelines using validated datasets, see [ANALYSIS_GUIDELINES.md](ANALYSIS_GUIDELINES.md).

For system scope and limitations, see [COORDINATE_SYSTEM_SCOPE.md](COORDINATE_SYSTEM_SCOPE.md).