#!/usr/bin/env python3
"""statistical_analysis.py - GenomeScope 2.0 estimates of both platforms: reads the summaries, removes outliers (IQR x 3 on each of
the four estimates, both platforms pooled), pairs the platforms by individual and runs paired Wilcoxon signed-rank tests. The
comparison between sexes and the correlations with depth are not reported in the paper.
Usage: python3 statistical_analysis.py --illumina_dir <dir> --aviti_dir <dir> --sex_file <tsv> --coverage_file <tsv> --output_dir <dir>
"""

import argparse
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')


# =============================================================================
# PARSING
# =============================================================================

def parse_genomescope_summary(file_path):
    """Parse a GenomeScope2 summary.txt file and return a dict of metrics."""
    data = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if 'Homozygous (aa)' in line:
            values = re.findall(r'(\d+\.?\d*%)', line)
            if len(values) >= 2:
                data['homozygous_min'] = float(values[0].replace('%', ''))
                data['homozygous_max'] = float(values[1].replace('%', ''))

        elif 'Heterozygous (ab)' in line:
            values = re.findall(r'(\d+\.?\d*%)', line)
            if len(values) >= 2:
                data['heterozygous_min'] = float(values[0].replace('%', ''))
                data['heterozygous_max'] = float(values[1].replace('%', ''))

        elif 'Genome Haploid Length' in line:
            values = re.findall(r'([\d,]+) bp', line)
            if len(values) >= 2:
                data['genome_length_min'] = int(values[0].replace(',', ''))
                data['genome_length_max'] = int(values[1].replace(',', ''))

        elif 'Read Error Rate' in line:
            values = re.findall(r'(\d+\.?\d*%)', line)
            if len(values) >= 1:
                data['error_rate'] = float(values[0].replace('%', ''))

        elif 'Model Fit' in line:
            values = re.findall(r'(\d+\.?\d*%)', line)
            if len(values) >= 2:
                data['model_fit_min'] = float(values[0].replace('%', ''))
                data['model_fit_max'] = float(values[1].replace('%', ''))

    return data


def collect_genomescope_data(illumina_dir, aviti_dir):
    """Collect GenomeScope data from both platform directories."""
    all_data = []

    # Process Illumina
    if os.path.isdir(illumina_dir):
        print(f"Scanning Illumina directory: {illumina_dir}")
        for sample_dir in sorted(os.listdir(illumina_dir)):
            sample_path = os.path.join(illumina_dir, sample_dir)
            if not os.path.isdir(sample_path):
                continue
            # Extract the sample number from the folder name
            sample_num = sample_dir.split('ICDG')[0] if 'ICDG' in sample_dir and sample_dir[0].isdigit() else sample_dir.replace('ICDG', '')
            summary_file = os.path.join(sample_path, "genomescope2_output", f"{sample_dir}_summary.txt")
            if not os.path.exists(summary_file):
                # Try alternative naming
                gs_dir = os.path.join(sample_path, "genomescope2_output")
                if os.path.isdir(gs_dir):
                    candidates = [f for f in os.listdir(gs_dir) if f.endswith('_summary.txt')]
                    if candidates:
                        summary_file = os.path.join(gs_dir, candidates[0])
                    else:
                        continue
                else:
                    continue
            try:
                data = parse_genomescope_summary(summary_file)
                data['sample'] = sample_dir
                data['sample_num'] = sample_num
                data['platform'] = 'Illumina'
                all_data.append(data)
                print(f"  read {sample_dir}")
            except Exception as e:
                print(f"  failed {sample_dir}: {e}")

    # Process Aviti
    if os.path.isdir(aviti_dir):
        print(f"Scanning Aviti directory: {aviti_dir}")
        for sample_dir in sorted(os.listdir(aviti_dir)):
            sample_path = os.path.join(aviti_dir, sample_dir)
            if not os.path.isdir(sample_path):
                continue
            sample_num = sample_dir.replace('ICDG', '') if sample_dir.startswith('ICDG') else sample_dir.split('ICDG')[0]
            summary_file = os.path.join(sample_path, "genomescope2_output", f"{sample_dir}_summary.txt")
            if not os.path.exists(summary_file):
                gs_dir = os.path.join(sample_path, "genomescope2_output")
                if os.path.isdir(gs_dir):
                    candidates = [f for f in os.listdir(gs_dir) if f.endswith('_summary.txt')]
                    if candidates:
                        summary_file = os.path.join(gs_dir, candidates[0])
                    else:
                        continue
                else:
                    continue
            try:
                data = parse_genomescope_summary(summary_file)
                data['sample'] = sample_dir
                data['sample_num'] = sample_num
                data['platform'] = 'Aviti'
                all_data.append(data)
                print(f"  read {sample_dir}")
            except Exception as e:
                print(f"  failed {sample_dir}: {e}")

    df = pd.DataFrame(all_data)
    if len(df) > 0:
        # Compute mean values
        df['heterozygosity'] = (df['heterozygous_min'] + df['heterozygous_max']) / 2
        df['genome_length'] = (df['genome_length_min'] + df['genome_length_max']) / 2
        df['genome_length_gb'] = df['genome_length'] / 1e9
        df['model_fit'] = (df['model_fit_min'] + df['model_fit_max']) / 2
    return df


# =============================================================================
# PAIRS OF READ SETS
# =============================================================================

# read-set number on the other platform, where it differs from the number on this platform
OTHER_PLATFORM_NUM = {'12': '21', '21': '12'}


def get_paired_sample_num(sample_num, platform):
    """Get the matching sample_num on the OTHER platform."""
    if platform == 'Illumina':
        return OTHER_PLATFORM_NUM.get(sample_num, sample_num)
    else:
        return OTHER_PLATFORM_NUM.get(sample_num, sample_num)


def build_paired_df(df):
    """Build a paired dataframe where each row is one individual with both platforms."""
    illumina = df[df['platform'] == 'Illumina'].copy()
    aviti = df[df['platform'] == 'Aviti'].copy()

    # For Illumina samples, determine which Aviti sample_num they pair with
    illumina['pair_aviti_num'] = illumina['sample_num'].apply(
        lambda x: OTHER_PLATFORM_NUM.get(x, x)
    )

    pairs = []
    for _, ill_row in illumina.iterrows():
        aviti_num = ill_row['pair_aviti_num']
        aviti_match = aviti[aviti['sample_num'] == aviti_num]
        if len(aviti_match) > 0:
            av_row = aviti_match.iloc[0]
            pairs.append({
                'individual': ill_row['sample_num'],
                'illumina_sample': ill_row['sample'],
                'aviti_sample': av_row['sample'],
                'illumina_heterozygosity': ill_row['heterozygosity'],
                'aviti_heterozygosity': av_row['heterozygosity'],
                'illumina_error_rate': ill_row['error_rate'],
                'aviti_error_rate': av_row['error_rate'],
                'illumina_model_fit': ill_row['model_fit'],
                'aviti_model_fit': av_row['model_fit'],
                'illumina_genome_length_gb': ill_row['genome_length_gb'],
                'aviti_genome_length_gb': av_row['genome_length_gb'],
            })

    return pd.DataFrame(pairs)


# =============================================================================
# OUTLIER REMOVAL
# =============================================================================

def remove_outliers_iqr(df, columns, multiplier=3.0):
    """Remove outliers using IQR x multiplier criterion. Returns mask of inliers."""
    mask = pd.Series(True, index=df.index)
    outlier_info = []

    for col in columns:
        if col not in df.columns:
            continue
        q1 = df[col].quantile(0.25)
        q3 = df[col].quantile(0.75)
        iqr = q3 - q1
        lower = q1 - multiplier * iqr
        upper = q3 + multiplier * iqr
        col_outliers = (df[col] < lower) | (df[col] > upper)
        n_outliers = col_outliers.sum()
        if n_outliers > 0:
            outlier_samples = df.loc[col_outliers, 'sample'].tolist() if 'sample' in df.columns else df.index[col_outliers].tolist()
            outlier_info.append(f"  {col}: {n_outliers} outlier(s) - {outlier_samples}")
        mask &= ~col_outliers

    return mask, outlier_info


# =============================================================================
# STATISTICAL TESTS
# =============================================================================

def paired_wilcoxon_tests(paired_df):
    """Perform paired Wilcoxon signed-rank tests for each metric."""
    metrics = [
        ('heterozygosity', 'illumina_heterozygosity', 'aviti_heterozygosity'),
        ('error_rate', 'illumina_error_rate', 'aviti_error_rate'),
        ('model_fit', 'illumina_model_fit', 'aviti_model_fit'),
        ('genome_length_gb', 'illumina_genome_length_gb', 'aviti_genome_length_gb'),
    ]

    results = []
    for name, ill_col, av_col in metrics:
        ill_vals = paired_df[ill_col].dropna()
        av_vals = paired_df[av_col].dropna()
        # Align
        valid = paired_df[[ill_col, av_col]].dropna()
        if len(valid) < 3:
            results.append({
                'metric': name,
                'n_pairs': len(valid),
                'illumina_median': np.nan,
                'aviti_median': np.nan,
                'statistic': np.nan,
                'p_value': np.nan,
                'effect_size_r': np.nan,
                'significant': 'insufficient data'
            })
            continue

        ill_v = valid[ill_col].values
        av_v = valid[av_col].values

        try:
            stat, p = stats.wilcoxon(ill_v, av_v, alternative='two-sided')
            # effect_size_r = Z / sqrt(N); the paper gives the matched-pairs rank-biserial r of genomescope_stats.py
            n = len(valid)
            # Approximate Z from p-value for reporting
            z = stats.norm.ppf(1 - p / 2) if p > 0 else np.inf
            r = z / np.sqrt(n) if n > 0 else np.nan
        except Exception:
            stat, p, r = np.nan, np.nan, np.nan

        results.append({
            'metric': name,
            'n_pairs': len(valid),
            'illumina_median': np.median(ill_v),
            'illumina_mean': np.mean(ill_v),
            'illumina_std': np.std(ill_v, ddof=1),
            'aviti_median': np.median(av_v),
            'aviti_mean': np.mean(av_v),
            'aviti_std': np.std(av_v, ddof=1),
            'statistic': stat,
            'p_value': p,
            'effect_size_r': r,
            'significant': 'Yes' if p < 0.05 else 'No'
        })

    return pd.DataFrame(results)


def sex_comparison(df, sex_df):
    """Compare GenomeScope metrics between sexes using Mann-Whitney U tests."""
    # Merge sex info
    # sex_df should have columns: sample_id, sex
    # Map sample_num to sex
    sex_map = dict(zip(sex_df['sample_id'].str.replace('ICDG', ''), sex_df['sex']))

    df = df.copy()
    df['sex'] = df['sample_num'].map(sex_map)

    results = []
    metrics = ['heterozygosity', 'error_rate', 'model_fit', 'genome_length_gb']

    for platform in df['platform'].unique():
        pdata = df[df['platform'] == platform]
        for metric in metrics:
            if metric not in pdata.columns:
                continue
            males = pdata[pdata['sex'] == 'M'][metric].dropna()
            females = pdata[pdata['sex'] == 'F'][metric].dropna()

            if len(males) < 2 or len(females) < 2:
                results.append({
                    'platform': platform,
                    'metric': metric,
                    'n_males': len(males),
                    'n_females': len(females),
                    'male_median': males.median() if len(males) > 0 else np.nan,
                    'female_median': females.median() if len(females) > 0 else np.nan,
                    'statistic': np.nan,
                    'p_value': np.nan,
                    'significant': 'insufficient data'
                })
                continue

            stat, p = stats.mannwhitneyu(males, females, alternative='two-sided')
            # Effect size r
            n = len(males) + len(females)
            z = stats.norm.ppf(1 - p / 2) if p > 0 else np.inf
            r = z / np.sqrt(n) if n > 0 else np.nan

            results.append({
                'platform': platform,
                'metric': metric,
                'n_males': len(males),
                'n_females': len(females),
                'male_median': males.median(),
                'male_mean': males.mean(),
                'male_std': males.std(),
                'female_median': females.median(),
                'female_mean': females.mean(),
                'female_std': females.std(),
                'statistic': stat,
                'p_value': p,
                'effect_size_r': r,
                'significant': 'Yes' if p < 0.05 else 'No'
            })

    return pd.DataFrame(results), df


def correlation_analysis(df, coverage_df):
    """Investigate correlations between coverage and GenomeScope metrics."""
    # Merge coverage into df
    # coverage_df has columns: platform, sample, coverage, cv
    # We need to match sample_num

    df = df.copy()

    # Build coverage map: (platform, sample_num) -> coverage
    cov_map = {}
    if coverage_df is not None and len(coverage_df) > 0:
        for _, row in coverage_df.iterrows():
            plat = row['platform']
            sample = str(row['sample'])
            cov = row['coverage']
            # Extract number
            num = re.sub(r'[^\d]', '', sample.split('ICDG')[0] if 'ICDG' in sample and sample[0].isdigit() else sample.replace('ICDG', ''))
            if plat == 'Illumina':
                cov_map[('Illumina', num)] = cov
            elif plat in ('AVITI', 'Aviti'):
                cov_map[('Aviti', num)] = cov

    df['coverage'] = df.apply(
        lambda row: cov_map.get((row['platform'], row['sample_num']), np.nan), axis=1
    )

    metrics = ['heterozygosity', 'error_rate', 'model_fit', 'genome_length_gb']
    results = []

    for platform in df['platform'].unique():
        pdata = df[df['platform'] == platform].dropna(subset=['coverage'])
        for metric in metrics:
            if metric not in pdata.columns:
                continue
            valid = pdata[['coverage', metric]].dropna()
            if len(valid) < 4:
                results.append({
                    'platform': platform,
                    'metric': metric,
                    'n': len(valid),
                    'spearman_rho': np.nan,
                    'spearman_p': np.nan,
                    'pearson_r': np.nan,
                    'pearson_p': np.nan,
                    'significant_spearman': 'insufficient data'
                })
                continue

            rho, p_rho = stats.spearmanr(valid['coverage'], valid[metric])
            r, p_r = stats.pearsonr(valid['coverage'], valid[metric])

            results.append({
                'platform': platform,
                'metric': metric,
                'n': len(valid),
                'spearman_rho': rho,
                'spearman_p': p_rho,
                'pearson_r': r,
                'pearson_p': p_r,
                'significant_spearman': 'Yes' if p_rho < 0.05 else 'No',
                'significant_pearson': 'Yes' if p_r < 0.05 else 'No'
            })

    return pd.DataFrame(results), df


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='GenomeScope2 Statistical Analysis')
    parser.add_argument('--illumina_dir', required=True, help='Path to Illumina GenomeScope results')
    parser.add_argument('--aviti_dir', required=True, help='Path to Aviti GenomeScope results')
    parser.add_argument('--sex_file', required=True, help='Path to sex_summary.tsv (sample_id\\tsex)')
    parser.add_argument('--coverage_file', default=None, help='Path to qc_metrics_collected_for_plots.tsv')
    parser.add_argument('--output_dir', default='./statistical_results', help='Output directory')
    parser.add_argument('--iqr_multiplier', type=float, default=3.0, help='IQR multiplier for outlier removal')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. Collect data
    # -------------------------------------------------------------------------
    print("=" * 70)
    print("STEP 1: Collecting GenomeScope2 data")
    print("=" * 70)
    df = collect_genomescope_data(args.illumina_dir, args.aviti_dir)
    if len(df) == 0:
        print("ERROR: No data collected. Check paths.")
        sys.exit(1)

    n_illumina = len(df[df['platform'] == 'Illumina'])
    n_aviti = len(df[df['platform'] == 'Aviti'])
    print(f"\nTotal samples: {len(df)} (Illumina: {n_illumina}, Aviti: {n_aviti})")

    # -------------------------------------------------------------------------
    # 2. Outlier removal
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print(f"STEP 2: Outlier removal (IQR x {args.iqr_multiplier})")
    print("=" * 70)

    outlier_cols = ['heterozygosity', 'error_rate', 'model_fit', 'genome_length_gb']
    mask, outlier_info = remove_outliers_iqr(df, outlier_cols, multiplier=args.iqr_multiplier)

    n_removed = (~mask).sum()
    if outlier_info:
        print("Outliers detected:")
        for info in outlier_info:
            print(info)
    else:
        print("No outliers detected.")

    df_clean = df[mask].copy()
    print(f"Samples after outlier removal: {len(df_clean)} (removed {n_removed})")

    # Save outlier report
    outlier_report = df[~mask][['sample', 'platform', 'heterozygosity', 'error_rate', 'model_fit', 'genome_length_gb']]
    outlier_report.to_csv(os.path.join(args.output_dir, 'outliers_removed.csv'), index=False)

    # -------------------------------------------------------------------------
    # 3. Paired Wilcoxon signed-rank tests
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 3: Paired Wilcoxon signed-rank tests (Illumina vs Aviti)")
    print("=" * 70)

    paired_df = build_paired_df(df_clean)
    print(f"Number of paired individuals: {len(paired_df)}")

    if len(paired_df) >= 3:
        paired_results = paired_wilcoxon_tests(paired_df)
        print("\nPaired Wilcoxon Signed-Rank Test Results:")
        print("-" * 60)
        for _, row in paired_results.iterrows():
            print(f"  {row['metric']:25s} | n={row['n_pairs']:2.0f} | "
                  f"Illumina median={row['illumina_median']:.4f} | "
                  f"Aviti median={row['aviti_median']:.4f} | "
                  f"W={row['statistic']:.1f} | p={row['p_value']:.4f} | "
                  f"r={row['effect_size_r']:.3f} | sig={row['significant']}")
        paired_results.to_csv(os.path.join(args.output_dir, 'paired_wilcoxon_results.csv'), index=False)
        paired_df.to_csv(os.path.join(args.output_dir, 'paired_data.csv'), index=False)
    else:
        print("  Insufficient paired samples for statistical testing (need >= 3).")
        paired_results = pd.DataFrame()

    # -------------------------------------------------------------------------
    # 4. Sex comparison
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 4: Sex comparison (Males vs Females)")
    print("=" * 70)

    if os.path.exists(args.sex_file):
        sex_df = pd.read_csv(args.sex_file, sep='\t', header=None, names=['sample_id', 'sex'])
        sex_df['sex'] = sex_df['sex'].str.strip().str.upper()
        print(f"Sex data loaded: {len(sex_df)} samples")
        print(f"  Males: {(sex_df['sex'] == 'M').sum()}, Females: {(sex_df['sex'] == 'F').sum()}")

        sex_results, df_with_sex = sex_comparison(df_clean, sex_df)
        if len(sex_results) > 0:
            print("\nMann-Whitney U Test Results (Males vs Females):")
            print("-" * 80)
            for _, row in sex_results.iterrows():
                print(f"  {row['platform']:10s} | {row['metric']:25s} | "
                      f"M(n={row['n_males']:.0f})={row['male_median']:.4f} | "
                      f"F(n={row['n_females']:.0f})={row['female_median']:.4f} | "
                      f"U={row['statistic']:.1f} | p={row['p_value']:.4f} | "
                      f"sig={row['significant']}")
            sex_results.to_csv(os.path.join(args.output_dir, 'sex_comparison_results.csv'), index=False)
        else:
            print("  No sex comparison results generated.")
    else:
        print(f"  Sex file not found: {args.sex_file}")
        print("  Skipping sex comparison.")
        sex_results = pd.DataFrame()
        df_with_sex = df_clean.copy()
        df_with_sex['sex'] = np.nan

    # -------------------------------------------------------------------------
    # 5. Correlation analysis (coverage vs metrics)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 5: Correlation analysis (coverage vs GenomeScope metrics)")
    print("=" * 70)

    coverage_df = None
    if args.coverage_file and os.path.exists(args.coverage_file):
        coverage_df = pd.read_csv(args.coverage_file, sep='\t')
        # Standardize column names
        coverage_df.columns = [c.strip().lower() for c in coverage_df.columns]
        print(f"  Coverage file columns: {list(coverage_df.columns)}")

        # Try to find the coverage column under various names
        cov_col = None
        for candidate in ['coverage', 'cov', 'mean_coverage', 'avg_cov', 'depth',
                          'mean_cov', 'meancoverage', 'avg_coverage']:
            if candidate in coverage_df.columns:
                cov_col = candidate
                break

        if cov_col is None:
            # Try partial match
            for col in coverage_df.columns:
                if 'cov' in col or 'depth' in col:
                    cov_col = col
                    break

        if cov_col is not None:
            if cov_col != 'coverage':
                coverage_df = coverage_df.rename(columns={cov_col: 'coverage'})
            coverage_df['coverage'] = pd.to_numeric(coverage_df['coverage'], errors='coerce')
            print(f"  Coverage column identified: '{cov_col}', {len(coverage_df)} entries")
        else:
            print(f"  WARNING: Could not identify coverage column in: {list(coverage_df.columns)}")
            print(f"  Skipping correlation analysis.")
            coverage_df = None
    else:
        print(f"  Coverage file not found or not specified. Skipping correlation analysis.")

    corr_results, df_with_cov = correlation_analysis(df_clean, coverage_df)
    if len(corr_results) > 0 and not corr_results['spearman_rho'].isna().all():
        print("\nSpearman Correlation Results (Coverage vs Metrics):")
        print("-" * 80)
        for _, row in corr_results.iterrows():
            print(f"  {row['platform']:10s} | {row['metric']:25s} | n={row['n']:.0f} | "
                  f"rho={row['spearman_rho']:.3f} (p={row['spearman_p']:.4f}) | "
                  f"r={row['pearson_r']:.3f} (p={row['pearson_p']:.4f}) | "
                  f"sig(rho)={row['significant_spearman']}")
        corr_results.to_csv(os.path.join(args.output_dir, 'correlation_coverage_results.csv'), index=False)
    else:
        print("  No valid correlation results (coverage data may be missing).")

    # -------------------------------------------------------------------------
    # 6. Summary output
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 6: Summary statistics")
    print("=" * 70)

    summary = df_clean.groupby('platform').agg({
        'heterozygosity': ['count', 'mean', 'std', 'median'],
        'error_rate': ['mean', 'std', 'median'],
        'model_fit': ['mean', 'std', 'median'],
        'genome_length_gb': ['mean', 'std', 'median']
    }).round(4)
    print("\nDescriptive Statistics (after outlier removal):")
    print(summary)
    summary.to_csv(os.path.join(args.output_dir, 'descriptive_statistics.csv'))

    # Save full cleaned data
    df_clean.to_csv(os.path.join(args.output_dir, 'genomescope_clean_data.csv'), index=False)
    df_with_cov.to_csv(os.path.join(args.output_dir, 'genomescope_data_with_coverage.csv'), index=False)

    print(f"\n{'=' * 70}")
    print(f"DONE. All results saved to: {args.output_dir}")
    print(f"{'=' * 70}")
    print("Files generated:")
    for f in sorted(os.listdir(args.output_dir)):
        print(f"  {f}")


if __name__ == '__main__':
    main()
