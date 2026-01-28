#!/usr/bin/env python3
"""
REGENIE Results Summarizer
This script summarizes GWAS results from REGENIE output files

Usage:
    python summarize_results.py <results_file> [options]

Arguments:
    results_file    : REGENIE output file (.regenie)
    
Options:
    --p-threshold   : P-value threshold for significance (default: 5e-8)
    --output        : Output file prefix (default: summary)
    --top-n         : Number of top hits to report (default: 100)
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path


def calculate_lambda_gc(pvalues):
    """
    Calculate genomic inflation factor (lambda)
    
    Args:
        pvalues: Array of p-values
        
    Returns:
        Lambda GC value
    """
    # Convert p-values to chi-square statistics
    chisq = -2 * np.log(pvalues)
    
    # Calculate lambda as ratio of median observed to expected
    median_chisq = np.median(chisq)
    expected_median = 0.4549364  # Median of chi-square(1) distribution
    
    lambda_gc = median_chisq / expected_median
    
    return lambda_gc


def summarize_results(results_file, p_threshold=5e-8, output_prefix="summary", top_n=100):
    """
    Summarize GWAS results
    
    Args:
        results_file: Path to REGENIE results file
        p_threshold: Significance threshold
        output_prefix: Prefix for output files
        top_n: Number of top hits to report
    """
    
    print(f"Reading results from: {results_file}")
    
    # Read REGENIE results
    try:
        df = pd.read_csv(results_file, sep=r'\s+')
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    # Determine p-value column
    if 'LOG10P' in df.columns:
        df['P'] = 10 ** (-df['LOG10P'])
    elif 'P' in df.columns:
        pass
    else:
        print("Error: No p-value column found (expected 'LOG10P' or 'P')")
        sys.exit(1)
    
    # Calculate summary statistics
    n_variants = len(df)
    lambda_gc = calculate_lambda_gc(df['P'])
    min_p = df['P'].min()
    
    # Count significant variants
    n_significant = (df['P'] < p_threshold).sum()
    n_suggestive = (df['P'] < 1e-5).sum()
    
    # Get top hits
    top_hits = df.nsmallest(top_n, 'P')
    
    # Print summary
    print("\n" + "=" * 70)
    print("GWAS RESULTS SUMMARY")
    print("=" * 70)
    print(f"Total variants tested:              {n_variants:,}")
    print(f"Genomic inflation factor (λ):       {lambda_gc:.3f}")
    print(f"Minimum p-value:                    {min_p:.2e}")
    print(f"Genome-wide significant (p<5e-8):   {n_significant}")
    print(f"Suggestive associations (p<1e-5):   {n_suggestive}")
    print("=" * 70)
    
    # Interpretation
    print("\nINTERPRETATION:")
    if lambda_gc < 1.05:
        print("✓ Lambda is within acceptable range (good population structure control)")
    elif lambda_gc < 1.10:
        print("⚠ Lambda is slightly elevated (check for cryptic relatedness or stratification)")
    else:
        print("✗ Lambda is elevated (population stratification or other issues present)")
    
    # Save detailed summary
    summary_file = f"{output_prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("GWAS RESULTS SUMMARY\n")
        f.write("=" * 70 + "\n")
        f.write(f"Results file: {results_file}\n")
        f.write(f"Total variants: {n_variants:,}\n")
        f.write(f"Lambda GC: {lambda_gc:.3f}\n")
        f.write(f"Min p-value: {min_p:.2e}\n")
        f.write(f"Genome-wide significant: {n_significant}\n")
        f.write(f"Suggestive: {n_suggestive}\n")
    
    print(f"\nSummary saved to: {summary_file}")
    
    # Save top hits
    top_hits_file = f"{output_prefix}_top_{top_n}_hits.txt"
    top_hits.to_csv(top_hits_file, sep='\t', index=False)
    print(f"Top {top_n} hits saved to: {top_hits_file}")
    
    # Save significant hits
    if n_significant > 0:
        sig_hits = df[df['P'] < p_threshold]
        sig_file = f"{output_prefix}_significant.txt"
        sig_hits.to_csv(sig_file, sep='\t', index=False)
        print(f"Significant hits saved to: {sig_file}")
    
    # Display top 10 hits
    print("\nTOP 10 ASSOCIATIONS:")
    print("-" * 70)
    display_cols = ['CHROM', 'GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'BETA', 'SE', 'P']
    available_cols = [col for col in display_cols if col in top_hits.columns]
    
    if available_cols:
        print(top_hits[available_cols].head(10).to_string(index=False))
    else:
        print("Could not display top hits (column names may differ)")
    
    return {
        'n_variants': n_variants,
        'lambda_gc': lambda_gc,
        'min_p': min_p,
        'n_significant': n_significant,
        'n_suggestive': n_suggestive
    }


def main():
    """Main function"""
    
    parser = argparse.ArgumentParser(
        description='Summarize REGENIE GWAS results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python summarize_results.py results.regenie
    python summarize_results.py results.regenie --p-threshold 1e-7 --top-n 50
        """
    )
    
    parser.add_argument('results_file', help='REGENIE results file')
    parser.add_argument('--p-threshold', type=float, default=5e-8,
                       help='P-value threshold (default: 5e-8)')
    parser.add_argument('--output', default='summary',
                       help='Output prefix (default: summary)')
    parser.add_argument('--top-n', type=int, default=100,
                       help='Number of top hits (default: 100)')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.results_file).exists():
        print(f"Error: File not found: {args.results_file}")
        sys.exit(1)
    
    # Run summarization
    summarize_results(
        args.results_file,
        p_threshold=args.p_threshold,
        output_prefix=args.output,
        top_n=args.top_n
    )


if __name__ == '__main__':
    main()
