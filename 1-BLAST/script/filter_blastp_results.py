#!/usr/bin/env python3
"""
Filter BLASTP results by identity and coverage thresholds.

Usage:
    python3 filter_blastp_results.py <input_file> <output_file> <identity_threshold> <coverage_threshold>
"""

import sys


def filter_blastp_results(input_file, output_file, identity_threshold, coverage_threshold):
    """
    Filter BLASTP results and calculate query coverage.
    
    Args:
        input_file: Raw BLASTP output file with source filename as first column (tab-separated)
        output_file: Filtered results file
        identity_threshold: Minimum percent identity (0-100)
        coverage_threshold: Minimum query coverage percentage (0-100)
    """
    results = []

    # Read and filter results
    with open(input_file, 'r') as f:
        # Skip header line
        header = f.readline()
        
        for line in f:
            fields = line.strip().split('\t')
            # source_file is the first column, followed by BLASTP output columns
            if len(fields) < 15:  # 1 source_file + 14 BLASTP columns
                continue
            
            source_file = fields[0]
            blastp_fields = fields[1:]
            
            pident = float(blastp_fields[2])
            length = int(blastp_fields[3])
            qlen = int(blastp_fields[12])
            
            # Calculate query coverage
            query_coverage = (length / qlen) * 100
            
            # Apply filters
            if pident >= identity_threshold and query_coverage >= coverage_threshold:
                results.append((source_file, blastp_fields))

    # Sort by bitscore (highest first)
    results.sort(key=lambda x: float(x[1][11]), reverse=True)

    # Write filtered results
    with open(output_file, 'w') as f:
        f.write('\t'.join([
            'source_file', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
            'qlen', 'slen', 'query_coverage(%)'
        ]) + '\n')
        
        for source_file, blastp_fields in results:
            qlen = int(blastp_fields[12])
            length = int(blastp_fields[3])
            query_coverage = (length / qlen) * 100
            f.write(source_file + '\t' + '\t'.join(blastp_fields) + '\t' + f'{query_coverage:.2f}' + '\n')


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print('Usage: filter_blastp_results.py <input_file> <output_file> <identity_threshold> <coverage_threshold>')
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    identity_threshold = float(sys.argv[3])
    coverage_threshold = float(sys.argv[4])
    
    filter_blastp_results(input_file, output_file, identity_threshold, coverage_threshold)
