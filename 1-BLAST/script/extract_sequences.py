#!/usr/bin/env python3
"""
Extract sequences from FAA files based on BLASTP results (Streaming version).

Usage:
    python3 extract_sequences.py <blastp_results> <faa_path_list> <output_dir>
"""

import sys
import os
from collections import defaultdict


def build_faa_path_map(faa_path_list):
    """
    Build a mapping of filename to full file path.
    
    Args:
        faa_path_list: Path to file containing FAA file paths
        
    Returns:
        Dictionary mapping filename to full path
    """
    faa_path_map = {}
    
    with open(faa_path_list, 'r') as f:
        for line in f:
            faa_file = line.strip()
            if faa_file and os.path.isfile(faa_file):
                filename = os.path.basename(faa_file)
                faa_path_map[filename] = faa_file
    
    return faa_path_map


def collect_required_sequences(blastp_results):
    """
    Collect which sequences we need from which files.
    
    Args:
        blastp_results: Path to filtered BLASTP results TSV file
        
    Returns:
        Dictionary mapping source_file to set of required sequence IDs
    """
    required = defaultdict(set)
    
    with open(blastp_results, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            
            source_file = fields[0]
            sseqid = fields[2]  # Subject sequence ID
            
            required[source_file].add(sseqid)
    
    return required


def extract_sequence_from_file(faa_file, target_ids, output_dir):
    """
    Stream through a single FAA file and extract required sequences.
    
    Args:
        faa_file: Path to FAA file
        target_ids: Set of sequence IDs to extract
        output_dir: Output directory for extracted sequences
        
    Returns:
        Number of sequences extracted
    """
    source_file = os.path.basename(faa_file)
    source_id = os.path.splitext(source_file)[0]
    
    extracted_count = 0
    current_id = None
    current_seq = []
    writing = False
    
    with open(faa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Process previous sequence if we were extracting it
                if writing and current_id and current_seq:
                    output_file = os.path.join(output_dir, f"{source_id}.faa")
                    with open(output_file, 'w') as out_f:
                        out_f.write(f">{source_id}|{current_id}\n")
                        seq = ''.join(current_seq)
                        for i in range(0, len(seq), 60):
                            out_f.write(seq[i:i+60] + '\n')
                    extracted_count += 1
                
                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
                writing = current_id in target_ids
            else:
                if writing:
                    current_seq.append(line)
        
        # Process last sequence
        if writing and current_id and current_seq:
            output_file = os.path.join(output_dir, f"{source_id}.faa")
            with open(output_file, 'w') as out_f:
                out_f.write(f">{source_id}|{current_id}\n")
                seq = ''.join(current_seq)
                for i in range(0, len(seq), 60):
                    out_f.write(seq[i:i+60] + '\n')
            extracted_count += 1
    
    return extracted_count


def extract_sequences(blastp_results, faa_path_list, output_dir):
    """
    Extract matching sequences from FAA files based on BLASTP results (streaming).
    
    Args:
        blastp_results: Path to filtered BLASTP results TSV file
        faa_path_list: Path to file containing FAA file paths
        output_dir: Output directory for extracted sequences
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Build filename to path mapping
    print("Building file path mapping...")
    faa_path_map = build_faa_path_map(faa_path_list)
    
    # Collect required sequences per file
    print("Analyzing required sequences...")
    required_seqs = collect_required_sequences(blastp_results)
    
    # Stream process each file
    print("Extracting sequences (streaming)...")
    total_extracted = 0
    
    for source_file, target_ids in required_seqs.items():
        if source_file not in faa_path_map:
            print(f"Warning: File {source_file} not found in path list", file=sys.stderr)
            continue
        
        faa_file = faa_path_map[source_file]
        count = extract_sequence_from_file(faa_file, target_ids, output_dir)
        total_extracted += count
        print(f"  Processed {source_file}: {count} sequences extracted")
    
    print(f"Total: Extracted {total_extracted} sequences to {output_dir}")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: extract_sequences.py <blastp_results> <faa_path_list> <output_dir>')
        sys.exit(1)
    
    blastp_results = sys.argv[1]
    faa_path_list = sys.argv[2]
    output_dir = sys.argv[3]
    
    extract_sequences(blastp_results, faa_path_list, output_dir)
