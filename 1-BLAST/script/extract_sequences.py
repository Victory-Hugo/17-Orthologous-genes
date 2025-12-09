#!/usr/bin/env python3
"""
Extract sequences from FAA files based on BLASTP results.

Usage:
    python3 extract_sequences.py <blastp_results> <faa_path_list> <output_dir>
"""

import sys
import os


def read_fasta(fasta_file):
    """
    Read FASTA file and return a dictionary of sequences.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Get first part of header
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def load_faa_files(faa_path_list):
    """
    Load all FAA files and create a mapping of filename to sequences.
    
    Args:
        faa_path_list: Path to file containing FAA file paths
        
    Returns:
        Dictionary mapping filename to sequences dictionary
    """
    faa_map = {}
    
    with open(faa_path_list, 'r') as f:
        for line in f:
            faa_file = line.strip()
            if faa_file and os.path.isfile(faa_file):
                filename = os.path.basename(faa_file)
                faa_map[filename] = read_fasta(faa_file)
    
    return faa_map


def extract_sequences(blastp_results, faa_path_list, output_dir):
    """
    Extract matching sequences from FAA files based on BLASTP results.
    
    Args:
        blastp_results: Path to filtered BLASTP results TSV file
        faa_path_list: Path to file containing FAA file paths
        output_dir: Output directory for extracted sequences
    """
    # Load FAA files
    print("Loading FAA files...")
    faa_map = load_faa_files(faa_path_list)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process BLASTP results
    print("Extracting sequences...")
    extracted_count = 0
    
    with open(blastp_results, 'r') as f:
        # Skip header
        next(f)
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            
            source_file = fields[0]
            sseqid = fields[2]  # Subject sequence ID
            
            # Extract source identifier (remove .faa extension)
            source_id = os.path.splitext(source_file)[0]
            
            # Get sequences from the source file
            if source_file in faa_map:
                sequences = faa_map[source_file]
                
                if sseqid in sequences:
                    # Create output filename with source information
                    # Format: {source_id}.faa
                    output_file = os.path.join(output_dir, f"{source_id}.faa")
                    
                    # Write sequence with source information in header
                    with open(output_file, 'w') as out_f:
                        # Header format: >{source_id}|{sseqid}
                        out_f.write(f">{source_id}|{sseqid}\n")
                        seq = sequences[sseqid]
                        # Write sequence in 60-character lines
                        for i in range(0, len(seq), 60):
                            out_f.write(seq[i:i+60] + '\n')
                    
                    extracted_count += 1
                else:
                    print(f"Warning: Sequence {sseqid} not found in {source_file}", file=sys.stderr)
            else:
                print(f"Warning: File {source_file} not found in FAA map", file=sys.stderr)
    
    print(f"Extracted {extracted_count} sequences to {output_dir}")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: extract_sequences.py <blastp_results> <faa_path_list> <output_dir>')
        sys.exit(1)
    
    blastp_results = sys.argv[1]
    faa_path_list = sys.argv[2]
    output_dir = sys.argv[3]
    
    extract_sequences(blastp_results, faa_path_list, output_dir)
