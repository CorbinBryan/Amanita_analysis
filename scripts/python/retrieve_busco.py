#!/usr/bin/env python3
"""
Script to concatenate BUSCO orthogroups found in at least 3/4 of genomes.
Creates a single FASTA file with sequences named after their source sample.
"""

import os
import sys
import tarfile
import tempfile
import shutil
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set, Tuple
import argparse


def extract_sample_name(tar_filename: str) -> str:
    """Extract sample name from tar filename."""
    # Remove .tar.gz extension
    name = tar_filename.rsplit('.tar.gz', 1)[0]
    # Remove 'run_' prefix if present
    if name.startswith('run_'):
        name = name[4:]
    return name


def process_busco_archives(busco_dir: str) -> Tuple[Dict[str, Dict[str, list]], int]:
    """
    Process BUSCO tar.gz archives and collect orthogroup sequences.
    
    Returns:
        Tuple of (Dict mapping orthogroup_id -> {sample_id -> [sequences]}, genome_count)
    """
    busco_dir = Path(busco_dir).resolve()
    orthogroups = defaultdict(lambda: defaultdict(list))
    genome_count = 0
    
    # Get all tar.gz files
    tar_files = sorted(busco_dir.glob('*.tar.gz'))
    total_files = len(tar_files)
    
    print(f"Found {total_files} BUSCO archives\n")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        for idx, tar_path in enumerate(tar_files, 1):
            sample_name = extract_sample_name(tar_path.name)
            print(f"[{idx}/{total_files}] Processing {tar_path.name}...", end=' ')
            
            try:
                # Extract tar.gz to temporary directory
                with tarfile.open(tar_path, 'r:gz') as tar:
                    tar.extractall(tmpdir)
                
                # Find the extracted directory (usually named after the tar file)
                extracted_dirs = list(tmpdir.glob('*busco*'))
                if not extracted_dirs:
                    print("ERROR: No BUSCO directory found")
                    continue
                
                busco_dir_extracted = extracted_dirs[0]
                seq_dir = busco_dir_extracted / 'single_copy_busco_sequences'
                
                if not seq_dir.exists():
                    print("ERROR: No single_copy_busco_sequences directory")
                    shutil.rmtree(busco_dir_extracted)
                    continue
                
                # Process only DNA sequence files (.fna) in this genome
                ortho_files = list(seq_dir.glob('*.fna'))
                orthogroup_ids = set()
                
                for seq_file in ortho_files:
                    ortho_id = seq_file.stem
                    orthogroup_ids.add(ortho_id)
                    
                    # Read the FASTA file
                    with open(seq_file, 'r') as f:
                        content = f.read()
                    
                    orthogroups[ortho_id][sample_name].append(content)
                
                print(f"OK ({len(orthogroup_ids)} orthogroups)")
                genome_count += 1
                
                # Clean up extracted directory
                shutil.rmtree(busco_dir_extracted)
                
            except Exception as e:
                print(f"ERROR: {e}")
                continue
    
    print(f"\nProcessed {genome_count} genomes successfully")
    print(f"Found {len(orthogroups)} unique orthogroups\n")
    
    return orthogroups, genome_count


def update_fasta_headers(fasta_content: str, sample_name: str, ortho_id: str) -> str:
    """
    Update FASTA headers to include sample name.
    Format: >ORTHO_ID|SAMPLE_NAME|ORIGINAL_HEADER_INFO
    """
    lines = fasta_content.strip().split('\n')
    updated_lines = []
    
    for line in lines:
        if line.startswith('>'):
            # Original header format: >ortho_id:sample:contig:coords
            # New format: >ortho_id|sample_name|original_info
            header_info = line[1:]  # Remove '>'
            parts = header_info.split(':', 1)
            if len(parts) > 1:
                original_info = parts[1]
                new_header = f">{ortho_id}|{sample_name}|{original_info}"
            else:
                new_header = f">{ortho_id}|{sample_name}"
            updated_lines.append(new_header)
        else:
            updated_lines.append(line)
    
    return '\n'.join(updated_lines)


def write_individual_fastas(orthogroups: Dict[str, Dict[str, list]], 
                            output_dir: str) -> int:
    """
    Write each orthogroup to a separate FASTA file in the output directory.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    total_sequences = 0
    
    # Sort orthogroup IDs for consistent output
    for ortho_id in sorted(orthogroups.keys()):
        samples = orthogroups[ortho_id]
        output_file = output_dir / f"{ortho_id}.fasta"
        
        with open(output_file, 'w') as out:
            # Sort samples for consistency
            for sample_name in sorted(samples.keys()):
                sequences = samples[sample_name]
                
                for seq_content in sequences:
                    # Update headers with sample name
                    updated_content = update_fasta_headers(seq_content, sample_name, ortho_id)
                    out.write(updated_content + '\n')
                    total_sequences += 1
    
    return total_sequences


def main():
    parser = argparse.ArgumentParser(
        description='Extract BUSCO orthogroups found in at least 3/4 of genomes and create individual FASTA files.'
    )
    parser.add_argument(
        'busco_dir',
        help='Directory containing BUSCO tar.gz files'
    )
    parser.add_argument(
        '-o', '--output',
        default='busco_orthogroups',
        help='Output directory for individual FASTA files (default: busco_orthogroups)'
    )
    parser.add_argument(
        '-p', '--proportion',
        type=float,
        default=0.75,
        help='Minimum proportion of genomes to include (default: 0.75 = 3/4)'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    busco_dir = Path(args.busco_dir)
    if not busco_dir.exists():
        print(f"ERROR: Directory not found: {busco_dir}")
        sys.exit(1)
    
    if not busco_dir.is_dir():
        print(f"ERROR: Not a directory: {busco_dir}")
        sys.exit(1)
    
    # Count available tar files
    tar_files = list(busco_dir.glob('*.tar.gz'))
    if not tar_files:
        print(f"ERROR: No .tar.gz files found in {busco_dir}")
        sys.exit(1)
    
    print(f"BUSCO Orthogroup Extraction Script")
    print(f"=" * 50)
    print(f"Directory: {busco_dir}")
    print(f"Output directory: {args.output}")
    print(f"Proportion threshold: {args.proportion * 100:.0f}%\n")
    
    # Process archives
    orthogroups, processed_genomes = process_busco_archives(str(busco_dir))
    
    # Calculate minimum count based on actually processed genomes
    min_count = max(1, int(processed_genomes * args.proportion))
    
    print(f"Recalculating threshold based on {processed_genomes} successfully processed genomes...")
    print(f"Minimum genomes required: {min_count}\n")
    
    # Re-filter orthogroups based on actual genome count
    filtered_orthogroups = {}
    for ortho_id, samples in orthogroups.items():
        if len(samples) >= min_count:
            filtered_orthogroups[ortho_id] = samples
    
    if not filtered_orthogroups:
        print("ERROR: No orthogroups found meeting the minimum criteria")
        sys.exit(1)
    
    # Write output
    print(f"Writing individual FASTA files to {args.output}...")
    total_seqs = write_individual_fastas(filtered_orthogroups, args.output)
    
    output_path = Path(args.output).resolve()
    num_files = len(list(output_path.glob('*.fasta')))
    
    print(f"\nSuccess!")
    print(f"=" * 50)
    print(f"Output directory: {output_path}")
    print(f"Total FASTA files created: {num_files}")
    print(f"Total sequences written: {total_seqs}")
    print(f"Total orthogroups: {len(filtered_orthogroups)}")



if __name__ == '__main__':
    main()
