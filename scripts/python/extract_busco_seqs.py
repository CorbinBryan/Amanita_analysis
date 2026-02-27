#!/usr/bin/env python3
"""
extract_busco_by_gene.py

Scan BUSCO run outputs (compressed archives or directories) in a folder,
collect sequences for each BUSCO gene present in >= fraction of taxa,
and write one FASTA file per BUSCO gene with sequence headers named by genome.

Usage: python3 extract_busco_by_gene.py --input-dir . --out-dir busco_gene_fastas
"""
import argparse
import os
import tarfile
import zipfile
import tempfile
import shutil
from collections import defaultdict
import sys


def iter_fasta_records(path):
    """Yield (header, seq) for a FASTA file at path."""
    header = None
    seq_chunks = []
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line[0] == '>':
                if header is not None:
                    yield header, ''.join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, ''.join(seq_chunks)


def is_nucleotide_sequence(seq, threshold=0.85):
    """Return True if seq looks like a nucleotide sequence by fraction of nucleotide chars."""
    if not seq:
        return False
    seq = seq.upper()
    nuc_chars = set('ACGTURYSWKMBDHVN')
    letters = [c for c in seq if c.isalpha()]
    if not letters:
        return False
    nuc_count = sum(1 for c in letters if c in nuc_chars)
    return (nuc_count / len(letters)) >= threshold


def find_fasta_files(root):
    """Recursively find likely BUSCO fasta files under root."""
    exts = {'.fa', '.faa', '.fasta', '.fas', '.fna'}
    candidates = []
    for dirpath, dirs, files in os.walk(root):
        for f in files:
            lf = f.lower()
            name, e = os.path.splitext(lf)
            if e in exts:
                candidates.append(os.path.join(dirpath, f))
            # also include files with typical BUSCO nucleotide labels
            if 'nuc' in lf or 'nucl' in lf or 'cds' in lf or 'transcripts' in lf:
                candidates.append(os.path.join(dirpath, f))
    return candidates


def safe_taxon_name(name):
    # make a filesystem- and fasta-safe taxon name
    return name.replace(' ', '_').replace('/', '_')


def extract_archive(archive_path, dest_dir):
    if tarfile.is_tarfile(archive_path):
        try:
            with tarfile.open(archive_path, 'r:*') as tf:
                tf.extractall(dest_dir)
            return True
        except Exception:
            return False
    elif zipfile.is_zipfile(archive_path):
        try:
            with zipfile.ZipFile(archive_path, 'r') as zf:
                zf.extractall(dest_dir)
            return True
        except Exception:
            return False
    else:
        return False


def get_gene_id_from_header(header):
    # Use the first token as the gene id; strip common prefixes
    token = header.split()[0]
    # remove leading 'lcl|' or similar
    if '|' in token:
        token = token.split('|')[-1]
    return token


def gather_taxon_sequences(path, taxon_name):
    """Return dict gene_id -> list of sequences (header, seq) found in this taxon."""
    seqs = defaultdict(list)
    if os.path.isdir(path):
        root = path
    else:
        # path might be an archive file
        tempd = tempfile.mkdtemp(prefix='busco_extract_')
        ok = extract_archive(path, tempd)
        if not ok:
            shutil.rmtree(tempd)
            return seqs
        root = tempd
    try:
        fasta_files = find_fasta_files(root)
        # also try to find any file names containing 'busco' to narrow scope
        if not fasta_files:
            for dirpath, dirs, files in os.walk(root):
                for f in files:
                    if 'busco' in f.lower() or 'single_copy' in f.lower():
                        fp = os.path.join(dirpath, f)
                        _, e = os.path.splitext(fp.lower())
                        if e in {'.fa', '.faa', '.fasta', '.fas', '.fna'}:
                            fasta_files.append(fp)

        for fp in fasta_files:
            for header, seq in iter_fasta_records(fp):
                # only keep nucleotide-like sequences
                if not is_nucleotide_sequence(seq):
                    continue
                gid = get_gene_id_from_header(header)
                seqs[gid].append((header, seq))
    finally:
        if not os.path.isdir(path):
            shutil.rmtree(root)
    return seqs


def main():
    p = argparse.ArgumentParser(description='Collect BUSCO gene FASTAs present in a fraction of taxa')
    p.add_argument('--input-dir', '-i', default='.', help='Directory with BUSCO output archives or folders')
    p.add_argument('--out-dir', '-o', default='busco_gene_fastas', help='Output directory for per-gene FASTAs')
    p.add_argument('--min-fraction', '-f', type=float, default=0.75, help='Minimum fraction of taxa that must contain the gene (default 0.75)')
    p.add_argument('--extensions', '-e', nargs='*', default=['.tar.gz', '.tgz', '.zip', '.gz', '.tar'], help='Archive extensions to consider')
    args = p.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # find candidate taxa entries: either directories in input_dir or archive files
    entries = []
    for name in os.listdir(input_dir):
        if name.startswith('.'):
            continue
        full = os.path.join(input_dir, name)
        if os.path.isdir(full):
            entries.append((name, full))
        else:
            lname = name.lower()
            for ext in args.extensions:
                if lname.endswith(ext):
                    # taxon name inferred from archive basename without ext
                    tname = name
                    # strip multiple known extensions
                    for e in ['.tar.gz', '.tgz', '.zip', '.gz', '.tar']:
                        if tname.lower().endswith(e):
                            tname = tname[: -len(e)]
                    entries.append((tname, full))
                    break

    if not entries:
        print('No archives or directories found in', input_dir, file=sys.stderr)
        sys.exit(1)

    print('Processing', len(entries), 'taxa...')

    # species -> gene -> [(header, seq), ...]
    gene_map = defaultdict(lambda: defaultdict(list))
    taxa_names = []
    for taxon_name, path in entries:
        taxa_names.append(taxon_name)
        print('  scanning', taxon_name)
        taxon_seqs = gather_taxon_sequences(path, taxon_name)
        for gid, recs in taxon_seqs.items():
            gene_map[gid][taxon_name].extend(recs)

    n_taxa = len(taxa_names)
    min_count = int((args.min_fraction * n_taxa) + 0.999999)
    print('Found', len(gene_map), 'unique gene IDs across taxa')
    print('Writing genes present in >=', min_count, '/', n_taxa, 'taxa')

    for gid, taxon_dict in gene_map.items():
        present = len(taxon_dict)
        if present < min_count:
            continue
        out_fp = os.path.join(out_dir, f'{gid}.fasta')
        with open(out_fp, 'w', encoding='utf-8') as outfh:
            for taxon in sorted(taxon_dict.keys()):
                recs = taxon_dict[taxon]
                if len(recs) == 1:
                    header, seq = recs[0]
                    outfh.write(f'>{safe_taxon_name(taxon)}\n')
                    outfh.write(seq + '\n')
                else:
                    # multiple copies: enumerate
                    for i, (header, seq) in enumerate(recs, start=1):
                        outfh.write(f'>{safe_taxon_name(taxon)}|copy{i}\n')
                        outfh.write(seq + '\n')

    print('Done. Per-gene FASTAs are in', out_dir)


if __name__ == '__main__':
    main()
