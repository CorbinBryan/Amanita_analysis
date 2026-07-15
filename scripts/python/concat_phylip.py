#!/usr/bin/env python3

from collections import OrderedDict
import os
from concurrent.futures import ProcessPoolExecutor


input_dir = "./buscos/alns/aln_for_hyde"
output_file = "./buscos/concatenated.phy"


def process_file(file_path):
    """Parse one PHYLIP file and return (n_sites, sequence_dict)."""
    with open(file_path) as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"Empty file: {file_path}")

    header_parts = lines[0].split()
    if len(header_parts) < 2:
        raise ValueError(f"Malformed header in {file_path}")

    n_sites = int(header_parts[1])

    seqs = {}
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        taxon, seq = parts[0], parts[1]
        seqs[taxon] = seq

    return n_sites, seqs


def main():
    phylip_files = sorted(
        [f for f in os.listdir(input_dir) if f.endswith(".phy")]
    )

    if not phylip_files:
        raise ValueError(f"No .phy files found in {input_dir}")

    file_paths = [os.path.join(input_dir, f) for f in phylip_files]

    # ✅ Parallel parsing
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_file, file_paths))

    # ✅ Collect ALL taxa
    all_taxa = OrderedDict()
    for _, seqs in results:
        for taxon in seqs:
            if taxon not in all_taxa:
                all_taxa[taxon] = []

    total_length = 0

    # ✅ Build concatenation
    for n_sites, seqs in results:
        for taxon in all_taxa:
            if taxon in seqs:
                all_taxa[taxon].append(seqs[taxon])
            else:
                all_taxa[taxon].append("N" * n_sites)
        total_length += n_sites

    # ✅ Validate lengths
    for taxon, parts in all_taxa.items():
        seq = "".join(parts)
        if len(seq) != total_length:
            raise ValueError(
                f"Length mismatch for {taxon}: {len(seq)} != {total_length}"
            )

    # ✅ Write PHYLIP
    with open(output_file, "w") as out:
        out.write(f"{len(all_taxa)} {total_length}\n")
        for taxon, parts in all_taxa.items():
            seq = "".join(parts)
            out.write(f"{taxon} {seq}\n")

    print(f"Done. Wrote {output_file}")


if __name__ == "__main__":
    main()