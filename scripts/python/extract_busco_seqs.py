#!/usr/bin/env python3


import os
import tarfile
from collections import defaultdict
from Bio import SeqIO
import math
from multiprocessing import Pool, cpu_count
from io import TextIOWrapper

INPUT_DIR = "./busco_outputs"
OUTPUT_DIR = "busco_combined_fastas"

os.makedirs(OUTPUT_DIR, exist_ok=True)


def process_archive(fname):
    busco_sequences = defaultdict(list)

    genome_name = fname.replace("run_", "").replace("_busco.tar.gz", "")
    tar_path = os.path.join(INPUT_DIR, fname)

    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            # ✅ Only grab nucleotide BUSCO files
            if "single_copy_busco_sequences/" in member.name and member.name.endswith(".fna"):
                
                busco_id = os.path.splitext(os.path.basename(member.name))[0]

                f = tar.extractfile(member)
                if f is None:
                    continue

                fasta_handle = TextIOWrapper(f)

                for record in SeqIO.parse(fasta_handle, "fasta"):
                    record.id = genome_name
                    record.description = ""
                    busco_sequences[busco_id].append(record)

    busco_count = len(busco_sequences)
    return genome_name, busco_sequences, busco_count


def main():
    archives = [f for f in os.listdir(INPUT_DIR) if f.endswith(".tar.gz")]

    print(f"Found {len(archives)} archives")

    n_threads = min(cpu_count(), len(archives))
    print(f"Using {n_threads} processes")

    combined_buscos = defaultdict(list)
    genome_counts = {}
    genomes = []

    with Pool(n_threads) as pool:
        results = pool.map(process_archive, archives)

    for genome_name, busco_dict, count in results:
        genomes.append(genome_name)
        genome_counts[genome_name] = count

        for busco_id, records in busco_dict.items():
            combined_buscos[busco_id].extend(records)

    # Print counts
    print("\nBUSCO counts per genome:")
    for g, c in sorted(genome_counts.items()):
        print(f"{g}\t{c}")

    # Threshold
    threshold = math.ceil(0.75 * len(genomes))
    print(f"\nThreshold: {threshold} / {len(genomes)} genomes")

    # Write outputs
    kept = 0
    for busco_id, records in combined_buscos.items():
        if len(records) >= threshold:
            outpath = os.path.join(OUTPUT_DIR, f"{busco_id}.fasta")
            with open(outpath, "w") as out_f:
                SeqIO.write(records, out_f, "fasta")
            kept += 1

    print(f"\nTotal BUSCOs written: {kept}")


if __name__ == "__main__":
    main()
