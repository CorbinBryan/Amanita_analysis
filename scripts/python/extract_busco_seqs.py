#!/usr/bin/env python3

import os
import tarfile
import tempfile
from collections import defaultdict
from Bio import SeqIO

INPUT_DIR = "PATH_TO_DIRECTORY_WITH_TARS"
OUTPUT_DIR = "busco_single_copy_nt_output"
MIN_FRACTION = 0.75  # 3/4 threshold

os.makedirs(OUTPUT_DIR, exist_ok=True)

# busco_id -> list of SeqRecords
busco_dict = defaultdict(list)
genomes = []

# ---------------------------------------------------------
# STEP 1: Parse tar.gz BUSCO runs
# ---------------------------------------------------------

for fname in os.listdir(INPUT_DIR):
    if fname.startswith("run_") and fname.endswith("_busco.tar.gz"):
        genome = fname.replace("run_", "").replace("_busco.tar.gz", "")
        genomes.append(genome)

        tar_path = os.path.join(INPUT_DIR, fname)

        with tempfile.TemporaryDirectory() as tmpdir:
            with tarfile.open(tar_path, "r:gz") as tar:
                tar.extractall(tmpdir)

            # Locate single copy nucleotide sequences
            for root, dirs, files in os.walk(tmpdir):
                if "single_copy_busco_sequences" in root:
                    for file in files:
                        if file.endswith(".fna"):  # ONLY nucleotide
                            busco_id = file.replace(".fna", "")
                            file_path = os.path.join(root, file)

                            for record in SeqIO.parse(file_path, "fasta"):
                                record.id = genome
                                record.name = genome
                                record.description = ""
                                busco_dict[busco_id].append(record)

# ---------------------------------------------------------
# STEP 2: Apply 3/4 threshold
# ---------------------------------------------------------

total_genomes = len(genomes)
min_genomes_required = int(total_genomes * MIN_FRACTION + 0.999)

print(f"Total genomes: {total_genomes}")
print(f"Minimum genomes required: {min_genomes_required}")

for busco_id, records in busco_dict.items():
    if len(records) >= min_genomes_required:
        output_path = os.path.join(OUTPUT_DIR, f"{busco_id}.fasta")
        SeqIO.write(records, output_path, "fasta")

print("Done.")