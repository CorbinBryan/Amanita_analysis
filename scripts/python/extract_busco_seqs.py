#!/usr/bin/env python3

import os
import tarfile
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------
INPUT_DIR = "/Users/corbinbryan/Desktop/Amanita_analysis/buscos/busco_outputs"
OUTPUT_DIR = "/Users/corbinbryan/Desktop/Amanita_analysis/buscos/busco_single_copy_nt_output"

MIN_FRACTION = 0.15
N_WORKERS = os.cpu_count() - 1  # adjust if needed

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------
# HELPERS
# ---------------------------------------------------------
def clean_genome_name(name):
    for ext in [".fna", ".fasta", ".fa"]:
        if name.endswith(ext):
            name = name[:-len(ext)]
    return name


def process_tarball(fname):
    """
    Extract BUSCO sequences from one tar.gz file.
    Returns:
        genome_name, dict(busco_id -> SeqRecord)
    """

    genome = fname.replace("run_", "").replace("_busco.tar.gz", "")
    genome = clean_genome_name(genome)

    tar_path = os.path.join(INPUT_DIR, fname)
    local_dict = {}

    with tarfile.open(tar_path, "r:gz") as tar:
        with tempfile.TemporaryDirectory() as tmpdir:
            tar.extractall(tmpdir)

            for root, _, files in os.walk(tmpdir):
                # Only process correct BUSCO directory
                if "single_copy_busco_sequences" not in root:
                    continue

                for file in files:
                    if not file.endswith(".fna"):
                        continue

                    fpath = os.path.join(root, file)

                    # Extract correct BUSCO ID
                    base = os.path.splitext(file)[0]
                    busco_id = base.split("_")[0]

                    for record in SeqIO.parse(fpath, "fasta"):
                        record.id = f"{genome}|{busco_id}"
                        record.description = ""

                        # Ensure one sequence per genome per BUSCO
                        if busco_id not in local_dict:
                            local_dict[busco_id] = record

    return genome, local_dict


# ---------------------------------------------------------
# MAIN PIPELINE
# ---------------------------------------------------------
def main():

    busco_dict = defaultdict(dict)
    genomes = []

    tar_files = [
        f for f in os.listdir(INPUT_DIR)
        if f.startswith("run_") and f.endswith("_busco.tar.gz")
    ]

    print(f"Found {len(tar_files)} BUSCO archives")
    print(f"Using {N_WORKERS} workers\n")

    # ---------------------------------------------------------
    # STEP 1: PARALLEL PROCESSING
    # ---------------------------------------------------------
    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {executor.submit(process_tarball, f): f for f in tar_files}

        for future in as_completed(futures):
            fname = futures[future]

            try:
                genome, result_dict = future.result()
            except Exception as e:
                print(f"ERROR processing {fname}: {e}")
                continue

            genomes.append(genome)

            for busco_id, record in result_dict.items():
                busco_dict[busco_id][genome] = record

            print(f"Finished: {genome}")

    # ---------------------------------------------------------
    # STEP 2: FILTERING
    # ---------------------------------------------------------
    total_genomes = len(genomes)
    min_required = int(total_genomes * MIN_FRACTION + 0.999)

    print(f"\nTotal genomes: {total_genomes}")
    print(f"Minimum required: {min_required}")
    print(f"Total BUSCOs found: {len(busco_dict)}")

    # ---------------------------------------------------------
    # STEP 3: WRITE OUTPUT
    # ---------------------------------------------------------
    kept = 0

    for busco_id, genome_map in busco_dict.items():
        if len(genome_map) < min_required:
            continue

        output_path = os.path.join(OUTPUT_DIR, f"{busco_id}.fasta")
        SeqIO.write(genome_map.values(), output_path, "fasta")
        kept += 1

    print(f"\nBUSCOs retained: {kept}")
    print("Done.")


# ---------------------------------------------------------
# ENTRY POINT (CRITICAL for macOS/Windows)
# ---------------------------------------------------------
if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()  # safe on all platforms
    main()