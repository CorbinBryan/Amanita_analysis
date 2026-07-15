#!/usr/bin/env python3

from Bio import AlignIO
import glob
import itertools
import math
from multiprocessing import Pool, cpu_count

# -----------------------------
# Parameters
# -----------------------------
ALIGNMENT_DIR = "./buscos/aln/tral_*.fasta"
MIN_BUSCOS = 200

OUTPUT_MATRIX = "ani_matrix.tsv"
OUTPUT_COUNTS = "busco_counts.tsv"

# Number of cores (adjust if needed)
N_CORES = max(1, cpu_count() - 1)


# -----------------------------
# Function: pairwise identity
# -----------------------------
def pairwise_identity(seq1, seq2):
    matches = 0
    length = 0

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            length += 1
            if a == b:
                matches += 1

    return matches / length if length > 0 else None


# -----------------------------
# Process ONE BUSCO file
# -----------------------------
def process_alignment(aln_file):
    local_results = {}

    try:
        aln = AlignIO.read(aln_file, "fasta")
    except Exception:
        return local_results

    if len(aln) < 2:
        return local_results

    for a, b in itertools.combinations(aln, 2):
        sp_a = a.id.strip()
        sp_b = b.id.strip()

        pair = tuple(sorted([sp_a, sp_b]))
        pid = pairwise_identity(str(a.seq), str(b.seq))

        if pid is None:
            continue

        if pair not in local_results:
            local_results[pair] = {"sum": 0.0, "count": 0}

        local_results[pair]["sum"] += pid
        local_results[pair]["count"] += 1

    return local_results


# -----------------------------
# Merge results
# -----------------------------
def merge_results(global_res, local_res):
    for pair, data in local_res.items():
        if pair not in global_res:
            global_res[pair] = {"sum": 0.0, "count": 0}

        global_res[pair]["sum"] += data["sum"]
        global_res[pair]["count"] += data["count"]


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":

    alignment_files = glob.glob(ALIGNMENT_DIR)

    if len(alignment_files) == 0:
        raise Exception("No alignment files found.")

    print(f"Found {len(alignment_files)} BUSCO alignments")
    print(f"Using {N_CORES} CPU cores")

    # Run in parallel
    with Pool(N_CORES) as pool:
        results_list = pool.map(process_alignment, alignment_files)

    # Combine results
    results = {}
    for res in results_list:
        merge_results(results, res)

    print("Finished processing BUSCOs")

    # Get species
    species = sorted(set([s for pair in results for s in pair]))
    print(f"Detected {len(species)} species")

    # Initialize matrices
    identity_matrix = {
        s: {t: 1.0 if s == t else math.nan for t in species}
        for s in species
    }

    count_matrix = {
        s: {t: 0 for t in species}
        for s in species
    }

    # Fill matrices
    for pair, data in results.items():
        a, b = pair
        count = data["count"]

        if count >= MIN_BUSCOS:
            avg = data["sum"] / count
        else:
            avg = math.nan

        identity_matrix[a][b] = avg
        identity_matrix[b][a] = avg

        count_matrix[a][b] = count
        count_matrix[b][a] = count

    # Write identity matrix
    with open(OUTPUT_MATRIX, "w") as out:
        out.write("\t" + "\t".join(species) + "\n")
        for s in species:
            row = [s]
            for t in species:
                val = identity_matrix[s][t]
                row.append("NA" if math.isnan(val) else f"{val:.6f}")
            out.write("\t".join(row) + "\n")

    # Write count matrix
    with open(OUTPUT_COUNTS, "w") as out:
        out.write("\t" + "\t".join(species) + "\n")
        for s in species:
            row = [s] + [str(count_matrix[s][t]) for t in species]
            out.write("\t".join(row) + "\n")

    print("\n✅ Done!")