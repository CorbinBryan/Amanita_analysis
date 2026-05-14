#!/usr/bin/env python3

from Bio import AlignIO
import glob
import itertools
import math

# -----------------------------
# Parameters (EDIT THESE)
# -----------------------------
ALIGNMENT_DIR = "./buscos/aln/al_*.fasta"
OUTPUT_MATRIX = "ani_matrix.tsv"
OUTPUT_COUNTS = "busco_counts.tsv"

MIN_BUSCOS = 150   # minimum number of shared BUSCOs required

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
# Main computation
# -----------------------------
results = {}

alignment_files = glob.glob(ALIGNMENT_DIR)
if len(alignment_files) == 0:
    raise Exception("No alignment files found.")

print(f"Found {len(alignment_files)} BUSCO alignments")

for aln_file in alignment_files:
    try:
        aln = AlignIO.read(aln_file, "fasta")
    except Exception as e:
        print(f"Skipping {aln_file}: {e}")
        continue

    # Skip BUSCOs with fewer than 2 species
    if len(aln) < 2:
        continue

    # Compute pairwise identities
    for a, b in itertools.combinations(aln, 2):
        sp_a = a.id.strip()
        sp_b = b.id.strip()

        pair = tuple(sorted([sp_a, sp_b]))

        pid = pairwise_identity(str(a.seq), str(b.seq))
        if pid is None:
            continue

        if pair not in results:
            results[pair] = {"sum": 0.0, "count": 0}

        results[pair]["sum"] += pid
        results[pair]["count"] += 1


# -----------------------------
# Get species list
# -----------------------------
species = sorted(set([s for pair in results for s in pair]))
print(f"Detected {len(species)} species")

# -----------------------------
# Initialize matrices
# -----------------------------
identity_matrix = {
    s: {t: 1.0 if s == t else math.nan for t in species}
    for s in species
}

count_matrix = {
    s: {t: 0 for t in species}
    for s in species
}

# -----------------------------
# Fill matrices
# -----------------------------
for pair, data in results.items():
    a, b = pair
    count = data["count"]

    if count >= MIN_BUSCOS:
        avg = data["sum"] / count
    else:
        avg = math.nan  # insufficient data

    identity_matrix[a][b] = avg
    identity_matrix[b][a] = avg

    count_matrix[a][b] = count
    count_matrix[b][a] = count

# -----------------------------
# Write identity matrix
# -----------------------------
with open(OUTPUT_MATRIX, "w") as out:
    out.write("\t" + "\t".join(species) + "\n")
    for s in species:
        row = [s]

        for t in species:
            val = identity_matrix[s][t]
            if math.isnan(val):
                row.append("NA")
            else:
                row.append(f"{val:.6f}")

        out.write("\t".join(row) + "\n")

# -----------------------------
# Write BUSCO count matrix
# -----------------------------
with open(OUTPUT_COUNTS, "w") as out:
    out.write("\t" + "\t".join(species) + "\n")
    for s in species:
        row = [s] + [str(count_matrix[s][t]) for t in species]
        out.write("\t".join(row) + "\n")

# -----------------------------
# Summary output
# -----------------------------
print(f"\n✅ Identity matrix written to: {OUTPUT_MATRIX}")
print(f"✅ BUSCO count matrix written to: {OUTPUT_COUNTS}")
print(f"✅ Minimum BUSCO threshold: {MIN_BUSCOS}")
