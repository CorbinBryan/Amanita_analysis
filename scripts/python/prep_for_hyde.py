import os
from Bio import SeqIO

# ==== USER SETTINGS ====
input_dir = "./buscos/alns/aln_filtered"      # folder with input alignments
output_dir = "./buscos/alns/aln_for_hyde"
keep_file = "./buscos/keep_for_hyde.txt"
file_extension = ".fasta"     # change if needed

# =======================

# Read taxa to keep
with open(keep_file) as f:
    taxa_to_keep = set(line.strip() for line in f if line.strip())

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Process each alignment file
for filename in os.listdir(input_dir):
    if not filename.endswith(file_extension):
        continue

    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, filename)

    kept_records = []

    for record in SeqIO.parse(input_path, "fasta"):
        # Match by ID (usually first word in FASTA header)
        if record.id in taxa_to_keep:
            kept_records.append(record)

    # Write filtered sequences
    SeqIO.write(kept_records, output_path, "fasta")

    print(f"Processed {filename}: kept {len(kept_records)} sequences")

print("✅ Done!")