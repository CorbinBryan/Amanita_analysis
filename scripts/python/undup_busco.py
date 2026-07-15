#!/usr/bin/env python3

import tarfile
import os
import io
from collections import defaultdict

# ---- USER INPUT ----
input_tar = "run_Amaapr_busco_copy.tar.gz"
output_tar = "./busco_outputs/run_Amaapr1_AssemblyScaffolds_Repeatmasked.fasta_busco.tar.gz"
genome_fasta = "../Amaapr1_AssemblyScaffolds_Repeatmasked.fasta"

# ---------------------

# ---- Load genome ----
print("Loading genome...")
genome = {}

with open(genome_fasta) as f:
    seq_id = None
    seq_chunks = []

    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id:
                genome[seq_id] = "".join(seq_chunks)
            seq_id = line[1:].split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line)

    if seq_id:
        genome[seq_id] = "".join(seq_chunks)

print(f"Loaded genome: {len(genome)} sequences")

# ---- Open BUSCO tar ----
tar_in = tarfile.open(input_tar, "r:gz")
members = tar_in.getmembers()

# ---- Find full_table.tsv ----
full_table_member = next(m for m in members if m.name.endswith("full_table.tsv"))

# ---- Parse duplicated BUSCOs + scores ----
duplicated_buscos = set()
busco_scores = defaultdict(list)

with tar_in.extractfile(full_table_member) as f:
    for line in f:
        line = line.decode("utf-8")
        if line.startswith("#"):
            continue

        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue

        busco_id = parts[0]
        status = parts[1]

        if "Duplic" in status:
            duplicated_buscos.add(busco_id)

            # Try collecting score (column varies slightly by BUSCO version)
            try:
                score = float(parts[5])
                busco_scores[busco_id].append(score)
            except:
                pass

print(f"Duplicated BUSCOs: {len(duplicated_buscos)}")

# ---- Reverse complement ----
def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

# ---- Parse GFF attributes ----
def parse_attrs(attr_string):
    d = {}
    for item in attr_string.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

# ---- Extract sequences from GFF (with transcript grouping) ----
def extract_transcripts(gff_text):
    transcripts = defaultdict(list)

    for line in gff_text.splitlines():
        if line.startswith("#"):
            continue

        cols = line.split("\t")
        if len(cols) < 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, attrs = cols
        start, end = int(start), int(end)

        attr_dict = parse_attrs(attrs)

        parent = attr_dict.get("Parent") or attr_dict.get("ID") or "unknown"

        transcripts[parent].append(
            (feature, seqid, start, end, strand)
        )

    results = []

    for parent, feats in transcripts.items():
        cds_feats = [f for f in feats if f[0] == "CDS"]
        use_feats = cds_feats if cds_feats else feats

        seq_chunks = []

        for feature, seqid, start, end, strand in use_feats:
            if seqid not in genome:
                continue

            seq = genome[seqid][start - 1:end]

            if strand == "-":
                seq = revcomp(seq)

            seq_chunks.append((start, seq))

        if not seq_chunks:
            continue

        seq_chunks.sort(key=lambda x: x[0])
        sequence = "".join([s for _, s in seq_chunks])

        results.append(sequence)

    return results

# ---- Collect sequences ----
busco_candidates = defaultdict(list)

for m in members:
    if "multi_copy_busco_sequences" in m.name and m.name.endswith(".gff"):

        busco_id = os.path.basename(m.name).split(".")[0]

        if busco_id not in duplicated_buscos:
            continue

        with tar_in.extractfile(m) as f:
            gff_text = f.read().decode("utf-8")

        sequences = extract_transcripts(gff_text)

        for seq in sequences:
            busco_candidates[busco_id].append(seq)

print(f"BUSCOs with candidates: {len(busco_candidates)}")

# ---- Select best sequence ----
selected_sequences = {}

for busco_id, seqs in busco_candidates.items():

    # Prefer highest BUSCO score if available
    if busco_scores[busco_id]:
        best_seq = max(seqs, key=lambda s: len(s))  # fallback still length-based
    else:
        best_seq = max(seqs, key=lambda s: len(s))

    selected_sequences[busco_id] = best_seq

print(f"Selected BUSCOs: {len(selected_sequences)}")

# ---- Build new tar ----
tar_out = tarfile.open(output_tar, "w:gz")

# detect root
root_prefix = next((m.name.split("/")[0] for m in members if "/" in m.name), "")

# copy original contents
for m in members:
    fileobj = tar_in.extractfile(m) if m.isfile() else None
    tar_out.addfile(m, fileobj)

# add new sequences
for busco_id, seq in selected_sequences.items():

    new_path = os.path.join(
        root_prefix,
        "busco_sequences",
        "single_copy_busco_sequences",
        f"{busco_id}.fna"
    )

    fasta = f">{busco_id}\n{seq}\n"
    data = fasta.encode()

    ti = tarfile.TarInfo(new_path)
    ti.size = len(data)

    tar_out.addfile(ti, io.BytesIO(data))

tar_in.close()
tar_out.close()

print(f"✅ Done: {output_tar}")
