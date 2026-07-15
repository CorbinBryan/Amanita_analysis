from Bio import SeqIO

seqs = list(SeqIO.parse("concatenated_for_hyde.fasta", "fasta"))

n_ind = len(seqs)
n_sites = len(seqs[0].seq)

with open("data.txt", "w") as out:
    out.write(f"{n_ind} {n_sites}\n")
    
    for s in seqs:
        seq = str(s.seq).upper()
        seq = seq.replace("N", "-")  # HyDe prefers no Ns
        out.write(f"{s.id} {seq}\n")

print("✅ data.txt created")