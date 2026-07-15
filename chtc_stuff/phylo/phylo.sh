ACCESS="$1"

mafft --auto ${ACCESS}.fasta --adjustdirection > al_${ACCESS}.fasta

trimal -automated1 -in al_${ACCESS}.fasta -out tral_${ACCESS}.fasta 

iqtree2 -s tral_${ACCESS}.fasta -pre ${ACCESS}