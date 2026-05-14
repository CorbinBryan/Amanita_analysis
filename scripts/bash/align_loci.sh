#!/bin/bash

cd ../../buscos/busco_single_copy_nt_output

if [[ -d ../aln ]]; then 
    rm -r ../aln
    mkdir ../aln
else 
    mkdir ../aln
fi

ls *.fasta > BUSCO_list.txt

while read FASTA; do
    NUM_ACC=$(grep -c ">" "${FASTA}")
    echo "Processing $FASTA: $NUM_ACC sequences"
    if [[ $NUM_ACC -gt 54 ]]; then
        /Users/corbinbryan/anaconda3/envs/bio_tools/bin/mafft --adjustdirection --thread 8 "${FASTA}" > ../aln/al_${FASTA}
    fi
done < BUSCO_list.txt

rm BUSCO_list.txt