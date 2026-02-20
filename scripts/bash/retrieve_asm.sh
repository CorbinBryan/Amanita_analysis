#!/bin/bash

# BASEDIR="$1" 

# for DIR in $(ls ${BASEDIR}); do 
#     if [[ -d ${DIR} ]]; then
#         ACC="$DIR"
#         mv ${BASEDIR}/${DIR}/assembly/contigs.fasta ${BASEDIR}/${ACC}_contigs.fasta
#         mv ${BASEDIR}/${DIR}/assembly/scaffolds.fasta ${BASEDIR}/${ACC}_scaffolds.fasta
#     fi 
# done

#!/bin/bash

BASEDIR="$1" 

for DIR in "${BASEDIR}"/*; do 
    if [[ -d "$DIR" ]]; then
        ACC="$(basename "$DIR")"
        mv "$DIR/assembly/contigs.fasta" "${BASEDIR}/${ACC}_contigs.fasta"
        mv "$DIR/assembly/scaffolds.fasta" "${BASEDIR}/${ACC}_scaffolds.fasta"
    fi 
done