#!/bin/bash

DIR="$1"

for FILE in $(ls ${DIR}); do
    BASE=$(echo "$FILE" | cut -f1 -d ".")
    EXTENSION=$(echo "$FILE" | cut -f2 -d ".")
    if [[ $EXTENSION == "fasta" ]]; then
        continue
    elif [[ $EXTENSION == "fna" ]]; then
        mv ${DIR}/${FILE} ${DIR}/${BASE}.fasta
        tar -czf ${DIR}/${BASE}.fasta.tar.gz ${DIR}/${BASE}.fasta
    elif [[ $EXTENSION == "tar" ]]; then
        mv ${DIR}/${FILE} ${DIR}/${BASE}.fasta.tar.gz
    fi 
done