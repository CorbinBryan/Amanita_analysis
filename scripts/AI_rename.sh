#!/bin/bash

DIR="$1"

shopt -s nullglob
for FILEPATH in "$DIR"/*; do
    FILE=$(basename "$FILEPATH")
    BASE=$(echo "$FILE" | cut -f1 -d ".")
    EXTENSION=$(echo "$FILE" | cut -f2 -d ".")
    if [[ $EXTENSION == "fasta" ]]; then
        continue
    elif [[ $EXTENSION == "fna" ]]; then
        TARGET="${DIR}/${BASE}.fasta"
        ARCHIVE="${DIR}/${BASE}.fasta.tar.gz"

        if [[ -e "$TARGET" || -e "$ARCHIVE" ]]; then
            continue
        fi

        mv -- "$FILEPATH" "$TARGET"
        tar -czf "$ARCHIVE" -C "$DIR" "${BASE}.fasta"
    elif [[ $EXTENSION == "tar" ]]; then
        TARGET="${DIR}/${BASE}.fasta.tar.gz"

        if [[ -e "$TARGET" ]]; then
            continue
        fi

        mv -- "$FILEPATH" "$TARGET"
    fi
done