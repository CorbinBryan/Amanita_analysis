#!/bin/bash

for NAME in $(awk 'BEGIN{FS = ","}{print $1}' pop_asign_hyde.txt | sed 's/.fasta//g' | sed 's/.fna//g'); do 
    NEW_NAME=$(echo "$NAME" | cut -f1 -d "_")
    sed -i '' "s/$NAME/$NEW_NAME/g" concat_for_quibl.treefile
done