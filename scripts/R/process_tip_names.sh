#!/bin/bash

if [ -f tips_species.tsv ]; then
    rm tips_species.tsv
fi

for NAME in $(cat tip_names.txt); do
    if [[ "$NAME" =~ ^SRR ]]; then
        SEARCH=$(echo "$NAME" | cut -f1 -d "_")
        KEY=$(grep "$SEARCH" ../../srr_doc.csv | awk 'BEGIN{FS=","}{print $1}')
        echo -e "$NAME\t$KEY" >> tips_species.tsv
    elif [[ "$NAME" =~ ^GCA ]] && ! [[ $NAME =~ GCA_001983385 ]]; then
        SEARCH=$(echo "$NAME" | cut -f1 -d ".") 
        echo "searching for $SEARCH"
        KEY=$(grep "$SEARCH" ../../ncbi_dataset.tsv | awk 'BEGIN{FS="\t"}{print $4}')
        echo -e "$NAME\t$KEY" >> tips_species.tsv
    elif [[ "$NAME" =~ GCA_001983385 ]]; then
        echo -e "$NAME\tAmanita phalloides" >> tips_species.tsv
    elif [[ "$NAME" =~ ^Volvo1_AssemblyScaffolds ]]; then
        echo -e "$NAME\tVolvariella volvacea" >> tips_species.tsv
    elif [[ "$NAME" =~ ^Amapyr1 ]]; then
        echo -e "$NAME\tAmanita aff. conicoverrucosa" >> tips_species.tsv
    elif [[ "$NAME" =~ ^Amagr1 ]]; then
        echo -e "$NAME\tAmanita aff. grandis" >> tips_species.tsv
    elif [[ "$NAME" =~ ^30280 ]]; then
        echo -e "$NAME\tAmanita populiphila" >> tips_species.tsv
    elif [[ "$NAME" =~ ^40019 ]]; then
        echo -e "$NAME\tAmanita citrina" >> tips_species.tsv
    elif [[ "$NAME" =~ ^40225 ]]; then
        echo -e "$NAME\tAmanita baningiana" >> tips_species.tsv
    elif [[ "$NAME" =~ ^NesPan3 ]]; then
        echo -e "$NAME\tAmanita pantherina" >> tips_species.tsv
    fi
done 