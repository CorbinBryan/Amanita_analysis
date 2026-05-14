#!/bin/bash


for FILE in /Users/corbinbryan/Desktop/Amanita_analysis/buscos/busco_outputs/*.tar.gz; do
    GENOME=$(basename ${FILE} | sed 's/.tar.gz//g' | sed 's/run_//g')
    C_VAL=$(tar -xOf ${FILE} run_${GENOME}/short_summary_${GENOME}.txt | grep 'C:')
    echo -e "$GENOME\t$C_VAL"
done