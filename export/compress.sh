#!/bin/bash

for i in $(ls); do
    if [[ "$i" =~ tar.gz$ ]]; then
        continue
    elif [[ "$i" =~ gz$ ]] && ! [[ "$i" =~ tar.gz$ ]]; then 
        NAME=$(echo "$i" | sed 's/.gz//g')
        gunzip ${i}
        tar -czf ${NAME}.tar.gz ${NAME}
    else
        tar -czf ${i}.tar.gz ${i}
    fi
done