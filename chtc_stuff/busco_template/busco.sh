#!/bin/sh
#
#unzipping the genome folder
set -e 
#tar -xzvf "$1".fna.tar.gz

#mv ./"$1"/"$1"_prot.faa ./
#TODO
access="$1"
#
#cp /staging/bryan7/agaricales_odb10.2020-08-05.tar.gz .
tar -xzvf agaricales_odb10.2020-08-05.tar.gz
#
cat /usr/local/bin/busco/config/config.ini | sed 's|/usr/local/bin/|/usr/bin/|'> /usr/local/bin/busco/config/temp
mv /usr/local/bin/busco/config/temp /usr/local/bin/busco/config/config.ini
#
run_BUSCO.py -i "$access".fna -o "$access"_busco -l ./agaricales_odb10 -m genome
#
tar -czf run_"$access"_busco.tar.gz run_"$access"_busco/

rm *