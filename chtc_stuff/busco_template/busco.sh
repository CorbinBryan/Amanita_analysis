#!/bin/sh
#
#unzipping the genome folder
set -e 

access="$1"

tar -xzf ${access}.tar.gz
#
tar -xzvf agaricales_odb10.2020-08-05.tar.gz
#
cp /usr/local/bin/busco/config/config.ini ./config.ini
sed -i 's|/usr/local/bin/|/usr/bin/|g' ./config.ini
export BUSCO_CONFIG_FILE="$PWD/config.ini"
#
run_BUSCO.py -i "$access" -o "$access"_busco -l ./agaricales_odb10 -m genome
#
tar -czf "run_${access}_busco.tar.gz" "run_${access}_busco/"