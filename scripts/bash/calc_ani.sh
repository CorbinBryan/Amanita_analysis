#!/bin/bash

mkdir tmp 

cd tmp 

for i in $(ls ../export/*.tar.gz); do
    cp ${i} ./
done

rm *AV* 

rm 40019_S37_L007_contigs*

rm 40225_S35_L007_contigs*

rm SRR*

rm *.sh

rm 23NC21_S36_L007_contigs*

rm 30280_S34_L007_contigs*

rm Volvo1_AssemblyScaffolds_Repeatmasked.fasta.tar.gz

rm Amagr* Amapyr*

mv GCA_000827485.1_Amanita_muscaria_Koide_BX008_v1.0_genomic.fna.tar.gz ..

rm GCA_*

mv ../GCA_000827485.1_Amanita_muscaria_Koide_BX008_v1.0_genomic.fna.tar.gz ./


for i in $(ls); do
    tar -xzf ${i}
done

rm *.tar.gz


#ls > q_list.txt

#s > ref_list.txt

#astANI --ql q_list.txt --rl ref_list.txt -o ../ani_output.tsv --fragLen 1000

skani triangle ./* > ../skani_dist_matrix.txt

skani triangle ./* -E > ../skani_ani_edge_list.txt

cd ..

rm -r tmp