# Samples Analyzed
* Downloaded all available NCBI genomes 
* Downloaded SRA for 11 novel chinese taxa (there were 25 in total but genomes were poor quality)
* Downloaded *Amanita* aff. *grandis* and *A.* aff. *conicoverrucosa* from Mycocosm
* Additionally: included 25 *A. muscaria* genomes, two of which correspond to *A. muscaria* var. *flavivolvata* 
* Additionally, I downloaded the *Volvariella volvaciens* genome to serve as an outgroup
* NCBI metadata information is stored in TSV format as `ncbi_dataset.tsv` in this directory. 
* Note that each genome downloaded from NCBI has it's own directory in `./data_set/ncbi_genomes`. This is also the case for the JGI genomes, which are stored in the same directory. 
* I transferred all files to the CHTC computing cluster (stored in a directory in staging).
* The file suffixes are not uniform--I have opted to write a script to change the names which I will transfer to the cluster as well. Just in case I screwed up the original script, I got AI to rewrite it for me. 
* I also wrote a script to retrieve each assembled (contigs and scaffolds) SRA genome and rename it according to the directory it was in. This was an oversight on my part--next time I would recommend renaming during the assembly script to prevent having to do this in the future. 
* I don't believe that I said this above, but I did assemble all of the SRA genomes using SPAdes. My script is available in `./chtc_stuff`.
* 