## Samples Analyzed
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
## Genome assembly of SRA files
* SRA files were assembled using SPAdes, with parameters specified in `./chtc_stuff/asm`. This particular script uses my Docker image for SPAdes, which also includes a couple of other tools for interacting with FASTQ files. The image is `multimuscaria/fastqtk:4.0`. 
* Although both the contigs and scaffolds were saved, only the contigs were used for downstream phylogenomic inference.
## BUSCO analysis
* BUSCO was run using an older version of the program packaged in a docker image created by a collaborator (`aflatoxing/maker3`). 
* Note that several changes had to be made to previous versions of this script for it run. The CHTC has once again changed how they treat Docker environments--for some reason, editing files in the Docker image is no longer permitted unless you copy them into the scratch directory made at run-time. 
* I've included the updated script and submission file for running BUSCO on the cluster (`./chtc_stuff/busco_template`). Note that I removed the syntax for transferring files as this directory is public. 
* An additional error--when making a tar ball, it seems the cluster prefers both the output and input arguments to be surrounded in double quotes. Accordingly, variables should be surrounded with curly brackets. This is also just good practice for clarity and disambiguation. 