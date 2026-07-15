# BUSCO analysis
* BUSCO was run on the CHTC at UW Madison. Details on the particular docker image used are available in `./busco.sub`. This file also contains the runtime parameters, including disk, CPU, and RAM usage.
* The script, `./busco.sh` contains the exact script used to run BUSCO at run-time.
* The provided script is reproducible in the docker image specified in `./busco_template/busco.sub`.
# Gene-tree analysis
* FASTAs were filtered for those containing 75% of taxa using an in-house script. 
* Alignment, trimming, and gene-tree inference was carried out as specied in `./phylo/phylo.sh`. The file `./phylo.sun` contains the docker image used. 
* MAFFT was run with `--auto` to automatically select the most appropriate alignment algorithm and `--adjustdirection` enabled. 
* trimAl was run with automated algorithm selection, `-automated`-. 
* IQTree was run with otherwise default settings. 
* The entire script is reproducible as is in the docker image specified in `./phylo/phylo.sub`
# Species tree analysis:
* filtered, concated alignments from which gene-trees were made were used to generate a species tree as such:
```sh 
docker run --rm -it -v $(pwd):$(pwd) -w $(pwd) /root/ASTAL/bin/astral4 -i concat_gene.trees -o astral-species.tree
```
# ICA analysis
* ICA was determined using RAxML-HPC: 
```sh
raxmlHPC -f i -t astral-species.tree -z concat_gene.trees -m GTRCAT -n ICA.tree
```