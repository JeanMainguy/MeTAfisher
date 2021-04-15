<!--
 # Table of contents
- [Program overview](#Program_overview)
- [Dependence files](#Dependence_files)
- [Files format requirement](#Files_format_requirement)
- [Output file](#Output_file)
- [How to use MeTAfisher](#How_to_use_MeTAfisher) -->


# Metafisher
Metafisher is a tool to retrieve toxin antitoxin (TA) systems of type II in genomic sequences.


## Quick start

1. Clone the metafisher repository:
```bash
git clone https://github.com/JeanMainguy/MeTAfisher
cd MeTAfisher/
```

2. Create an environement with python3 and HMMER installed. You can create this environment using conda:

```bash
conda env create -f env/metafisher.yml
conda activate metafisher
```

3. Run metafisher on the genome of Desulfovibrio vulgaris:

```bash
./metafisher/metafisher.py --gff data_test/GCF_000070465.1/GCF_000070465.1_ASM7046v1_genomic.gff.gz \
                          --faa data_test/GCF_000070465.1/GCF_000070465.1_ASM7046v1_protein.faa.gz \
                          --outdir metafisher_results -v

```


## Licence
This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bionitio-team/bionitio/master/LICENSE).


## Identification of TA genes

To identify potential Toxin and Antitoxin genes, metAfisher uses a list of domains known to be specific of TA systems. These domains are searched in the protein sequences by the tool HHMER.
On top of the domain search, potential TA genes can be identified by diamond search based on all sequences of TADB (https://bioinfo-mml.sjtu.edu.cn/TADB2).

### Create diamond database

To use diamond search strategy, a diamond database with the TADB sequences need to be created.

The protein sequences of Toxin and Antitoxin can be downloaded on the TADB website: https://bioinfo-mml.sjtu.edu.cn/TADB2/download.html

1. Download TADB protein sequences

```bash
wget https://bioinfo-mml.sjtu.edu.cn/TADB2/download/TADB2/20171013/protein/type_II_pro_T.fas
wget https://bioinfo-mml.sjtu.edu.cn/TADB2/download/TADB2/20171013/protein/type_II_pro_AT.fas
```

2. concat fasta files and build diamond db

```bash

mkdir TA_data
cat type_II_pro_T.fas type_II_pro_AT.fas > TA_data/type_II_TA.fasta

diamond makedb --in type_II_TA.fasta -d TA_data/type_II_TA

```

3. Generate stat files

These file are needed to score the potential TA systems. It computes how often a domain is associated with another one in a TA system of TADB.   

```bash
python metafisher/compute_tadb_stat.py --toxin_faa TA_data/type_II_pro_T.fas --antitoxin_faa TA_data/type_II_pro_AT.fas -v
```

### Launch MeTAfisher with diamond search

```bash

./metafisher/metafisher.py --gff data_test/GCF_000070465.1/GCF_000070465.1_ASM7046v1_genomic.gff.gz \
                         --faa data_test/GCF_000070465.1/GCF_000070465.1_ASM7046v1_protein.faa.gz\
                         --outdir metafisher_results \
                         --diamond_db TA_data/type_II_TA.dmnd -v

```

## Scoring of TA systems



## Output files

MeTAfisher provides different output file: 4 result files showing TA systems found and a csv file gathering statistics information of the analysis is also provided.
All the output files start with 5 lines starting with `##`.

```

## result_S
## Rescue lonely gene : False
## Distance threshold from -100nt to 300nt
## Length threshold from 30aa to 1500aa
```

The first line tell us the kind of output file it is. Then the next line indicates if the `rescue` mode have been used and finnally the last two lines tell which distance and length threshold have been used during the analysis.

The name of the output file follow a regular pattern, they all start with the name of the metagenome or genome provided when launching the program follow by the name of the output file and the appropriate extension.

For genome Desulfovibrio vulgaris DP4 we have:

The result files:

* Desulfovibrio_vulgaris_DP4_result_H.txt
* Desulfovibrio_vulgaris_DP4_result_S.txt
* Desulfovibrio_vulgaris_DP4_result_T.csv
* Desulfovibrio_vulgaris_DP4_result_GFF.gff

The stat file:

* Desulfovibrio_vulgaris_DP4_contig_stat.csv


### Result files
Four result ﬁles currently exist, they provide the result of the analysis in a different way.

* Short result
* Human readable result
* Tabular result
* Gff file


### Short result
One line for each systems found.
One line is composed of the two id of the genes, the strand and finally the score of the system.
This output is just there to give an overview of the result and to get the gene id.  

### Tabular result:
This result file may be open on a spreadsheet.
Each TA gene is displayed on a line. Information are displayed within 9 columns:

1. Contig: contig or chromosome names.
2. Gene number: number used by the program to identify gene. The same one found in the input gff and faa/fna files.
3. Gene id: The genes id used by genbank when the input files are made from a genbank file
4. start: postion of the start
5. end: end position
6. length
7. length_score: not display yet
8. strand: defined as + (forward) or - (reverse).
9. feature: CDS or ORF. if ORF the gene is from a home-made prediction done during the Rescue step
10. domain; info about the best domain of the gene
11. Neighbor gene: info about the neighbor genes. Down if the it is a downstream neigbor and Up if it is a upstream neighbor


The program also provides a file gathering quantitative information about the analysis of the contigs/chromosomes analyzed.
