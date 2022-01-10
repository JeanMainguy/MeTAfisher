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

## Output files



## Scoring of TA systems
