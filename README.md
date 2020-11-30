<!--
 # Table of contents
- [Program overview](#Program_overview)
- [Dependence files](#Dependence_files)
- [Files format requirement](#Files_format_requirement)
- [Output file](#Output_file)
- [How to use MeTAfisher](#How_to_use_MeTAfisher) -->


# Overview
Program to retrieve toxin antitoxin (TA) systems in genomes and metagenomes

# Quick start

1. Clone the metafisher repository: `git clone https://github.com/JeanMainguy/MeTAfisher`

2. Create an environement with python3 and HMMER install. You can create an environment using conda.

3. Run metafisher on the data test:

```bash



```


# Licence
This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bionitio-team/bionitio/master/LICENSE).



<a name="Dependence_files"/>
## Dependence files
MeTAfisher requires specific files to work. These files need to be placed in the same folder and the path to this folder is given as an argument to the program.
The dependences folder has to have:

* `distance.csv` and `length_TA.csv`: used to build the score
* `ALL_plus_MET_curatted` : the database of TA HMM profile.
* `domaines_METAfisher.csv` : a table with information about the profiles of the data base (i.e TA family, description, source..)

Additionally MeTAfisher needs HMMER to be installed as the program retrieve TA systems according hmmsearch output. And obviously it requires python 2.7.

<a name="Files_format_requirement"/>

## Files format requirement
Metafisher takes 4 different files as an input. The files need be in the same folder and they need to have the same name with only a specific extension for each and every one of them.

### Fasta files
Metafisher requires 3 files in the fasta files format:

1. The chromosome or scaffold sequences with the extension **.fasta**  
2. The protein sequences in amino acide with the extension **.faa**
3. The protein sequences in nucleotide with the extension **.fna**

For the chromosome or scaffold sequences the headers of the sequences need to start with the name of the chromosome or scaffold followed by a space or a new line.    


 `>CP000697.1 Acidiphilium cryptum JF-5, complete genome` Here the name of the chromosome is `CP000697.1` and it is followed by a white space, the other information are ignored by the program.    ```>ICM0007MP0313_1000001```  Here the name of the metagenome is `ICM0007MP0313_1000001` and it is followed by a new line.



For the protein sequences in amino acid and nucleotide (respectively .faa and .fna) the header of the sequence should have the name of the chrm/scaffold followed by an | and the gene number which is used by the program to identify the genes.\\


<!-- One typically found on NCBI:

 `>lcl|CP000697.1_prot_Acry_0001_1 [gene=Acry_0001] [protein.....` with the  ```>lcl|``` followed by the name of the chromosome or scaffold  and the id name of the gene (here ```_prot_Acry_0001_1```), in this id name there is the protein number 0001 wich is used by the program to identify the genes.
 -->
Exemple : \\
 `>ICM0007MP0313_1000001|5`

<a name="gff_file"/>

### GFF file
MeTAfisher needs a gff file containing information of the genes under analysed.
> The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines.

The columns used by MeTAfisher:
* 1. **seqname** - name of the chromosome or scaffold
* 3. **feature** - feature type name, e.g. Gene, Variation, Similarity
* 4. **start** - Start position of the feature, with sequence numbering starting at 1.
* 5. **end** - End position of the feature, with sequence numbering starting at 1.
* 7. **strand** - defined as + (forward) or - (reverse).
* 8. **frame** - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
* 9. **attribute** - A semicolon-separated list of tag-value pairs, providing additional information about each feature.   


The program looks only at the line with feature equal to `CDS`.
 MeTAfisher needs to match the protein hit by hmmsearch with the proper gff line, and to do so the protein number should be identical in the gff file and in the header of the fasta files.
 In order to do that, the attribute section (columns 9) of the gff file should follow the following pattern:
`ID=< chromosome or scaffold>|<protein number>`.

Example of correct gff line :    ```ICM0007MP0313_1000310	GPF	CDS	13787	14128	.	-	0	ID=ICM0007MP0313_1000310|19;partial=00;sta.....```

<a name="Output_file"/>

## Output file

MeTAfisher provides different output file: 4 result files showing TA systems found and a csv file gathering statistique information of the analysis is also provided.
All the output files start with 5 lines starting with `##`.

```
## result_S
## Rescue lonely gene : False
## Resize gene : True
## Distance threshold from -100nt to 300nt
## Length threshold from 30aa to 1500aa
```

The first line tell us the kind of output file it is. Then the two following lines tell if the Rescue and Resize have been used and finnally the last two lines tell which distance and length threshold have been used during the analysis.

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

* Short result:
* Human readable result
* Table csv-like result


### Short result
One line for each systems found.
One line is composed of the two id of the genes, the strand and finally the score of the system.
This output is just there to give an overview of the result and to get the gene id.  

### Table csv-like result:
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

<a name="How_to_use_MeTAfisher"/>

## How to use MeTAfisher

First of all you should download MeTAfisher files and folders from gitHub.

### File format challenges
The file format requirement of MeTAfisher are quite heavy and may require a lot of effort if you start with some files that do not follow the required pattern. However if the sequence you want to analyze is on genbank, you can use the genbank file of the sequence and process it through the script genbank_parser.py. This script takes as argument the genbank file and creates 3 files needed for the TA analysis (.faa, fna and .gff) in the same folder as the genbank file. `genbank_parser.py` is only able to parse a genbank file with all sequences. When downloading the file on the genbank web page, you should be certain that the sequence of the genes are displayed. In order to do so, you may click in the "Customize" view left panel on the option "Show sequence" and then on "Update view". Finally, to download the file, click on the upper left "Send to" button > Complete Record > Choose Destination: File > Format: Genbank > Create File.
Then be sure to place the genbank file in a known and specific directory. Then you can launch `genbank_parser.py`:

`python <path to the script>/genbank_parser.py <path and name of the genbank file>`

So for example if you are at the root of MeTAfisher project and you want to parse a genbank file placed in `./data/Desulfovibrio_vulgaris_DP4/` you should enter the following command:

`python program/genbank_parser.py data/Desulfovibrio_vulgaris_DP4/sequence.gb`

### Launching the program
Once you have all the files required and correctly shaped you may want to launch MeTAfisher.
To launch MeTAfisher you can either launch it directly with the script `main_MetaF.py` and provide correctly all the required argument in the command line - or you may also use the MetaF_launcher.sh bash script as a template.
You need to provide different information in the script:
* general_output_pathway: the folder where the result folder of the sequence analyzed will be created.
* dependency_pathway: the pathway of the dependence folder where the dependences of MeTAfisher are stored.
* Sequence_name: corresponds to the name of the Metagenome or Genome analyzed
* data_pathway: Path to the data folder where fasta and gff files are stored.
* data_name: the common name without the extension of the data files. the 3 fasta files and the gff file should be in the same folder and with a unique name, only the file's extension are different.  

Then, the script takes care of the rest: it creates a folder named after the Sequence name given in the general output folder provided. Then, it launches hmmsearch which searches TA domains in the amino acid sequence. If the hmmsearch table output already exists, it skips this step in order not to launch twice the same thing. And finally it launches the python script `main_MetaF.py` with the required arguments. If you want to activate the "Resize" and/or the "Rescue" step, you can then add the flag --Resize or --Rescue in the command line.

In the output folder that you gave to the program you will find the hmmsearch table result and the different output files.




## Score



REF
http://www.ensembl.org/info/website/upload/gff.html
