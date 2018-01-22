# MeTAfisher
Program to retrieve toxin antitoxin (TA) systems in metagenomes

[comment]: <> (General intro to the program)

## Program overview
MeTAfisher is written in Python 2.7. It is made up of four files which need to be in the same folder.
* `main_MetaF.py` : It initializes the variable and launches the different functions of the program inside a loop iterating the different contigs. This file takes argument to be able to correctly initialyze the variable.
* `Function_MetaF.py` : All the general functions of the program
* `Object_MetaF.py` :  Classes and methods.
* `Orf_MetaF.py` : Functions specific to ORF process, used when the option `--Rescue` is called.


`main_MetaF.py` needs argument to run the program :
blballablabalb


## Dependence files
MeTAfisher requires specific files to work. These files need to be placed in the same folder and the path to this folder is given as an argument to the program.
The dependences folder has to have:

* `distance.csv` and `length_TA.csv`: used to build the score
* `ALL_plus_MET_curatted` : the database of TA HMM profile.
* `domaines_METAfisher.csv` : a table with information about the profiles of the data base (i.e TA family, description, source..)

## Files format requierement
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

### GFF file
> The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines.

The columns used by MeTAfisher:
1. **seqname** - name of the chromosome or scaffold
3. **feature** - feature type name, e.g. Gene, Variation, Similarity
4. **start** - Start position of the feature, with sequence numbering starting at 1.
5. **end** - End position of the feature, with sequence numbering starting at 1.
7. **strand** - defined as + (forward) or - (reverse).
8. **frame** - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
9. **attribute** - A semicolon-separated list of tag-value pairs, providing additional information about each feature.   


The program looks only at the line with feature equal to `CDS`.
 MeTAfisher needs to retrieve the protein number to be able to match the protein hit by hmmsearch with the proper gff line. 
 In order to do that, the attribute section should follow the following pattern:
1. One with `ID=< chromosome or scaffold>|<protein number>`.

Example of correct gff line :    ```ICM0007MP0313_1000310	GPF	CDS	13787	14128	.	-	0	ID=ICM0007MP0313_1000310|19;partial=00;sta.....```

 <!-- CP000569.1	Genbank	gene	13754	15256	.	+	.	ID=gene12;Name=murE;gbkey=Gene;gene=murE;gene_biotype=protein_coding;locus_tag=APL_0013 -->

## Output file
Diﬀerent output ﬁles currently exist.
* Short output :
* Human readable file :
* Table csv like file:

The program also provides a file gathering quantitative information about the analysis of the contigs/chromosomes analyzed.    

## How to use MeTAfisher

### File format challenges
The file format requirement of MeTAfisher are quite heavy and may require a lot of effort if you start with some files that do not follow the required pattern. However if the sequence you want to analyze is on genbank, you can use the genbank file of the sequence and process it through the script genbank_parser.py. This script takes as argument the genbank file and creates 3 files needed for the TA analysis (.faa, fna and .gff) in the same folder as the genbank file. Genbank_parser.py is only able to process the full genbank file. When downloading the file on the genbank web page, you should be certain that the sequence of the genes are displayed. In order to do so, you may click in the "Customize" view left panel on the option "Show sequence" and then on "Update view". Finally, to download the file, click on the upper left "Send to" button > Complete Record > Choose Destination: File > Format: Genbank > Create File. 

### Launching the program
Once you have all the files required and correctly shaped you may want to launch MeTAfisher. 
To launch MeTAfisher you can either launch it directly with the script `main_MetaF.py`and provide correctly all the required argument in the command line - or you may also use the MetaF_launcher.sh bash script as a template. 
You need to provide different information in the script:
* general_output_pathway: the folder where the result folder of the sequence analyzed will be created. 
* dependency_pathway: the pathway of the dependence folder where the dependences of MeTAfisher are stored. 
* Sequence_name: 
* data_pathway:
* data_name:

Then, the script takes care of the rest: it creates a folder named after the Sequence name given in the general output folder provided. Then, it launches hmmsearch which searches TA domains in the amino acid sequence. If the hmmsearch table output already exists, it skip this step in order not to launch twice the same thing. And finally it launches the python script main_MetaF.py with the required argument. If you want to activate the "Resize" and/or the "Rescue" step, you can then add the flag --Resize or --Rescue in the command line. 








REF
http://www.ensembl.org/info/website/upload/gff.html
