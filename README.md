# MeTAfisher
Program to retrieve toxin antitoxin (TA) systems in metagenomes

[comment]: <> (General intro to the program)

## Program overview
MeTAfisher is written in Python 2.7. It is made of four files which need to be in the same folder.
* `main_MetaF.py` : It initialies the variable and launchs the different functions of the program inside a loop iterating the different contigs. This file takes argument to be able to correctly initilyse the variable
* `Function_MetaF.py` : All the general function of the program
* `Object_MetaF.py` :  Classes and methods.
* `Orf_MetaF.py` : Fucntions specific to ORF process, used when the option `--Rescue` is called.


`main_MetaF.py` needs argument to run the program :
blballablabalb


## Dependence files
MeTAfisher requires specific files to work. These files need to be place in the same folder and the path to this folder given as an argument to the program.
The dependeces folder has to have:

* `distance.csv` and `length_TA.csv`: used to build the score
* `ALL_plus_MET_curatted` : the database of TA HMM profile.
* `domaines_METAfisher.csv` : a table with information about the profiles of the data base (i.e TA family, description, source..)

## Files format requierement
Metafisher take 4 different files as an input. The files need be in the same folder and to have the same name with only a specif  extension for each.
### Fasta files
Metafisher requires 3 files in the fasta files format:
1. The chromosome or scaffold sequences with the extension **.fasta**  
2. The protein sequences in amino acide with the extension **.faa**
3. The protein sequences in nucleotide with the extension **.fna**

For the chromosome or scaffold sequences the headers of the sequence need to start with the name of the  chromosome or scaffold follow by a space or a new line.    


 `>CP000697.1 Acidiphilium cryptum JF-5, complete genome` here the name of the chromosome is `CP000697.1` and it is followed by a white space, the other information are ignored by the program.    ```>ICM0007MP0313_1000001```  here the name of the metagenome is `ICM0007MP0313_1000001` and followed by a new line



For the protein sequences in amino acid and nucleotide (respectively .faa and .fna) the header of the sequence should have the name of the chrm/scaffold followed by an | and the gene number used by the program to identify the genes.\\


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
 MeTAfisher needs to retrieve the protein number to be able to match the protein hited by hmmsearch with the proper gff line.  To do so the attribute section needs to follow to possible formats:
1. One with `ID=< chromosome or scaffold>|<protein number>`.
2.  One other possible format is with ????

Example of correct gff line :    ```ICM0007MP0313_1000310	GPF	CDS	13787	14128	.	-	0	ID=ICM0007MP0313_1000310|19;partial=00;sta.....```

 <!-- CP000569.1	Genbank	gene	13754	15256	.	+	.	ID=gene12;Name=murE;gbkey=Gene;gene=murE;gene_biotype=protein_coding;locus_tag=APL_0013 -->

## Output file
Diﬀerent possible output ﬁles exist.
* Short output :
* Human readable file :
* ...

The program also provide a file gathering quantitative information about the analysis of the contigs/chromosomes analysed.    



REF
http://www.ensembl.org/info/website/upload/gff.html
