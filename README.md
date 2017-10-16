# MeTAfisher
Program to retrieve toxin antitoxin (TA) systems in metagenomes 

[comment]: <> (General intro to the program)


## Files format requierement
Metafisher take 4 different files as an input. The files need be in the same folder and to have the same name with only a specif  extension for each.
### Fasta files
Metafisher requires 3 files in the fasta files format: 
1. The genome or metagenome sequences with the extension **.fasta**  
2. The protein sequences in amino acide with the extension **.faa**
3. The protein sequences in nucleotide with the extension **.fna**

For the genome or metagenome sequences the headers of the sequence need to start with the name of the genome or metagenome follow by a space or a new line    
**Example:**   
 ```>CP000697.1 Acidiphilium cryptum JF-5, complete genome ``` here the name of the genome is CP000697.1 and it is followed by a white space the other information are ignore by the program.  
 ```>ICM0007MP0313_1000001```  here the name of the metagenome is ICM0007MP0313_1000001 and follwed by a new line  
For the protein sequences in amino acid and nucleotide (respectively .faa and .fna) there are currently two different possibile format for the header.    
*one typically found on ncbi:
 ```>lcl|CP000697.1_prot_Acry_0001_1 [gene=Acry_0001] [protein..... ``` with the  ```>lcl|``` followed by the name of the genome/metagenome and the id name of the gene (here  ```_prot_Acry_0001_1```), in this id name there is the protein number 0001 wich is used by the program to identify the genes.   

*an other one simpler: ```>ICM0007MP0313_1000001|5 ``` with the name of the sequence analysed followed by an | and the protein number used by the program to identify the genes. 

#### GFF file
Gff file: The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines.

2 possibilities:

CP000569.1	Genbank	gene	13754	15256	.	+	.	ID=gene12;Name=murE;gbkey=Gene;gene=murE;gene_biotype=protein_coding;locus_tag=APL_0013 

need to have everything 


REF 
http://www.ensembl.org/info/website/upload/gff.html
