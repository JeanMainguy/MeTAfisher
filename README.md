# MeTAfisher
Program to retrieve toxin antitoxin (TA) systems in metagenomes 

[comment]: <> (General intro to the program)


## Files format requierement
### Fasta files
Metafisher requires 3 files in the fasta files format: 
#### the genome or metagenome sequences with the extension **.fasta**
The header of this file need to start with the name of the genome or metagenome follow by a space or line returned
example: 
">CP000697.1 Acidiphilium cryptum JF-5, complete genome" in this examaple the name of the genome is CP000697.1 and it is followed by a white space.
">ICM0007MP0313_1000001" and here the name of the metagenome is ICM0007MP0313_1000001 and follwed by a new line


>CP000697.1 Acidiphilium cryptum JF-5, complete genome
and here the name of the metagenome is ICM0007MP0313_1000001 and follwed by a new line
>ICM0007MP0313_1000001

For the files containing amino acid and nucleotid sequences (respectively .faa and .fna)
2 differents format:
one typically found on ncbi:
>lcl|CP000697.1_prot_Acry_0001_1 [gene=Acry_0001] [protein.....
and an other one simpler: 
>ICM0007MP0313_1000001|5 
with the name of the sequence analysed followed by an | and the gene number

Gff file: The GFF (General Feature Format) format consists of one line per feature, each containing 9 columns of data, plus optional track definition lines. 

2 possibilities:

CP000569.1	Genbank	gene	13754	15256	.	+	.	ID=gene12;Name=murE;gbkey=Gene;gene=murE;gene_biotype=protein_coding;locus_tag=APL_0013 

need to have everything 


REF 
http://www.ensembl.org/info/website/upload/gff.html
