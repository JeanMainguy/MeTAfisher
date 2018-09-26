#!/bin/bash

##### TEMPLATE for MeTAfisher #####

##### You may change the following lines ###

  general_output_pathway=./output # the location of the result folder specific to the sequence analysed will be created
  meTAfisher_pathway=.  # the path to the folder of MeTAfisher containing dependence and program folders (it is the current directory if you launch this script from there)

  MetaG_name=Desulfovibrio_vulgaris_DP4_test # Name of the chromosome or Metagenome. Result files will be named after it
  data_pathway=./data/Desulfovibrio_vulgaris_DP4 # path to the data where the data files are stored
  data_name=sequence #common name of the data files (here we have sequence.faa sequence.fna sequence.gff sequence.fasta in the data folder)

#######

#### The next lines don't need to be be changed  ####
dependency_pathway=${meTAfisher_pathway}/dependence

# Creation of the output folders if they don't exist
output_pathway=${general_output_pathway}/${MetaG_name}
echo creation of folder ${output_pathway} if it does not exist
mkdir ${general_output_pathway}  2> /dev/null
mkdir ${output_pathway}  2> /dev/null

#HMMSEARCH parameters below should normally not be changed
#HMM result table python script expects : table_hmm = output_way + '/' + metaG_name + '_output_tableHMM.txt'
table_hmm=${output_pathway}/${MetaG_name}_output_tableHMM.txt
HMM_DB=${dependency_pathway}/ALL_plus_MET_curatted.hmm
faa_data=${data_pathway}/${data_name}.faa
echo ${table_hmm}

# Test if table HMM already exist
if [ -e ${table_hmm} ]
then
    echo "Hmm table already exist"
else
    echo "Hmmsearch"
    hmmsearch -E 0.5 --domtblout ${table_hmm}  ${HMM_DB} ${faa_data} > /dev/null
fi


# usage: MeTAfisher [-h] [--Resize] [--Rescue] [--contig_name CONTIG_NAME]
#                   [--HMM_db HMM_DB]
#                   MetaG_name output_pathway data_pathway data_name
#                   dependency_pathway

echo "MeTAfisher"
python2 ${meTAfisher_pathway}/program/main_MetaF.py  ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
# python program/main_MetaF.py ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
