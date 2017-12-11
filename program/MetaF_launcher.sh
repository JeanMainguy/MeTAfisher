#!/bin/bash

general_output_pathway=./output
dependency_pathway=./dependence

# MP0313
MetaG_name=MP0313
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/Malaspina_MetaG/MP0313
data_name=MP0313.549

# Acry Genome
MetaG_name=Acry
data_pathway=./data/Acry_Genome
data_name=sequence

#Magnetococcus_marinus
MetaG_name=Magnetococcus_marinus
data_pathway=./data/Magnetococcus_marinus
data_name=sequence


output_pathway=${general_output_pathway}/${MetaG_name}
mkdir ${output_pathway}  2> /dev/null


#HMMSEARCH parameter below should normally not be changed
#HMM result table python script expects : table_hmm = output_way + '/' + metaG_name + '_output_tableHMM.txt'
table_hmm=${output_pathway}/${MetaG_name}_output_tableHMM.txt
HMM_DB=${dependency_pathway}/ALL_plus_MET_curatted.hmm
faa_data=${data_pathway}/${data_name}.faa
#!/bin/bash
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
python program/main_MetaF.py --Rescue ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
# python program/main_MetaF.py ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
