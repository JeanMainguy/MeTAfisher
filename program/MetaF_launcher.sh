#!/bin/bash

# HMMSEARCH

general_output_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/output
# PYTHON script
# MP0313
MetaG_name=MP0313
output_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/output/MP0313
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/Malaspina_MetaG/MP0313
data_name=MP0313.549


# Acry Genome
MetaG_name=Acry
output_pathway=${general_output_pathway}/${MetaG_name}
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/data/Acry_Genome
data_name=sequence_Fchecked

dependency_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/dependence

mkdir ${output_pathway}


#HMMSEARCH parameter below should normally not be changed
#HMM result table python script expects : table_hmm = output_way + '/' + metaG_name + '_output_tableHMM.txt'
table_hmm=${output_pathway}/${MetaG_name}_output_tableHMM.txt
HMM_DB=${dependency_pathway}/ALL_plus_MET_curatted.hmm
faa_data=${data_pathway}/${data_name}.faa

# hmmsearch -E 0.5 --domtblout ${table_hmm}  ${HMM_DB} ${faa_data} > /dev/null

# usage: MeTAfisher [-h] [--Resize] [--Rescue] [--contig_name CONTIG_NAME]
#                   [--HMM_db HMM_DB]
#                   MetaG_name output_pathway data_pathway data_name
#                   dependency_pathway


# python program/main_MetaF.py --Resize --Rescue ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
python program/main_MetaF.py ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
