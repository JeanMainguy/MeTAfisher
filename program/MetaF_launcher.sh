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

#Alkalilimnicola_ehrlichii_MLHE-1
MetaG_name=Alkalilimnicola_ehrlichii_MLHE-1
data_pathway=./data/Alkalilimnicola_ehrlichii_MLHE-1
data_name=sequence

#Acaryochloris_marina_MBIC11017
MetaG_name=Acaryochloris_marina_MBIC11017
data_pathway=./data/Acaryochloris_marina_MBIC11017
data_name=sequence

# #Acaryochloris_marina_MBIC11017 FROM GENBANK directly
# MetaG_name=Acaryochloris_marina_MBIC11017_test
# data_pathway=./data/Acaryochloris_marina_MBIC11017_test
# data_name=sequence
# gb_file=data/Acaryochloris_marina_MBIC11017_test/sequence.gb
# python program/genbank_parser.py ${gb_file}

#data/Nitrobacter_hamburgensis_X14_plsm3
MetaG_name=Nitrobacter_hamburgensis_X14_plsm3
data_pathway=./data/Nitrobacter_hamburgensis_X14_plsm3
data_name=sequence

#data/Bordetella_avium_197N
MetaG_name=Bordetella_avium_197N
data_pathway=./data/Bordetella_avium_197N
data_name=sequence

output_pathway=${general_output_pathway}/${MetaG_name}
echo creation of folder ${output_pathway}
mkdir ${general_output_pathway}  2> /dev/null
mkdir ${output_pathway}  2> /dev/null

#HMMSEARCH parameter below should normally not be changed
#HMM result table python script expects : table_hmm = output_way + '/' + metaG_name + '_output_tableHMM.txt'
table_hmm=${output_pathway}/${MetaG_name}_output_tableHMM.txt
HMM_DB=${dependency_pathway}/ALL_plus_MET_curatted.hmm
faa_data=${data_pathway}/${data_name}.faa
echo ${table_hmm}
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
python program/main_MetaF.py --Resize ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
# python program/main_MetaF.py ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
