#!/bin/bash

general_output_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/output

# MP0313
MetaG_name=MP0313
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/Malaspina_MetaG/MP0313
data_name=MP0313.549

# Acry Genome
MetaG_name=Acry
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/data/Acry_Genome
data_name=sequence

#Magnetococcus_marinus
MetaG_name=Magnetococcus_marinus
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/data/Magnetococcus_marinus
data_name=sequence

dependency_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/dependence
output_pathway=${general_output_pathway}/${MetaG_name}
mkdir ${output_pathway}  2> /dev/null

# is it necessary to reformat the input file?
reformat_path=${output_pathway}/reformated_files
mkdir ${reformat_path}  2> /dev/null
python program/format_check.py ${data_pathway} ${data_name} ${reformat_path}
if [[ $? = 0 ]]; then
    echo "Reformating is a success"
    #symbolic link of the chromosome/contig seq beause it is not changed by the refomarting step
    ln -s ${data_pathway}/${data_name}.fasta ${reformat_path}/${data_name}.fasta  2> /dev/null
    #We then replace the data_pathway with the patway of the reformated files
    data_pathway=${reformat_path}

else
    echo "Reformating has failed: $?"
    exit 1
fi

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
