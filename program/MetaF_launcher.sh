#!/bin/bash

# HMMSEARCH


# PYTHON script
# MP0313
MetaG_name=MP0313
output_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/output/MP0313
data_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/data/MP0313
data_name=MP0313.549
dependency_pathway=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/dependence


# usage: MeTAfisher [-h] [--Resize] [--Rescue] [--contig_name CONTIG_NAME]
#                   [--HMM_db HMM_DB]
#                   MetaG_name output_pathway data_pathway data_name
#                   dependency_pathway


python main_MetaF.py --Resize --Rescue ${MetaG_name} ${output_pathway} ${data_pathway} ${data_name} ${dependency_pathway}
