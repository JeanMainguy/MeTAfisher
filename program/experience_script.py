import argparse
from subprocess import call
import Function_MetaF as fct2
import Object_MetaF as obj

obj.Gene.output_way = "/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/output/MP0313/test"
obj.Gene.hmmdb = "/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/dependence/ALL_plus_MET_curatted.hmm"
faa_file = "/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher/data/MP0313/MP0313.549.faa"

fct2.HMM_launcher(faa_file, add_to_name='test')

# home = "/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/MeTAfisher"
#
#
# bash_commande = "hmmsearch -E 0.5  --domtblout {}/output/MP0313/test/adjOrf_table_hmm.txt  {}/dependence/ALL_plus_MET_curatted.hmm {}/data/MP0313/MP0313.549.faa".format(home, home, home)
# print bash_commande
#
# call(bash_commande, shell=True)
