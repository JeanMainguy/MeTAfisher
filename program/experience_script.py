import argparse
from subprocess import call
import Function_MetaF as fct2
import Object_MetaF as obj
# parser = argparse.ArgumentParser(prog='MeTAfisher', description='Identification of Toxin Antitoxin Systems', epilog="ADDITIONAL INFORMATION:")
# parser.add_argument("MetaG_name", help="Name of the Metagenome or Genome")
# parser.add_argument("table_hmm", help="Table result output of --domtblout of the hmmsearch")
# parser.add_argument("output_pathway", help="Pathway of the output folder")
# parser.add_argument("data_pathway", help="Pathway of the data folder where all the data files are stored : fna_file, faa_file, scaffold_file, gff_file")
# parser.add_argument("data_name", help="Common name of all data files, the name without the extension.  fna_file, faa_file, scaffold_file, gff_file are different only by their extansions : .fna, .faa, .fasta, .gff respectively")
# parser.add_argument("dependency_pathway", help="Pathway of the dependency folder where all te depedencies are stored")
# parser.add_argument('--Resize', dest='resize', action='store_true',
#                     help="Resize the genes if they are too big for the thresholds and take into account the possible start along the sequence. To do only if the gene prediction is not trustable")
# parser.add_argument('--Rescue', dest='rescue', action='store_true',
#                     help='To do the rescue step of lonely genes, by default it is False')
# parser.add_argument('--contig_name', default=False, help='Name of a specific contig to analysed. The program will analysed only this conitg')
#
#
# # parser.add_argument('--home_path', default='./', help='pathway of MeTAfisher program')
# # parser.add_argument("--fna_file", help="Multifasta file of predicted genes sequence in nucleotide")
# # parser.add_argument("--faa_file", help="Multifasta file of predicted genes sequence in amino acid")
# # parser.add_argument("--scaffold_file", help="File of the complete sequence of the contigs or of the genome")  # OPTIONEL SEULEMENT QUAND RESCUE IS ON
# # parser.add_argument("--gff_file", help="Information file of the gene in a gff file format")
#
# args = parser.parse_args()
# print 'rescue ', args.rescue
# print 'resize ', args.resize
# print args.MetaG_name
# print "contig name", args.contig_name
# if args.contig_name:
#     print 'args.contig_name is not False'
# else:
#     print "args.contig_name is FALSE"
# # print args.accumulate(args.integers)
#
# # metaG_name
# # input_way
# # gff_file
# # scaffold_file
# # fna_file
# # table_hmm

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
