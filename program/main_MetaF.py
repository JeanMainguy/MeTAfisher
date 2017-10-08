# coding: utf-8
import Object_MetaF as obj
import Function_MetaF as fct2
import Orf_MetaF as orf2
import sys
import resource
import csv
import argparse


def using(point=""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           ''' % (point, usage[0], usage[1],
                  (usage[2] * resource.getpagesize()) / 1000000.0)


parser = argparse.ArgumentParser(prog='MeTAfisher', description='Identification of Toxin Antitoxin Systems', epilog="ADDITIONAL INFORMATION:")
parser.add_argument("metaG_name", help="Name of the Metagenome or Genome")
# parser.add_argument("table_hmm", help="Table result output of --domtblout of the hmmsearch")
parser.add_argument("output_pathway", help="Pathway of the output folder")
parser.add_argument("data_pathway", help="Pathway of the data folder where all the data files are stored : fna_file, faa_file, scaffold_file, gff_file")
parser.add_argument("data_name", help="Common name of all data files, the name without the extension.  fna_file, faa_file, scaffold_file, gff_file are different only by their extansions : .fna, .faa, .fasta, .gff respectively")
parser.add_argument("dependency_pathway", help="Pathway of the dependency folder where all te depedencies are stored")
parser.add_argument('--Resize', dest='resize', action='store_true',
                    help="Resize the genes if they are too big for the thresholds and take into account the possible start along the sequence. To do only if the gene prediction is not trustable")
parser.add_argument('--Rescue', dest='rescue', action='store_true',
                    help='To do the rescue step of lonely genes, by default it is False')
parser.add_argument('--contig_name', default=False, help='Name of a specific contig to analysed. The program will analysed only this conitg')
parser.add_argument("--HMM_db", default="ALL_plus_MET_curatted.hmm", help="name of the HMM database")


args = parser.parse_args()
metaG_name = args.metaG_name
output_way = args.output_pathway
data_way = args.data_pathway
data_name = args.data_name
dependence_way = args.dependency_pathway
HMM_db = dependence_way + '/' + args.HMM_db

gff_file = data_way + '/' + data_name + '.gff'
scaffold_file = data_way + '/' + data_name + '.fasta'
fna_file = data_way + '/' + data_name + '.fna'
faa_file = data_way + '/' + data_name + '.faa'
table_hmm = output_way + '/' + metaG_name + '_output_tableHMM.txt'
hmm_adjorf_launcher = output_way + '/' + data_name + '.gff'


# Flag
resize = args.resize
rescue = args.rescue

info_contig_stat = True
output_human = True
output_synthetic = True

# Storing information as Gene class attribut to be use when we launch hmmsearch
obj.Gene.output_way = output_way
obj.Gene.hmmdb = HMM_db

# DISTANCE AND LENGTH DICO :
# Dist and length of TA from TADB to mmake a proba
file_len = dependence_way + "/length_TA.csv"
file_dist = dependence_way + "/distance.csv"

# CSV FILE DOMAINS
csv_domain = dependence_way + '/domaines_METAfisher.csv'
with open(csv_domain, 'r') as csvdo:
    reader = csv.DictReader(csvdo)
    for row in reader:
        obj.Gene.domain_dict[row['hmm_name']] = {k: v for k, v in row.iteritems() if k in ['acc', 'family', 'type']}
# print obj.Gene.domain_dict
# SETTING THE THRESHOLD :
# There are the first the thrshold there very large then the proabilty score step is going to define the better conf
# treshold size of gene


lenMin = 10 * 3  # VERY IMPORTANT LENGTH HAVE TO BE GIVEN IN AA and then transform in nt
lenMax = 500 * 3

obj.Gene.length_min = lenMin
obj.Gene.length_max = lenMax

# threshold distance for tandem
distanceMin = -100
distanceMax = 300

obj.Gene.distanceMin = distanceMin
obj.Gene.distanceMax = distanceMax

# Output
# TODO put some flag here to has the choice of the ouput
file_out = '{}/{}_result_rescue_{}_resize_{}.txt'.format(output_way, metaG_name, rescue, resize)
header = "## Rescue lonely gene : {}\n".format(rescue)
header += "## Resize gene : {}\n".format(resize)
header += "## Distance threshold from {}nt to {}nt\n".format(distanceMin, distanceMax, lenMin, lenMax)
header += "## Length threshold from {}aa to {}aa\n".format(lenMin, lenMax)
if output_human:
    file_out_H = '{}/{}_result_H_rescue_{}_resize_{}.txt'.format(output_way, metaG_name, rescue, resize)
    floutH = open(file_out_H, "w")
    floutH.write("## {}: human readable output\n")
    floutH.write(header)
if output_synthetic:
    file_out_S = '{}/{}_result_S_rescue_{}_resize_{}.txt'.format(output_way, metaG_name, rescue, resize)
    floutS = open(file_out_S, "w")
    floutS.write("## {}: synthetic output\n".format(metaG_name))
    floutS.write(header)
# Output files of adj ORF :
# out_faa_file = output_way + metaG_name + '_adjacent.faa'
# out_fna_file = output_way + metaG_name + '_adjacent.fna'
# out_gff_file = output_way + metaG_name + '_adjacent.gff'

tmp_adj_orf_faa = output_way + '/temporary_adjOrf.faa'

if info_contig_stat:
    header = ['contig', 'gene with TA domain', 'lonely gene', 'linked gene']
    if rescue:
        header += ['adjacent orf', 'rescue flag', 'orf with TA domain', 'lonely gene rescue']
    fl_stat = open('{}/{}_contig_stat.csv'.format(output_way, metaG_name), 'a')
    writer_stat = csv.DictWriter(fl_stat, fieldnames=header, delimiter='\t')
    writer_stat.writeheader()
    obj.Gene.metaG_stat = dict.fromkeys(header, 0)
    obj.Gene.metaG_stat['contig'] = metaG_name


# Open step of the scaffold file ! used by fast_fasta function in orf file
# the file won't be close until the end.Each time the algo won't have to
# read the whle file again until the scaffold needed. We can do that
# because the scaffold list has been sorted
# New version use this dictionnary:
dico_orf = {}
dico_orf['fl'] = open(scaffold_file, 'r')
dico_orf['line'] = next(dico_orf['fl'])

# Set the stops and starts codon as an atribut of the class Gene
# Parameter for the orfinder

id_genetCode = 11
table = orf2.getGeneticCode(id_genetCode)
obj.Gene.codon_starts = table['start']
obj.Gene.codon_stops = table['stop']

# GFF dico use to delete the orf that correspond to a predicted gene
# The program retrieve the positions of the end of the predicted gene
# And don't proceess the orf that have the same end position

if rescue:
    gff_dico = {}
    fl_csv = open(gff_file, 'r')
    gff_dico['csv'] = csv.reader(fl_csv, delimiter='\t')
    gff_dico['line'] = next(gff_dico['csv'])

# Open file fna to retrieve sequence of the predicted gene
# PUT info in a dico to not have to retourned it every time !!
dico_seq = {}
dico_seq["fl"] = open(fna_file, 'r')
dico_seq["line"] = dico_seq["fl"].readline()
# PUT info in a dico to not have to retourned it every time !! dico is used in check_size
dico_seq["codon_start"] = table['start']
print "start", table['start']


# SCORE PREPARaTION
k = 20
dist_mltp = 5
len_mltp = 6
dico_len, Ntot_len = fct2.from_file_to_dict(file_len)
dico_dist, Ntot_dist = fct2.from_file_to_dict(file_dist)


obj.Gene.length_proba = fct2.give_proba_dict(int(lenMin / 3), int(lenMax / 3), dico_len, k, Ntot_len)
obj.Gene.distance_proba = fct2.give_proba_dict(distanceMin, distanceMax, dico_dist, k, Ntot_dist)
obj.Gene.distance_proba = fct2.transformation(obj.Gene.distance_proba, dist_mltp)
obj.Gene.length_proba = fct2.transformation(obj.Gene.length_proba, len_mltp)


# Take every scaffold present in the HMM output and sort them in order to
# be able to retrieve their sequence correctly in the input file
scaffold_list = sorted(fct2.get_list_scaffold(table_hmm))
print 'nb sca', len(scaffold_list)
# Loop : each scaffold is treated independntly here
# scaffold_list = [
#     'ICM0007MP0313_1000008', 'ICM0007MP0313_1000022', 'ICM0007MP0313_1000126',
#     'ICM0007MP0313_1000131', 'ICM0007MP0313_1000288', 'ICM0007MP0313_1000321']
# scaffold_list = sorted(['ICM0007MP0313_1000103', 'ICM0007MP0313_1000207', 'ICM0007MP0313_1000073', 'ICM0007MP0313_1000346'])
# scaffold_list = ['ICM0007MP0313_1000346']
# scaffold_list = ['ICM0007MP0313_1000085', 'ICM0007MP0313_1000086', 'ICM0007MP0313_1000089', 'ICM0007MP0313_1000408']
# scaffold_list = ['ICM0007MP0313_1000073']
for scaffold in scaffold_list:
    print '* *' * 25
    print ' * ' * 25
    print scaffold
    # header = ['contig', 'gene with TA domain', 'lonely gene', 'linked gene', 'adjacent orf', 'rescue flag', 'hmm orf', 'lonely  gene rescue']
    # Reset obj.TA_gene and Orf class attribut
    obj.TA_gene.genes = []
    obj.TA_gene.genes_strand = {'+': [], '-': []}
    obj.TA_gene.linked = set()
    obj.TA_gene.lonely = None
    obj.Orf.adj_orf = {}
    obj.Orf.adj_orf_index = 0
    obj.Orf.hmm_orf = {}

    fct2.get_hmm_genes(scaffold, table_hmm, gff_file)

    if resize:
        fct2.get_start_po(dico_seq)  # resize and calculate start position
    else:
        fct2.check_size(obj.TA_gene.genes_strand)  # eliminate te genes that have a length > threshold
    # for g in obj.TA_gene.genes_strand['-']:
        # print g
    # print "==" * 20
    # STEP : GENE PAIR ORGANISATION CHECKING
    fct2.get_adj()   # set the adj gene for each gene which has

    fct2.create_lonely_gene_list()
    initial_nb_lonely = len(obj.TA_gene.genes) - len(obj.TA_gene.linked)

    if obj.TA_gene.lonely is not None and rescue is True:

        # print 'THE LONELY GENE ARE :', len(obj.TA_gene.lonely['+'] + obj.TA_gene.lonely['-'])
        # for lone in obj.TA_gene.lonely['+'] + obj.TA_gene.lonely['-']:
            # print lone
        orf2.rescue_lonely_gene(dico_orf, gff_dico, scaffold, tmp_adj_orf_faa)

    # index to give a number of gene in the new gff file
    # Index start where the predicted gene end in order to not have same id for two gene in the round 2orf.
    fct2.score_TA_list(obj.TA_gene.genes_strand)
    # counter
    if info_contig_stat:
        fct2.contig_stat_manager(writer_stat, scaffold, initial_nb_lonely, rescue)

    # output
    if obj.TA_gene.linked:  # If there is some gene linked meaning if tere is TA system
        if output_human:
            pass
        if output_synthetic:
            floutS.write('==' * 5 + scaffold + '==' * 5 + '\n')
            fct2.write_result(obj.TA_gene.linked, floutS)
# Writing of stat information about the contig
if info_contig_stat:
    fl_stat.write("#Rescue lonely gene : {}\n".format(rescue))
    fl_stat.write("#Resize gene : {}\n".format(resize))
    writer_stat.writerow(obj.Gene.metaG_stat)

    fl_stat.close()
print using()

floutS.close()
