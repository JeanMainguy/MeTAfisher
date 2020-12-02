#!/usr/bin/env python3

"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 27 nov. 2020
License     : MIT
Maintainer  : jean.mainguy@outlook.fr


Program to retrieve toxin antitoxin (TA) systems in genomes or metagenomes.
"""

import Object_MetaF as obj
import Function_MetaF as fct
import Orf_MetaF as orf
import OutputFct_MetaF as out
import Score_MetaF as score
import csv
import argparse
import os
import logging
import sys


def init_logging(verbose_flag):
    """Initialise logging."""
    if verbose_flag:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')
        logging.info('command line: %s', ' '.join(sys.argv))

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")


def parse_arguments():
    """Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    """
    project_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    default_tadb_stat_dir = os.path.join(project_dir, "TADB_stat")
    default_hmm_db = os.path.join(default_tadb_stat_dir, 'ALL_plus_MET_curatted.hmm')

    parser = argparse.ArgumentParser(
        prog='MeTAfisher',
        description='Identification of Toxin Antitoxin Systems')

    parser.add_argument('--gff',
                        help="Path to the gff file.", required=True)
    parser.add_argument('--faa',
                        help="Path to the faa file. Fasta file of the annotated sequences proteins.",
                        required=True)

    parser.add_argument("-n", "--name",
                        help="Name of the Metagenome or Genome",
                        default='metafisher')

    parser.add_argument("-o", '--outdir',
                        help="Path to the result directory",
                        default='metafisher_results')

    parser.add_argument("--tadb_stat_dir",
                        help="Pathway of the tadb dtat dir. default is in the TADB_stat dir of the tool",
                        default=default_tadb_stat_dir)

    help_resize = "Resize the genes if they are too big for the thresholds by taking into account the possible starts along the sequence."
    help_resize += "To do only if the gene prediction is not trustable. This option requires the argument --fna."
    parser.add_argument('--resize', action='store_true',
                        help=help_resize)
    parser.add_argument('--fna',
                        help="Path to the fna file, which is a fasta file of annotated genes dna sequences. This file is required when the option --resize is given.",
                        required=False)

    rescue_help = "Applying the rescue step on lonely genes."
    rescue_help += "When a gene has TA domain but no other annotated gene around have one,"
    rescue_help += "or around are identified and searched for potential TA domain."
    rescue_help += "This option requires the argument --genomic_seq."
    parser.add_argument('--rescue', action='store_true',
                        help=rescue_help)

    parser.add_argument('--genomic_seq',
                        help="Path to the genomic sequence file, which is genomic sequence of the genome/metagenome. This file is required when the option --rescue is given.",
                        required=False)

    parser.add_argument('--contig_name', dest='contig_name',
                        help='Name of a specific contig to analysed. The program will analysed only this conitg')

    parser.add_argument("--hmm_db", default=default_hmm_db,
                        help="name of the HMM database")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()

    if args.rescue and args.genomic_seq is None:
        parser.error('argument --genomic_seq is required when rescue mode is on.')

    if args.resize and args.fna is None:
        parser.error('argument --fna is required when resize mode is on.')

    return args


def main():
    """Orchestrate the execution of the program.s"""
    args = parse_arguments()

    init_logging(args.verbose)

    metaG_name = args.name
    outdir = args.outdir
    # data_way = args.data_dir

    TADB_stat = args.tadb_stat_dir

    HMM_db = args.hmm_db

    # Name of a specific contig or False by default if False all contig will be analysed
    contig_name = args.contig_name

    # Flag
    resize = args.resize
    rescue = args.rescue

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    gff_file = args.gff
    faa_file = args.faa

    if resize:
        logging.info('Resize mode is on.')
        fna_file = args.fna

    if rescue:
        logging.info('Rescue mode is on.')
        genomic_seq_file = args.genomic_seq
        # prot sequence of adj orf will be written there
        tmp_adj_orf_faa = outdir + '/temporary_adjOrf.faa'

    table_hmm = outdir + '/output_HMM_table_' + metaG_name + '.txt'

    info_contig_stat = True
    output_human = True
    output_table = True
    output_gff = True
    dict_output = {'result_H': output_human, "result_S": output_human,
                   'result_T': output_table, 'result_GFF': output_gff}

    # Storing information as Gene class attribut to be use when we launch hmmsearch
    obj.Gene.outdir = outdir
    obj.Gene.hmmdb = HMM_db

    if not os.path.isfile(table_hmm):
        logging.warning('HMM launcher')
        print(table_hmm)
        table_hmm = fct.HMM_launcher(faa_file, metaG_name)
        print(table_hmm)

    # DISTANCE AND LENGTH DICO :
    # Dist and length of TA from TADB to mmake a proba
    file_len = TADB_stat + "/length_TA.csv"
    file_dist = TADB_stat + "/distance.csv"

    # domain vs domain : occurence of domain association in TADB inA pair
    file_domain_association = os.path.join(TADB_stat, "domain_domain_association.json")
    # domain type how often the domain is found in a toxin and in antitoxin
    file_domain_gene_type = os.path.join(TADB_stat, "/domain_gene_type.json")

    # CSV FILE DOMAINS
    csv_domain = os.path.join(TADB_stat, '/domaines_METAfisher.csv')

    with open(csv_domain, 'r') as csvdo:
        reader = csv.DictReader(csvdo)
        for row in reader:
            obj.Gene.domain_dict[row['hmm_name']] = {
                k: v for k, v in row.items() if k in ['acc', 'family', 'type']}
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

    # dictionnary gathering the 4 thresholds
    thresholds = {"lenMin": lenMin, "lenMax": lenMax,
                  "distanceMin": distanceMin, "distanceMax": distanceMax}

    """
    # OUTPUT
    #  dict_output is a dictionnary with the file handler of each kind of output
    # key : name of the output | value : fl or False if not wanted
    # extra key is "is_output" is True when there is at least one output required
    # dico output is update so no need to return it
    # writer stat is a csv file handler or it is False when the flag info_contig_stat is False as well
    total stat is a dictionnary with the sum of all the colonne. It is written at the end
    """
    writer_stat, total_stat, fl_stat = out.output_manager(
        outdir, metaG_name, thresholds, dict_output, info_contig_stat, rescue, resize)

    # Set the stops and starts codon as an atribut of the class Gene
    # Parameter for the orfinder

    id_genetCode = 11
    table = orf.getGeneticCode(id_genetCode)
    obj.Gene.codon_starts = table['start']
    obj.Gene.codon_stops = table['stop']

    # GFF dico use to delete the orf that correspond to a predicted gene
    # The program retrieve the positions of the end of the predicted gene
    # And don't proceess the orf that have the same end position

    if rescue:
        # Open step of the scaffold file ! used by fast_fasta function in orf file
        # the file won't be close until the end.Each time the algo won't have to
        # read the whle file again until the scaffold needed. We can do that
        # because the scaffold list has been sorted
        # New version use this dictionnary:
        dico_orf = {}
        dico_orf['fl'] = open(genomic_seq_file, 'r')
        dico_orf['line'] = next(dico_orf['fl'])

        gff_dico = {}
        fl_csv = open(gff_file, 'r')
        gff_dico['csv'] = csv.reader(fl_csv, delimiter='\t')
        gff_dico['line'] = next(gff_dico['csv'])

    if resize:
        # Open file fna to retrieve sequence of the predicted gene
        # variable are stored in a dico to not have to retourned it every time !!
        dico_seq = {}
        dico_seq["fl"] = open(fna_file, 'r')

        dico_seq["line"] = dico_seq["fl"].readline()
        # PUT info in a dico to not have to retourned it every time !! dico is used in check_size
        dico_seq["codon_start"] = table['start']
        # print "start", table['start']

    # SCORE PREPARaTION
    k = 20
    # dist_mltp = 7
    # len_mltp = 7
    bonus_start = 7  # +7 is given to configuration that start with their initial start and not a start determined by the program

    obj.Gene.length_proba = score.score_manager(int(lenMin / 3), int(lenMax / 3), file_len, k)
    obj.Gene.distance_proba = score.score_manager(distanceMin, distanceMax, file_dist, k)
    obj.Gene.dict_domain_association = score.decoder(file_domain_association)
    obj.Gene.dict_domain_gene_type = score.decoder(file_domain_gene_type)

    # Take every scaffold present in the HMM output and sort them in order to
    # be able to retrieve their sequence correctly in the input file
    scaffold_list = fct.get_list_scaffold(table_hmm)
    if not contig_name:  # if the contig name option to give only one contig is not provided then contig_name is False and all contig are analysed
        scaffold_list = sorted(fct.get_list_scaffold(table_hmm))
    else:  # contig name is only is the name of one contig
        if contig_name not in scaffold_list:
            raise Exception('The contig name is incorrect or no hit have been found by hmmsearch')
        scaffold_list = [contig_name]

    logging.info(f'number of sequence to analysed {len(scaffold_list)}')

    # Loop : each scaffold is treated independntly here
    for scaffold in scaffold_list:

        logging.info(f"Analysing {scaffold}")

        # Reset obj.TA_gene and Orf class attribut
        obj.TA_gene.genes = []
        obj.TA_gene.genes_strand = {'+': [], '-': []}
        obj.TA_gene.linked = set()
        obj.TA_gene.lonely = None
        obj.Orf.adj_orf = {}
        obj.Orf.adj_orf_index = 0
        obj.Orf.hmm_orf = {}

        fct.get_hmm_genes(scaffold, table_hmm, gff_file)

        if resize:
            fct.get_start_po(dico_seq)  # resize and calculate start position
        else:
            # eliminate te genes that have a length > threshold
            fct.check_size(obj.TA_gene.genes_strand)

        # STEP : GENE PAIR ORGANISATION CHECKING
        fct.get_adj()   # set the adj gene for each gene which has

        fct.create_lonely_gene_list()

        initial_nb_lonely = len(obj.TA_gene.genes) - len(obj.TA_gene.linked)

        if obj.TA_gene.lonely is not None and rescue:
            orf.rescue_lonely_gene(dico_orf, gff_dico, scaffold, tmp_adj_orf_faa)

        # score.score_TA_list(obj.TA_gene.genes_strand, bonus_start)

        score.score_TA_list(obj.TA_gene.linked, bonus_start)

        # Write stat
        if info_contig_stat:
            out.contig_stat_manager(writer_stat, scaffold, initial_nb_lonely, rescue, total_stat)

        # write output
        # If there is some gene linked meaning if tere is TA system
        if obj.TA_gene.linked and dict_output['is_output']:
            out.write_result(obj.TA_gene.linked, dict_output, scaffold)

    # Total Stat information about the Metagenome
    if info_contig_stat:
        writer_stat.writerow(total_stat)
        fl_stat.close()

    for kfl in dict_output:
        try:
            dict_output[kfl].close()
        except AttributeError:
            pass


if __name__ == '__main__':
    main()
