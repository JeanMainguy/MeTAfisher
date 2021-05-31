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
import gzip


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
    default_hmm_db = os.path.join(default_tadb_stat_dir, 'TA_domains.hmm')

    parser = argparse.ArgumentParser(
        prog='MeTAfisher',
        description='Identification of Toxin Antitoxin Systems', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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
                        help="Pathway of the tadb stat dir",
                        default=default_tadb_stat_dir)

    # help_resize = "Resize the genes if they are too big for the thresholds by taking into account the possible starts along the sequence."
    # help_resize += "To do only if the gene prediction is not trustable. This option requires the argument --fna."
    # parser.add_argument('--resize', action='store_true',
    #                     help=help_resize)
    # parser.add_argument('--fna',
    #                     help="Path to the fna file, which is a fasta file of annotated genes dna sequences. This file is required when the option --resize is given.",
    #                     required=False)

    rescue_help = "Applying the rescue step on lonely genes."
    rescue_help += "When a gene has a TA domain but no other annotated gene around have one,"
    rescue_help += "Adjacent ORF around are identified and searched for potential TA domain."
    rescue_help += "This option requires the argument --genomic_seq."
    parser.add_argument('--rescue', action='store_true',
                        help=rescue_help)

    parser.add_argument('--genomic_seq',
                        help="Path to the genomic sequence file, which is genomic sequence of the genome/metagenome. This file is required when the option --rescue is given.",
                        required=False)

    parser.add_argument('--contig_name', dest='contig_name',
                        help='Name of a specific contig to analysed. The program will analysed only this conitg')

    parser.add_argument("--hmm_db", default=default_hmm_db,
                        help="path of the HMM database.")

    parser.add_argument("--diamond_db",
                        help="path of the diamond database.")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()

    if args.rescue and args.genomic_seq is None:
        parser.error('argument --genomic_seq is required when rescue mode is on.')

    # if args.resize and args.fna is None:
    #     parser.error('argument --fna is required when resize mode is on.')

    return args


def main():
    """Orchestrate the execution of the program"""
    args = parse_arguments()

    init_logging(args.verbose)

    metaG_name = args.name
    outdir = args.outdir
    # data_way = args.data_dir

    tadb_stat_dir = args.tadb_stat_dir

    hmm_db = args.hmm_db
    diamond_db = args.diamond_db

    # Name of a specific contig or False by default if False all contig will be analysed
    contig_name = args.contig_name

    # Flag
    # resize = args.resize
    rescue = args.rescue

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    gff_file = args.gff
    faa_file = args.faa

    # if resize:
    #     logging.info('Resize mode is on.')
    #     fna_file = args.fna

    if rescue:
        logging.info('Rescue mode is on.')
        genomic_seq_file = args.genomic_seq

    info_contig_stat = True
    output_human = True
    output_gene_tsv = True
    output_pair_tsv = True
    output_gff = True
    dict_output = {'result_H': output_human, "result_S": output_human,
                   'result_TA_genes': output_gene_tsv, 'result_TA_pairs': output_pair_tsv, 'result_GFF': output_gff}

    # DISTANCE AND LENGTH DICT :
    # Dist and length of TA from TADB to mmake a proba
    file_len = os.path.join(tadb_stat_dir, "length_TA.csv")
    file_dist = os.path.join(tadb_stat_dir, "distance.csv")

    # domain vs domain : occurence of domain association in TADB inA pair
    file_domain_association = os.path.join(tadb_stat_dir, "domain_domain_association.json")
    # domain type how often the domain is found in a toxin and in antitoxin
    file_domain_gene_type = os.path.join(tadb_stat_dir, "domain_gene_type.json")

    # CSV FILE DOMAINS
    csv_domain = os.path.join(tadb_stat_dir, 'TA_domains_info.tsv')
    info_domains = {}
    with open(csv_domain, 'r') as csvdo:
        reader = csv.DictReader(csvdo, delimiter='\t')
        for row in reader:
            info_domains[row['hmm_name']] = {
                k: v for k, v in row.items() if k in ['acc', 'family', 'type']}

    # THRESHOLD :
    # There are the first the thrshold there very large then the proabilty score step is going to define the better conf
    # treshold size of gene
    lenMin = 10 * 3  # VERY IMPORTANT LENGTH HAVE TO BE GIVEN IN AA and then transform in nt
    lenMax = 500 * 3

    obj.Gene.length_min = lenMin
    obj.Gene.length_max = lenMax

    # threshold distance for tandem in nucleotides
    distanceMin = -100
    distanceMax = 300

    obj.Gene.distanceMin = distanceMin
    obj.Gene.distanceMax = distanceMax

    # dictionnary gathering the 4 thresholds
    thresholds = {"lenMin": lenMin, "lenMax": lenMax,
                  "distanceMin": distanceMin, "distanceMax": distanceMax}

    writer_stat, total_stat, fl_stat = out.output_manager(
        outdir, metaG_name, thresholds, dict_output, info_contig_stat, rescue)

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
        orf_dict = {}
        proper_open = gzip.open if genomic_seq_file.endswith('.gz') else open
        orf_dict['fl'] = proper_open(genomic_seq_file, 'rt')
        orf_dict['line'] = next(orf_dict['fl'])

        gff_dict = {}
        proper_open = gzip.open if gff_file.endswith('.gz') else open
        fl_csv = proper_open(gff_file, 'rt')
        gff_dict['csv'] = csv.reader(fl_csv, delimiter='\t')
        gff_dict['line'] = next(gff_dict['csv'])

    # if resize:
    #     # Open file fna to retrieve sequence of the predicted gene
    #     # variable are stored in a dico to not have to retourned it every time !!
    #     fna_seq_dict = {}
    #     fna_seq_dict["fl"] = open(fna_file, 'r')
    #
    #     fna_seq_dict["line"] = fna_seq_dict["fl"].readline()
    #     # PUT info in a dico to not have to retourned it every time !! dico is used in check_size
    #     fna_seq_dict["codon_start"] = table['start']

    # Score preparation
    k = 20
    length_proba = score.score_manager(int(lenMin / 3), int(lenMax / 3), file_len, k)
    distance_proba = score.score_manager(distanceMin, distanceMax, file_dist, k)
    dict_domain_association = score.decoder(file_domain_association)
    dict_domain_gene_type = score.decoder(file_domain_gene_type)

    score_dict = {"length_proba": length_proba,
                  "distance_proba": distance_proba,
                  "domain_association": dict_domain_association,
                  "domain_gene_type": dict_domain_gene_type}

    element_to_rm = -2 if faa_file.endswith('.gz') else -1
    simple_faa_name = '.'.join(os.path.basename(faa_file).split('.')[:-element_to_rm])
    hmmsearch_result_file = os.path.join(outdir, f"{simple_faa_name}.hmmsearch")
    fct.hmmsearch(faa_file, hmm_db, hmmsearch_result_file)

    gene_to_hits = fct.get_ta_genes_from_hmmsearch(hmmsearch_result_file)

    if diamond_db:
        diamond_result_file = os.path.join(outdir, f"{simple_faa_name}.diamond")
        fct.diamond_blastp(faa_file, diamond_db, diamond_result_file)
        gene_to_hits = fct.get_ta_genes_from_diamond(
            diamond_result_file, gene_to_hits, min_coverage=95, min_pident=95)

    fct.annotate_ta_hits(gene_to_hits.values(), info_domains, dict_domain_gene_type)

    contig_to_genes = fct.get_genes_by_contigs(gene_to_hits, gff_file)

    for contig, genes in contig_to_genes:
        if contig_name and contig != contig_name:
            continue

        #logging.info(f'{contig}: {len(genes)} TA genes.')

        fct.check_size(genes)

        fct.compute_gene_adjacency(genes)

        initial_nb_lonely = len(genes) - len(list(fct.get_linked_genes(genes)))

        if initial_nb_lonely and rescue:
            adj_orfs = orf.get_adjacent_orfs(orf_dict, gff_dict, contig, genes)
            ta_orfs = orf.identify_ta_orfs(adj_orfs, genes, outdir,
                                           hmm_db, info_domains, dict_domain_gene_type)
        else:
            ta_orfs = []
            adj_orfs = []

        score.score_TA_list(genes, score_dict)
        # Write stat
        if info_contig_stat:
            out.contig_stat_manager(writer_stat, contig, initial_nb_lonely,
                                    rescue, total_stat, genes, adj_orfs)

        # write output
        if dict_output['is_output']:
            #out.write_result(fct.get_linked_genes(genes), dict_output, contig)
            out.write_result(genes, dict_output, contig)

    # Total Stat information about the Metagenome
    if info_contig_stat:
        logging.info(
            f'MeTAfisher has identified {total_stat["linked gene"]} genes belonging to putative TA systems.')
        writer_stat.writerow(total_stat)
        fl_stat.close()

    logging.info(f"Results have been written in {args.outdir}")

    for kfl in dict_output:
        try:
            dict_output[kfl].close()
        except AttributeError:
            pass


if __name__ == '__main__':
    main()
