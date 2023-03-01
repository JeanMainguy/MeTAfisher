#!/usr/bin/env python3

"""
Module      : Main
Description : Build stat file required by metafisher.
Copyright   : (c) Jean Mainguy, 18 jan. 2021
License     : MIT
Maintainer  : jean.mainguy@outlook.fr


Program to build stat files of TADB to be used in metafisher.
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
import re
import json


def encoder(filename, dict):
    logging.info(f'writing json file {filename}')
    with open(filename, 'w') as file:
        json.dump(dict, file, indent=4)


def write_table(filename, domain_domain_dict):
    logging.info(f'writing tsv file {filename}')
    domain_list = list(domain_domain_dict)
    with open(filename, 'w') as fl:
        fl.write('domains\t'+'\t'.join(domain_list)+'\n')
        for d in domain_list:
            line = [d]
            for d_n in domain_list:
                nb_pairs_in_common = '0' if d_n not in domain_domain_dict[d] else str(
                    domain_domain_dict[d][d_n])
                line.append(nb_pairs_in_common)
            fl.write('\t'.join(line)+'\n')


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
        description='Identification of Toxin Antitoxin Systems',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--toxin_faa',
                        help="Path to the reference toxin protein sequence.", required=True)

    parser.add_argument('--antitoxin_faa',
                        help="Path to the reference antitoxin protein sequence.", required=True)

    parser.add_argument("-o", '--outdir',
                        help="Path to the directory where files will be written",
                        default=default_tadb_stat_dir)

    parser.add_argument("--hmm_db", default=default_hmm_db,
                        help="path of the HMM database.")

    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()

    return args


def domain_domain_pair_association(domain_type_dict, opposite_type_dict={'T': 'AT', 'AT': 'T'}):
    """
    Compute domain domain association.

    domain_type_dict is a {domain_name:{T:[gene_ids], AT:[gene_ids]} ... }

    """

    domain_domain_dict = {}
    for domain, type2genes in domain_type_dict.items():

        domain_dict = domain_domain_dict.setdefault(domain, {})

        for domain_next, type2genes_next in domain_type_dict.items():
            if domain_next in domain_dict:
                continue
            domain_dict_next = domain_domain_dict.setdefault(domain_next, {})

            pairs = []

            for type, opposite_type in opposite_type_dict.items():
                genes = type2genes.setdefault(type, [])
                genes_next = type2genes_next.setdefault(opposite_type, [])

                pairs += list(set(genes_next) & set(genes))

            if len(pairs) > 0:
                domain_dict[domain_next] = pairs
                domain_dict_next[domain] = pairs

    return domain_domain_dict


def extract_db_and_gene_number(gene_id):
    # gene_number_pattern = re.compile(r'[^\d]+(\d+)')
    gene_number_pattern = re.compile(r'\|(AT|T)(.*+)')

    try:
        gene_number = gene_number_pattern.match(gene_id).group(1)
    except AttributeError:
        raise AttributeError(f'regex pattern failed to extract gene number in id {gene_id}.')

    db_name = gene_id.split('|')[0]
    return db_name, gene_number


def parse_tadb_ids(seq_file):
    gene_number_to_id = {}
    with open(seq_file) as fl:

        gene_ids = (l.split()[0][1:] for l in fl if l.startswith('>'))

        for i, id in enumerate(gene_ids):
            
            
            db_and_nb = extract_db_and_gene_number(id)
            # db_and_nb = tuple(id.split('|'))
            if db_and_nb in gene_number_to_id:
                #raise ValueError(f'db and gene number {db_and_nb} are used twice to identify a sequence in {seq_file}')
                logging.critical(f'Gene id {id} is used more than once in {seq_file}')

            gene_number_to_id[db_and_nb] = id

    # assert i == len(gene_number_to_id), "Same gene number is used at least twice. Check consistency of sequence id of {seq_file}"
    return gene_number_to_id


def get_genes_association(toxin_file, antitoxin_file):
    toxin_id_to_number = parse_tadb_ids(toxin_file)
    antitoxin_id_to_number = parse_tadb_ids(antitoxin_file)

    if len(toxin_id_to_number) != len(antitoxin_id_to_number):
        logging.critical(
            f"Not the same number of toxin genes ({len(toxin_id_to_number)}) and antitoxin genes ({len(antitoxin_id_to_number)}). check files: {toxin_file} and {antitoxin_file}")

    gene_numbers = set(toxin_id_to_number) | set(antitoxin_id_to_number)
    genes_association = {}
    genes_type = {}

    for db_name, gene_number in gene_numbers:

        try:
            antitoxin_id = antitoxin_id_to_number[(db_name, gene_number)]
            genes_type[antitoxin_id] = {'AT': 1, 'T': 0}

        except KeyError:

            logging.critical(f'No antitoxin gene with id {db_name}|AT{gene_number}')
            antitoxin_id = None
        try:
            toxin_id = toxin_id_to_number[(db_name, gene_number)]
            genes_type[toxin_id] = {'AT': 0, 'T': 1}
        except KeyError:
            logging.critical(f'No toxin with id {db_name}|T{gene_number}')
            toxin_id = None

        if toxin_id and antitoxin_id:
            genes_association[antitoxin_id] = {toxin_id: 1}
            genes_association[toxin_id] = {antitoxin_id: 1}

    return genes_association, genes_type


def domains_genes_association(hmm_result, domain_type_dict):

    gene_id_parser = re.compile(r"(?P<db_name>[^|]+)\|(?P<type>[A,T]{1,2})(?P<gene_number>\d+)")

    type_name = set()

    for hmmhit in fct.hmm_result_parser(hmm_result):

        domain = hmmhit.query_name
        re_result = gene_id_parser.match(hmmhit.target_name)
        type = re_result.group("type")  # T or AT
        gene_number = re_result.group("gene_number")

        domain_type_dict.setdefault(domain, {}).setdefault(type, []).append(gene_number)
        # domain_gene_dict.setdefault(domain, []).append(gene_number)
        type_name.add(type)

    assert len(type_name) == 1

    return type_name.pop()


def run_hmmsearch(faa_file, hmm_db, outdir):
    element_to_rm = -2 if faa_file.endswith('.gz') else -1
    simple_name = '.'.join(os.path.basename(faa_file).split('.')[:-element_to_rm])
    hmm_result = os.path.join(outdir, f"{simple_name}.hmmsearch")

    fct.hmmsearch(faa_file, hmm_db, hmm_result)

    return hmm_result


def main():
    """Orchestrate the execution of the program"""
    args = parse_arguments()

    init_logging(args.verbose)

    toxin_seq_file = args.toxin_faa
    antitoxin_seq_file = args.antitoxin_faa

    outdir = args.outdir
    hmm_db = args.hmm_db

    domain_type_dict = {}

    # get_genes_association(toxin_seq_file, 'T', domain_type_dict)
    genes_association, genes_type = get_genes_association(toxin_seq_file, antitoxin_seq_file)

    hmm_result_T = run_hmmsearch(toxin_seq_file, hmm_db, outdir)
    type_T = domains_genes_association(hmm_result_T, domain_type_dict)

    hmm_result_AT = run_hmmsearch(antitoxin_seq_file, hmm_db, outdir)
    type_AT = domains_genes_association(hmm_result_AT, domain_type_dict)

    types = {type_T: type_AT, type_AT: type_T}

    domain_domain_dict = domain_domain_pair_association(domain_type_dict, types)

    # count
    domain_domain_count_asso = {}
    for d, d_n in domain_domain_dict.items():
        domain_domain_count_asso[d] = {k: len(v) for k, v in d_n.items()}

    table_file = os.path.join(outdir, 'domain_domain_association.tsv')
    write_table(table_file, domain_domain_count_asso)

    domain_domain_count_asso.update(genes_association)

    json_d_d_file = os.path.join(outdir, 'domain_domain_association.json')
    encoder(json_d_d_file, domain_domain_count_asso)

    for domain, type_dict in domain_type_dict.items():
        domain_type_dict[domain] = {type_gene: len(set(gene_ids))
                                    for type_gene, gene_ids in type_dict.items()}

    domain_type_dict.update(genes_type)

    json_d_t_file = os.path.join(outdir, 'domain_gene_type.json')
    encoder(json_d_t_file, domain_type_dict)


if __name__ == '__main__':
    main()
