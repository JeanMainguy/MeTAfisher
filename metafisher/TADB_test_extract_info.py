import argparse
from Bio import SeqIO
import re
import logging
from os import path, makedirs, listdir
import prepare_data_from_accession as gb_download
import prepare_data_from_gbfile as gb_parser
from Bio import Entrez
import json


def build_genome_TA_info_dict(prot_acc_list, gb_prot_dir, info_dict):

    gb_download.download_many_gb_files(
        prot_acc_list, gb_prot_dir, db="Protein", rettype="gb", create_folder=False, batch_size=50)

    for acc, prot_acc_info in prot_acc_list.items():
        filename = path.join(gb_prot_dir, acc+".gb")
        TADB_nb = prot_acc_info["TADB_nb"]
        # gb_download.download_gb_file(acc, filename, db="Protein", rettype="gb")
        genome_accession = parse_prot_gb(filename, prot_acc_info)
        if genome_accession:
            # info_dict.setdefault(genome_accession, {}).setdefault(TADB_nb, []).append(prot_acc_info)
            info_dict.setdefault(genome_accession, []).append(prot_acc_info)


def parse_prot_gb(gb_file, prot_acc_info):
    # get genome accesion found in the header of the gb file
    # Example: DBSOURCE    REFSEQ: accession NC_000854.2

    with open(gb_file, "r") as input_handle:
        for i, record in enumerate(SeqIO.parse(input_handle, "genbank")):
            if 'db_source' in record.annotations:
                accession = record.annotations['db_source'].split(' ')[-1]

            strand = '+'
            for f in record.features:
                if f.type == "CDS" and "coded_by" in f.qualifiers:
                    info_genome = f.qualifiers["coded_by"].pop()

                    if info_genome.startswith("complement(") and info_genome.endswith(")"):
                        strand = '-'
                        info_genome = info_genome.replace('complement(', '')[:-1]
                    try:
                        genome_acc, positions = info_genome.split(':')
                    except ValueError:
                        logging.warning(
                            'Problem while parsing genome_accession in protein gb {}'.format(info_genome))
                        return False
                    prot_acc_info['strand'] = strand
                    prot_acc_info['positions'] = positions
                    prot_acc_info['start'], prot_acc_info['end'] = positions.split("..")
                    return genome_acc


def extract_protein_acc(file):
    prot_accessions = {}
    pattern = re.compile("TADB\|([AT]{1,2})([\d]*).*ref\|([^\|]*)\|([^[]*)\[([^]]*)\]")
    seqs = SeqIO.parse(file, "fasta")
    for seq in seqs:
        result = pattern.match(seq.description)
        if result:
            type = result.group(1)
            TADB_nb = result.group(2)
            acc = result.group(3)
            name = result.group(4)
            specie = result.group(5)

            prot_accessions[acc] = {'acc': acc, "name": name,
                                    'specie': specie, 'TADB_nb': TADB_nb, 'type': type}
        else:
            logging.warning('regex failed for seq {}'.format(str(seq).replace('\n', ' | ')))

    return prot_accessions


def check_consistency(info_dict, gb_genome_dir):
    genomes_acc = [version.split('.')[0] for version in info_dict]
    redundant_genomes = {genome for genome in genomes_acc if genomes_acc.count(genome) > 1}

    print(redundant_genomes)


def encode_TAT_info(info_dict, output_file):
    with open(output_file, 'w') as file:
        json.dump(info_dict, file, indent=4, sort_keys=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='Extract TA info from fasta files of toxin and antitoxin', description='')

    parser.add_argument("toxin_file", help="fasta file of toxin")
    parser.add_argument("antitoxin_file", help="fasta file of antitoxin")

    parser.add_argument(
        "email", help="mail address to be identified by the ncbi when using the module Entrez")
    parser.add_argument("genome_dir", help="directory where genomes are going to be store")
    parser.add_argument("output_result_file", help="json file where results are stored")
    args = parser.parse_args()

    toxin_fl = args.toxin_file
    antitoxin_fl = args.antitoxin_file
    Entrez.email = args.email

    gb_prot_dir = path.join(path.dirname(toxin_fl), 'gb_proteins')
    gb_genome_dir = args.genome_dir
    output_file = args.output_result_file

    print(gb_prot_dir, gb_genome_dir, output_file)
    try:
        makedirs(gb_prot_dir)
    except OSError:
        pass

    try:
        makedirs(gb_genome_dir)
    except OSError:
        pass

    print(gb_prot_dir)
    print(gb_genome_dir)

    ATs = extract_protein_acc(antitoxin_fl)
    Ts = extract_protein_acc(toxin_fl)

    print('AT:', len(ATs))
    print('T:', len(Ts))
    info_dict = {}

    build_genome_TA_info_dict(ATs, gb_prot_dir, info_dict)
    build_genome_TA_info_dict(Ts, gb_prot_dir, info_dict)

    print("genomes:", len(info_dict))

    gb_download.download_many_gb_files(
        list(info_dict.keys()), gb_genome_dir, db="nucleotide", rettype="gbwithparts")

    for gb_dir in listdir(gb_genome_dir):
        if gb_dir in info_dict:
            data_path_base = path.join(gb_genome_dir, gb_dir, gb_dir)
            filename = data_path_base + '.gb'

            if len(listdir(path.join(gb_genome_dir, gb_dir))) < 2:
                gb_parser.from_gb_to_required_format(data_path_base, filename)
            else:
                # print('formated files seems to exist already')
                logging.info('formated files seems to exist already')

    check_consistency(info_dict, gb_genome_dir)

    encode_TAT_info(info_dict, output_file)
