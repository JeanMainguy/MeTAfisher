import argparse
from Bio import SeqIO
import re
import logging
from os import path, makedirs, listdir
import prepare_data_from_accession as gb_download
import prepare_data_from_gbfile as gb_parser
from Bio import Entrez
import json


def build_genome_TA_info_dict(prot_acc_list, type, gb_prot_dir, info_dict):

    gb_download.download_many_gb_files(
        prot_acc_list, gb_prot_dir, db="Protein", rettype="gb", create_folder=False, batch_size=50)

    for i, acc in enumerate(prot_acc_list):
        # if round(i/len(prot_acc_list)*100) > last + 5:
        #     print("{}% processed".format(round(i/len(prot_acc_list))))
        #     last = round(i/len(prot_acc_list)*100)
        filename = path.join(gb_prot_dir, acc+".gb")
        # gb_download.download_gb_file(acc, filename, db="Protein", rettype="gb")

        genome_accession = get_genome_accession_from_prot_gb(filename)
        if genome_accession:
            info_dict.setdefault(genome_accession, {}).setdefault(type, []).append(acc)


def get_genome_accession_from_prot_gb(gb_file):
    # get genome accesion found in the header of the gb file
    # Example: DBSOURCE    REFSEQ: accession NC_000854.2

    with open(gb_file, "r") as input_handle:
        for i, record in enumerate(SeqIO.parse(input_handle, "genbank")):
            try:
                accession = record.annotations['db_source'].split(' ')[-1]
                return accession
            except KeyError:
                logging.warning(
                    'DBSOURCE line is missing in gb_file {} record {}'.format(gb_file, i+1))


def extract_protein_acc(file):
    prot_accessions = []
    pattern = re.compile(".*ref\|([^\|]*)")
    seqs = SeqIO.parse(file, "fasta")
    for seq in seqs:
        result = pattern.match(seq.description)
        if result:
            acc = result.group(1)
            prot_accessions.append(acc)
        else:
            logging.warning('regex failed for seq {}'.format(str(seq).replace('\n', ' | ')))

    return prot_accessions


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
    args = parser.parse_args()

    toxin_fl = args.toxin_file
    antitoxin_fl = args.antitoxin_file
    Entrez.email = args.email

    gb_prot_dir = path.join(path.dirname(toxin_fl), 'gb_proteins')
    gb_genome_dir = args.genome_dir
    output_file = path.join(args.genome_dir, 'TADB_TAT_info.json')

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

    build_genome_TA_info_dict(ATs, "antitoxin", gb_prot_dir, info_dict)
    build_genome_TA_info_dict(Ts, "toxin", gb_prot_dir, info_dict)

    encode_TAT_info(info_dict, output_file)

    print("genomes:", len(info_dict))

    gb_download.download_many_gb_files(
        info_dict.keys(), gb_genome_dir, db="nucleotide", rettype="gbwithparts")

    for gb_dir in listdir(gb_genome_dir):
        if gb_dir in info_dict:
            data_path_base = path.join(gb_genome_dir, gb_dir, gb_dir)
            filename = data_path_base + '.gb'
            print(data_path_base)
            print(filename)
            print("len(listdir(path.join(gb_genome_dir, gb_dir)))",
                  len(listdir(path.join(gb_genome_dir, gb_dir))))
            if len(listdir(path.join(gb_genome_dir, gb_dir))) < 2:
                gb_parser.from_gb_to_required_format(data_path_base, filename)
            else:
                print('formated files seems to exist already')
                logging.info('formated files seems to exist already')