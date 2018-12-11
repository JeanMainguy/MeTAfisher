import argparse
from Bio import SeqIO
import re
import logging
from os import path, makedirs, listdir
import prepare_data_from_accession as gb_download
import prepare_data_from_gbfile as gb_parser
from Bio import Entrez


def build_genome_TA_info_dict(prot_acc_list, type, gb_prot_dir, info_dict):
    for acc in prot_acc_list:
        filename = path.join(gb_prot_dir, acc+".gb")
        gb_download.download_gb_file(acc, filename, db="Protein", rettype="gb")
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='Extract TA info from fasta files of toxin and antitoxin', description='')

    parser.add_argument("toxin_file", help="fasta file of toxin")
    parser.add_argument("antitoxin_file", help="fasta file of antitoxin")

    parser.add_argument(
        "email", help="mail address to be identified by the ncbi when using the module Entrez")

    args = parser.parse_args()

    toxin_fl = args.toxin_file
    antitoxin_fl = args.antitoxin_file
    Entrez.email = args.email

    gb_prot_dir = path.join(path.dirname(toxin_fl), 'gb_proteins_single')
    gb_genome_dir = path.join(path.dirname(toxin_fl), 'gb_genomes_single')
    print(gb_prot_dir)
    print(gb_genome_dir)

    ATs = extract_protein_acc(antitoxin_fl)
    Ts = extract_protein_acc(toxin_fl)

    print(len(ATs))
    print(len(Ts))
    info_dict = {}

    build_genome_TA_info_dict(ATs, "antitoxin", gb_prot_dir, info_dict)
    build_genome_TA_info_dict(Ts, "toxin", gb_prot_dir, info_dict)

    print("genomes:", len(info_dict))
    # print(info_dict)
    try:
        makedirs(gb_genome_dir)
    except OSError:
        pass

    gb_download.download_many_gb_files(
        info_dict.keys(), gb_genome_dir, db="nucleotide", rettype="gbwithparts")

    for gb_dir in listdir(gb_genome_dir):
        if gb_dir in info_dict:
            data_path_base = path.join(gb_genome_dir, gb_dir, gb_dir)
            filename = data_path_base + '.gb'
            print(data_path_base)
            print(filename)
            gb_parser.from_gb_to_required_format(data_path_base, filename)

    #
    # for genome in info_dict:
    #
    #     print(genome)
    #     gb_download.download_gb_file(genome, filename)
    #     data_path_base = path.abspath(filename).split(".gb")[0]
    #     gb_parser.from_gb_to_required_format(data_path_base, filename)
    # gb_download.download_many_gb_files(ATs, gb_prot_dir, db='Protein', rettype='gb')
