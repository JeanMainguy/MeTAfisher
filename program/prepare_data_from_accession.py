from Bio import Entrez
import argparse
from os import path, makedirs
import prepare_data_from_gbfile as gb_parser


def isAccValid(acc):
    try:
        record = Entrez.read(Entrez.esummary(db="nucleotide", id=acc))
    except RuntimeError:
        return False
    except Exception as e:
        print "Unknown error:", e
        return e


def download_gb_file(acc, filename):
    if not path.isfile(filename):
        # Downloading...
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gbwithparts", retmode="text")
        except Exception as e:
            if isAccValid(acc):
                print "Unknown error:", e
                raise e
            else:
                raise ValueError('Accession is not correct')

        if not path.exists(path.dirname(filename)):
            makedirs(path.dirname(filename))

        with open(filename, "w") as fl:
            fl.write(handle.read())
        handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='prepare_data_from_accession', description='Prepare data to be analysed by MeTAfisher starting with the accession id of the genome')

    parser.add_argument("accession", help="accession of the genome. example:NC_009467")
    parser.add_argument(
        "email", help="mail address to be identified by the ncbi when using the module Entrez")
    args = parser.parse_args()
    accession = args.accession

    Entrez.email = args.email

    data_path = 'data/'

    gb_file = path.join(data_path, accession, accession + '.gb')

    download_gb_file(accession, gb_file)

    data_path_base = path.abspath(gb_file).split(".gb")[0]
    gb_parser.from_gb_to_required_format(data_path_base, gb_file)
