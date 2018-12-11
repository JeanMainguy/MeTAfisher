from Bio import Entrez, SeqIO
import argparse
from os import path, makedirs, listdir
import prepare_data_from_gbfile as gb_parser
import logging
import time

try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


def download_many_gb_files(acc_list, output, db='Protein', rettype='gb'):
    files = {f.replace(".gb", '') for f in listdir(output)}  # remove ext
    acc_list_to_download = set(acc_list) - {f for f in files if f in acc_list}
    print(files)
    print(acc_list)
    print(acc_list_to_download)
    if not acc_list_to_download:
        print("Nothing to download")
        return
    search_results = Entrez.read(Entrez.epost(db=db, id=",".join(acc_list_to_download)))
    count = len(acc_list_to_download)
    downloaded_acc = []
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    batch_size = 3
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                fetch_handle = Entrez.efetch(db=db,
                                             rettype=rettype, retmode="text",
                                             retstart=start, retmax=batch_size,
                                             webenv=webenv, query_key=query_key,
                                             idtype="acc")
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(15)
                else:
                    raise

        for record in SeqIO.parse(fetch_handle, "genbank"):
            acc = "{}.{}".format(record.annotations['accessions']
                                 [0], record.annotations['sequence_version'])

            print(acc)
            downloaded_acc.append(acc)
            filename = path.join(output, acc, "{}.gb".format(acc))
            makedirs(path.dirname(filename))

            SeqIO.write(record, filename, "gb")

        # with open(filename, "w") as fl:
        #     fl.write(data)
        #     fetch_handle.close()
        print('downloaded_acc', downloaded_acc)
        print('acc_list_to_download', acc_list_to_download)
        if downloaded_acc == list(acc_list_to_download):
            print('Every acc have been downloaded')
        else:
            missed_acc = acc_list_to_download - set(downloaded_acc)
            logging.warning('{} accessions have not been downloaded: {}'.format(
                len(missed_acc), missed_acc))


def isAccValid(acc, db):
    try:
        record = Entrez.read(Entrez.esummary(db=db, id=acc))
    except RuntimeError:
        return False
    except Exception as e:
        print("Unknown error:", e)
        return e


def download_gb_file(acc, filename, db="nucleotide", rettype="gbwithparts"):
    if not path.isfile(filename):
        # Downloading...
        try:
            handle = Entrez.efetch(db=db, id=acc, rettype=rettype, retmode="text")
        except Exception as e:
            if isAccValid(acc, db):
                print("Unknown error:", e)
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
