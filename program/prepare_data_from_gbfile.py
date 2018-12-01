import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from os import path
import argparse


def giveQualifiersInfo(qualifiers, keys):
    liste = []
    for key in keys:
        if key in qualifiers:
            liste.append("{}={}".format(key, qualifiers[key][0]))
    return liste


def write_gff(feat, fl_gff, chrm, attributes, nbgene):
    strand = '+' if feat.strand == 1 else '-'
    attributes = "{}|{};".format(chrm, str(nbgene)) + ';'.join(attributes)
    # attributes += giveQualifiersInfo(feat.qualifiers, ["locus_tag", 'protein_id', 'product'])

    list_gff = [chrm, 'Genbank', feat.type, str(
        feat.location.start.real + 1), str(feat.location.end.real), ".", strand, ".", attributes+'\n']
    fl_gff.write('\t'.join(list_gff))


def write_faa(feat, faa_fl, chrm, seq, attributes, nbgene):
    trans_table = int(feat.qualifiers['transl_table'][0])
    # if translation is available we use it prioritarly
    if 'translation' in feat.qualifiers:
        seqAA = Seq(feat.qualifiers['translation'][0])
    else:
        try:
            seqAA = seq.translate(table=trans_table, to_stop=True, cds=True)
        except:
            return False

    record = SeqRecord(seqAA, id="{}|{} ".format(chrm, str(nbgene)),
                       name="", description=' '.join(attributes))
    SeqIO.write(record, faa_fl, "fasta")


def write_fna(feat, fna_fl, chrm, seq, attributes, nbgene):

    # print type(seqAA)
    # attributes = giveQualifiersInfo(feat.qualifiers, ["locus_tag", 'protein_id', 'product'], sep = ' ')
    record = SeqRecord(seq, id="{}|{} ".format(chrm, str(nbgene)),
                       name="", description=' '.join(attributes))
    SeqIO.write(record, fna_fl, "fasta")


def see_objet(obj):
    print 'SEE OBJET ', type(obj)

    for attr in dir(obj):
        if attr[0] != '_':
            print attr, " ", getattr(obj, attr)


def from_gb_to_required_format(data_path_base, gb_file):
    gff_fl = open(data_path_base + '.gff', 'w')
    faa_fl = open(data_path_base + '.faa', 'w')
    fna_fl = open(data_path_base + '.fna', 'w')

    input_handle = open(gb_file, "rU")
    for record in SeqIO.parse(input_handle, "genbank"):
        # see_objet(record)
        SeqIO.write(record, data_path_base+".fasta", "fasta")
        i = 0
        CDS_not_translated = 0
        for f in record.features:
            # if i > 150:
            #     break
            if f.type == 'CDS':  # and 'translation' in f.qualifiers:
                i += 1
                seq = f.extract(record.seq)
                attributes = giveQualifiersInfo(
                    f.qualifiers, ["old_locus_tag", 'protein_id', "locus_tag", 'product'])
                if write_faa(f, faa_fl, record.id, seq,  attributes, i) is not False:
                    write_gff(f, gff_fl, record.id, attributes, i)
                    write_fna(f, fna_fl, record.id, seq, attributes, i)
                else:
                    CDS_not_translated += 1
        print('{} CDS have encountered a problem during translation over {}'.format(CDS_not_translated, i))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='Prepare_Data_from_gbfile', description='Prepare data to be analysed by MeTAfisher starting with the genbank file of the sequence')

    parser.add_argument("gb_file")

    args = parser.parse_args()
    gb_file = args.gb_file
    if path.isfile(path.join('~', gb_file)):
        raise ValueError('The provided file ({}) does not exist'.format(gb_file))
    if gb_file[-3:] != ".gb":
        raise ValueError(
            'The provided file ({}) does not have the correct extension .gb'.format(gb_file))

    # directory where formated data are stored is the one of the gbfile...
    data_path_base = path.abspath(gb_file).split(".gb")[0]

    from_gb_to_required_format(data_path_base, gb_file)
