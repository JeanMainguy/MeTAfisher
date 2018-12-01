import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from os import path


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


def from_gb_to_required_format(data_way, gb_file):
    gff_fl = open(data_way + '.gff', 'w')
    faa_fl = open(data_way + '.faa', 'w')
    fna_fl = open(data_way + '.fna', 'w')

    input_handle = open(gb_file, "rU")
    for record in SeqIO.parse(input_handle, "genbank"):
        # see_objet(record)
        SeqIO.write(record, data_way+".fasta", "fasta")
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
    try:
        gb_file = sys.argv[1]  # 'data/Acaryochloris_marina_MBIC11017/sequence.gb'
    except IndexError:
        raise IndexError(
            "Please provide the genbank file when you launch the script: `python genbank_parser.py <genbankFile.gb>`")
    if path.isfile(gb_file) and gb_file[-3:] != ".gb":
        raise ValueError('The provided file does not have the correct extension .gb')
    data_way = gb_file.split(".gb")[0]  # 'data/Acaryochloris_marina_MBIC11017/'

    from_gb_to_required_format(data_way, gb_file)
