from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def giveQualifiersInfo(qualifiers, keys, sep=';'):
    line = ''
    for key in keys:
        if key in qualifiers:
            line+= "{}={}{}".format(key, qualifiers[key][0], sep)
    return line

def write_gff(feat, fl_gff, chrm, nbgene):
    strand = '+' if feat.strand == 1 else '-'
    attributes = "{}|{};".format(chrm,str(nbgene))
    attributes += giveQualifiersInfo(feat.qualifiers, ["locus_tag", 'protein_id', 'product'])

    list_gff =[chrm, 'Genbank', feat.type, str(feat.location.start.real), str(feat.location.end.real), ".", strand, ".",attributes+'\n']
    fl_gff.write('\t'.join(list_gff))

def write_faa(feat, faa_fl, chrm, seq, nbgene):
    trans_table = int(feat.qualifiers['transl_table'][0])

    # seqAA = seq.translate(table=trans_table)
    seqAA = Seq(feat.qualifiers['translation'][0])
    # print type(seqAA)
    attributes = giveQualifiersInfo(feat.qualifiers, ["locus_tag", 'protein_id', 'product'], sep = ' ')
    record = SeqRecord(seqAA, id="{}|{} ".format(chrm,str(nbgene)), name="", description=attributes)
    SeqIO.write(record, faa_fl, "fasta")

def write_fna(feat, fna_fl, chrm, seq, nbgene):

    # print type(seqAA)
    attributes = giveQualifiersInfo(feat.qualifiers, ["locus_tag", 'protein_id', 'product'], sep = ' ')
    record = SeqRecord(seq, id="{}|{} ".format(chrm,str(nbgene)), name="", description=attributes)
    SeqIO.write(record, fna_fl, "fasta")


def see_objet(obj):
    print 'SEE OBJET ', type(obj)

    for attr in dir(obj):
        if attr[0] != '_':
            print attr, " ", getattr(obj, attr)



gk_file = '../data/Acaryochloris_marina_MBIC11017/sequence.gb'
gk_file = '../data/Acaryochloris_marina_MBIC11017/sequencenofull.gb'
data_way = '../data/Acaryochloris_marina_MBIC11017/'

gff_fl = open(data_way + 'sequenceNO.gff', 'w')
faa_fl = open(data_way + 'sequenceNO.faa', 'w')
fna_fl = open(data_way + 'sequenceNO.fna', 'w')

input_handle = open(gk_file, "rU")
for record in SeqIO.parse(input_handle, "genbank"):
    # see_objet(record)
    # SeqIO.write(record, data_way+"sequence.fasta", "fasta")
    i=0
    c = 0
    for f in record.features:
        # if i > 150:
        #     break
        if f.type == 'CDS' and 'translation' in f.qualifiers:
            i += 1
            seq = f.extract(record.seq)
            write_gff(f,gff_fl, record.id, i)
            write_faa(f, faa_fl, record.id, seq, i)
            write_fna(f, fna_fl, record.id, seq, i)

            if 'tox' in f.qualifiers['product'][0]:
                print f.qualifiers
                c += 1
print c
