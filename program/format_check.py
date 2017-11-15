import sys
import re
import csv


def is_gff_ok(gff):
    # TODO
    return False


def is_faa_ok(faa):
    # TODO
    return False


def is_fna_ok(fna):
    # TODO
    return False


def get_contig_name(fasta_fl):
    """
    retrieve the contig names of the fasta file
    contig name is the first world of the header

    >CP000471.1 Magnetococcus marinus MC-1, complete genome
    here the conig_name would be CP000471.1

    """
    liste = []
    try:
        with open(fasta_fl, 'r') as fl:
            for l in fl:
                if l[0] == '>':
                    try:
                        index = l.index(" ")
                    except ValueError:
                        index = len(l)
                    liste.append(l[1:index])
            return liste
    except IOError:
        return []


def extract_position(header):

    if not header[0] == '>':
        raise ValueError("the line given to extract postion is not a header\n", header)
    pattern = "\[location=[complent(]{0,11}(\d+)..(\d+)[)]?\]"
    p = re.compile(pattern)
    s = p.search(header)
    if s:
        start = s.group(1)
        end = s.group(2)
    else:
        print header
        raise Exception("There is no [location=XXX...XXX] in the header. :-(")
    return [start, end]


def search_gff_line(postion, gff_handler):
    for l in gff_handler:
        if len(l) < 9:
            continue
        # if l[2] == 'gene' and l[3] == position[0] and l[4] == position[1]:
        #     l[8].find('locus_tag')
        if l[2] == 'gene' and l[3] == position[0] and l[4] == position[1]:  # check if position of fna match with the line
            return l


def rewrite_seq(seqfile, newfl, gene_counter, line, chrm_name):
    header = '>{}|{} {}'.format(chrm_name, gene_counter, line[1:])
    newfl.write(header)
    for l in seqfile:
        if l[0] == '>':
            return l
        newfl.write(l)
    return ''  # when the seqfile is finish the return line is empty then the while loop will be able to end


def rewrite_gff(gff_newfl, gff_line, gene_counter):
    gff_line[8] = 'ID={}|{};{}'.format(gff_line[0], gene_counter, gff_line[8])
    gff_newfl.write("\t".join(gff_line) + '\n')


if __name__ == '__main__':
    """
    Normally the files are sorted according the postion of the genes in the genome.
    And also they are sorted in the same way between the different files
    """
    # print sys.argv  # ${data_pathway} ${data_name} ${output_pathway}
    data_pathway = sys.argv[1] + '/' + sys.argv[2]  # ${data_pathway}/${data_name}
    outpath = sys.argv[3] + '/' + sys.argv[2]  # output_pathway + data_name
    fasta_file = data_pathway + '.fasta'
    gff_file = data_pathway + '.gff'
    faa_file = data_pathway + '.faa'

    fna_file = data_pathway + '.fna'
    fna_existence = True

    contig_names = get_contig_name(fasta_file)

    faa_fl = open(faa_file, 'r')
    gff_fl = csv.reader(open(gff_file, 'r'), delimiter='\t')
    try:
        fna_fl = open(fna_file, 'r')
    except IOError:
        fna_existence = False

    faa_newfl = open(outpath + '.faa', 'w')

    gff_newfl = open(outpath + '.gff', 'w')
    gff_newfl.write("##GFF_file rewrite to match MeTAfisher output\n")
    gene_counter = 0
    faa_line = faa_fl.readline()
    if fna_existence:
        fna_newfl = open(outpath + '.fna', 'w')
        fna_line = fna_fl.readline()

    while faa_line:  # while there is line to read in the fna and faa files do..
        gene_counter += 1
        position = extract_position(faa_line)  # from the header

        gff_line = search_gff_line(position, gff_fl)
        if gff_line[0] not in faa_line:
            raise Exception("The contig name from gff line seems to be not the same in the header of faa file:\nGFF:\n{}\nFAA\n{}".format(gff_line[0], faa_line))
        if contig_names and gff_line[0] not in contig_names:  # if fasta file don' exist then contig name is an empty list so we skip this condition
            raise Exception("Contig or chromosome name from gff file don't match the name of the header of the chromomosome sequence(.fasta file)")

        faa_line = rewrite_seq(faa_fl, faa_newfl, gene_counter, faa_line, chrm_name=gff_line[0])
        rewrite_gff(gff_newfl, gff_line, gene_counter)

        if fna_existence:
            if not position == extract_position(fna_line):
                print "position from the fna file and faa file do no match"
                print "\nFNA line:\nfna_line\n\nFAA line:\nfaa_line\n"
                raise Exception('position from the fna file and faa file do not match')
            fna_line = rewrite_seq(fna_fl, fna_newfl, gene_counter, fna_line, chrm_name=gff_line[0])
        # print 'CONTIG NAME', gff_line[0]
        # print faa_line
