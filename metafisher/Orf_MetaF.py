#!/usr/bin/env python3

"""
Module      : Orf
Description : Orf specific functions.
Copyright   : (c) Jean Mainguy, 27 nov. 2020
License     : MIT
Maintainer  : jean.mainguy@outlook.fr


Tool to retrieve toxin antitoxin (TA) systems in genomes or metagenomes.
"""


from itertools import chain
import Object_MetaF as obj
import Function_MetaF as fct
import logging
import os


# def rescue_lonely_genes(orf_dict, gff_dict, contig, tmp_adjorf_faa, genes):
def get_adjacent_orfs(orf_dict, gff_dict, contig, genes):
    """
    Manage all the step of rescuing lonely gene.

    Argument :
        orf_dict: has the contig sequence file open and the current line. allow to retrieve the sequence
        gff_dict: gff file open and curent line. used by get_gff_ends() to get stop postion of the predicted gene in order
        to not take into acount the related orf.
        contig: name of the current contig
        tmp_adjorf_faa: name of the file where adj orf will be write
    """

    lonely_genes_on_strand_plus = [gene for gene in fct.get_lonely_genes(genes) if gene.strand == "+"]
    lonely_genes_on_strand_minus = [gene for gene in fct.get_lonely_genes(genes) if gene.strand == "-"]

    genes_plus = [gene for gene in genes if gene.strand == '+']
    genes_minus = [gene for gene in genes if gene.strand == '-']

    # Extraction of the sequence of the contig in order to find the orf inside.
    contig_seq, orf_dict["line"] = fct.get_fast_fasta(orf_dict['fl'], orf_dict['line'], contig)
    obj.Orf.seq = contig_seq
    # retrieve the end position of the predicted gene to then skip the related ORF:
    gff_ends = get_gff_ends(gff_dict, contig)
    adj_orfs = []
    # TEST IF THERE ARE LONELY GENE IN STRAND plus
    if lonely_genes_on_strand_plus:
        orf_generator = findORF(contig, contig_seq, rev=1)
        # read the orfs: check if the orf is a predicted gene
        # then check if the orf is adjacent with a TA_gene
        # If so then the faa ORF file is written...
        # orf manager will do everything !!

        adj_orfs += get_adjacent_orfs_by_strand(orf_generator, '+', genes_plus, gff_ends)
    # TEST IF THERE ARE LONELY GENES IN STRAND minus
    if lonely_genes_on_strand_minus:
        orf_generator = findORF(contig, contig_seq, rev=-1)

        adj_orfs += get_adjacent_orfs_by_strand(orf_generator, '-', genes_minus, gff_ends)
    return adj_orfs

def identify_ta_orfs(adj_orfs, genes, outdir):
    adj_orf_dict = {}
    gene_index = max([gene.gene_number for gene in genes])

    logging.info(f'{len(adj_orfs)} adjacent orfs have been found! ')

    adjorf_faa_file = os.path.join(outdir , 'temporary_adjOrf.faa')
    with open(adjorf_faa_file, 'w') as fl:
        for i, orf in enumerate(adj_orfs):
            orf.write_faa(fl, i)
            adj_orf_dict[i] = orf

    # launch hmmsearch
    hmm_result_file = os.path.join(outdir , 'output_HMM_table_adj_orf.txt')
    table_hmm = fct.HMM_launcher(adjorf_faa_file, hmm_result_file)
    # parsing result and store it in obj.Orf.hmm_orf
    ta_orfs = adjOrf_HMM(table_hmm, adj_orf_dict)
    logging.info(f'{len(ta_orfs)} have TA domains')
    # taking into account the domain in the possible start of the ORF
    # give gene number to hmm orf that fit with genes from gff file
    adjust_orf_attribut(ta_orfs, gene_index)

    # Adding adjOrf to list of TA gene and finding with who there are adjacent
    hmm_orf_get_adj(ta_orfs, genes)

    return ta_orfs


def hmm_orf_get_adj(hmm_orfs, genes):
    """Get adjacent gene."""

    for orf in hmm_orfs:
        strand = orf.strand

        for gene in (gene for gene in genes if gene.strand == strand):
            # if the orf have been already process and added to the list then we don't check it against itself
            if orf.gene_number == gene.gene_number:
                continue

            if gene.is_pre_adj_to(orf):
                add_post_pre_attr(strand, g_post=orf, g_prev=gene)

            elif orf.is_pre_adj_to(gene):
                add_post_pre_attr(strand, g_post=gene, g_prev=orf)

        # obj.TA_gene.genes_strand[strand].append(orf)
        genes.append(orf)



def add_post_pre_attr(strand, g_post, g_prev):
    """Add post and pre gene to gene object."""
    # obj.TA_gene.linked.add(g_post)
    # obj.TA_gene.linked.add(g_prev)
    if strand == '+':
        g_prev.post.append(g_post)
        g_post.prev.append(g_prev)
    else:  # in minus strand the post rev are inverted
        g_post.post.append(g_prev)
        g_prev.prev.append(g_post)


def adjust_orf_attribut(orfs, gene_index):
    """
    Adjust possible start to domain.

    Give post prev attribute.
    Give gene number to hmm orf that fit with genes from gff file.
    """
    for orf in orfs:
        orf.possible_start = [s for s in orf.possible_start if s <= orf.domain_Ct_border]
        orf.distanceMin = obj.Gene.distanceMin - abs(orf.possible_start[-1] - orf.possible_start[0])

        # WARNING absolute value of distance min should not be greater than the length of the gene !!
        # because it would allow overlap of more than the length of the gene
        # and then could give a post gene before a pre gene so a non sense
        if abs(orf.distanceMin) >= len(orf):
            orf.distanceMin = -len(orf) + 1  # if so distance min is -abs(len -1)

        orf.post = []
        orf.prev = []
        gene_index += 1
        orf.gene_number = gene_index


def codon_finder(liste, seq, frame=1, inf=0, sup='Not defined'):
    """Extract position of given codons in a sequence and in a given frame

    By default this function search codon of the given list in the whole sequence. However it possible to give the position in the sequence where codon has to be found with sup and inf parameter.

    Args:
        liste   : liste of codons
        seq     : DNA sequence
        frame   : give the frame of the search
        inf     : position in sequence where the function has to start to search
        sup     : position in sequence where the function has to stop to search

    Returns:
        Return a list of the position of the given codon in the sequence
    """
    position = []
    if sup == 'Not defined':
        sup = len(seq) - 2
    for i in range(inf + frame - 1, sup, 3):  # -1 to get frame from 0 to 2 and not from 1 to 3
        if seq[i:i + 3] in liste:
            position.append(i)
    return position

def get_adjacent_orfs_by_strand(orfs, strand, lonelyGenes, gff_ends):
    """
    Manage ORF.

    Lonely gene are a list of lonely gene from the same strand
    check to see if orf is a predicted gene.
    """
    adj_orfs = []
    for orf in orfs:
        if orf.real_end() in gff_ends[strand]:
            gff_ends[strand].remove(orf.real_end())
            continue
        # EDIT now orfinder size max and min correctly...
        orf.distanceMin = obj.Gene.distanceMin - (orf.possible_start[-1] - orf.possible_start[0])
        # WARNING absolute value of distance min should not be greater than the length of the gene !!
        # because it would allow overlap of more than the length of the gene
        # and then could give a post gene before a pre gene so a non sense

        if abs(orf.distanceMin) >= len(orf):
            orf.distanceMin = -len(orf) + 1  # if so distance min is len -1

        if orf.is_adj(lonelyGenes):
            adj_orfs.append(orf)

    if gff_ends[strand]:
        logging.warning(orf.scaffold)
        logging.warning(strand)
        logging.warning(f"len du saffold {len(orf.seq['data'])}")
        logging.warning("THERE ARE SOME PEDICTED GENE THAt DIDNT FOUND THEIR ORF :-(")

        logging.warning("scaffold " + orf.scaffold + ' strand :' + strand + '\n' + str(gff_ends[strand]))
        for g in gff_ends[strand]:
            logging.warning(str(g) + '\n')
    return adj_orfs

def adjOrf_HMM(table_hmm, adj_orf_dict):
    """Do something."""
    hmm_orfs = []
    with open(table_hmm, 'r') as fl:
        for line in fl:
            if line.startswith("#"):
                continue
            gene_number, domain = fct.hmmtable_parser(line)
            hmm_orf = adj_orf_dict[gene_number]

            try:
                hmm_orf.domain.append(domain)
                hmm_orf.domain_Ct_border = max(hmm_orf.domain_Ct_border, domain.ali_from * 3)
            except AttributeError:
                hmm_orf.domain_Ct_border = domain.ali_from * 3
                hmm_orf.domain = [domain]

                # obj.Orf.hmm_orf.setdefault(hmm_orf.strand, []).append(hmm_orf)
                hmm_orfs.append(hmm_orf)
    return hmm_orfs

def get_gff_ends(gff_dict, scaffold):
    """
    Get gff ends.

    gff_dict is composed of the line and the fl of the gff file
    Work in the same way as get_fast_fasta
    Yield the gff stop according the frame
    Then return a dictionary composed of the orfs by frame
    """
    gff_ends = {'+': [], '-': []}
    line = gff_dict["line"]
    fl = gff_dict["csv"]
    while line and line[0] != scaffold:
        line = next(fl)

    while line and line[0] == scaffold:
        if line[0] != scaffold:
            break
        if line[2] not in ['CDS', 'gene']:
            try:
                line = next(fl)
                continue
            except StopIteration:
                break
        # this number is used when we write the gff into file to have a number that doesn't overlap with other gene in round2
        # highestGeneNumber = int(line[8].split(";")[0].split('|')[1])

        if line[6] == '+':
            gff_ends['+'].append(int(line[4]))
        else:
            gff_ends['-'].append(int(line[3]))
        try:
            line = next(fl)
        except StopIteration:
            break

    gff_dict["line"] = line
    return gff_ends


def findORF(scaffold, seq, rev):
    """Determine the orf of a nucleic sequence.

    This function call the orf_by_frame() function  for each frame and give it to it different parameter
    especially the frame and the information if the sequence is reverse or not and finally gather the result in a dictionary

    Args:
        seq : sequence where the orf will be determine.
        threshold : the minimum size of orf
        nub_genTable : id Genetic code Table
    Returns:
        A list of ORFs in the class ListGene
    """
    threshold = obj.Orf.length_min
    starts = obj.Gene.codon_starts
    stops = obj.Gene.codon_stops

    if (rev == -1 and seq["rev"] is False) or (rev == 1 and seq['rev'] is True):
        seq = complement_reverse(seq)

    generator1 = orf_by_frame(seq, threshold, starts, stops, frame=1, rev=rev)
    generator2 = orf_by_frame(seq, threshold, starts, stops, frame=2, rev=rev)
    generator3 = orf_by_frame(seq, threshold, starts, stops, frame=3, rev=rev)

    return chain(generator1, generator2, generator3)


def orf_by_frame(seq, threshold, starts, stops, frame, rev):
    """Find the orf in a sequence for a given frame.

    First the function determine all the stop codon of the sequence in frame by calling codon_finder().
    Then it will find all the start codon flan by two stop codon by calling codon_finder
    And finally will build the dictionaries with the iformation relative to each orf by calling the remplissage() function

    Args:
        seq : sequence where the orf will be determine
        frame : give the frame of the search
        threshold : the minimum size of orf
        orf_dict: store start and stop codon of the geneic code 11
        rev : 1 or -1 to know if the sequence has been completed reverted in order to fill correctly the frame.
        EDIT: reverse info is stored in seq dictionary

    Returns:
        The list of the Orf of the frame as a list of dictionaries
    """
    length_seq = len(seq['data'])
    pos_stop = codon_finder(stops, seq['data'], frame=frame)
    if not pos_stop:  # Check the case if there is no stop codon in te sequence
        logging.debug("NO STOP IN THE SEQUENCE")
        if length_seq > threshold:
            twoStopLen = length_seq
            twoStopLen = twoStopLen - twoStopLen % 3
            add_to_inf = to_fit_len_max(twoStopLen)

            pos_start = []
            # to have a len that is a multiple of 3
            add_to_inf = to_fit_len_max(twoStopLen=length_seq)
            # research of starts in the whole sequence !
            pos_start = codon_finder(starts, seq['data'], inf=add_to_inf, frame=frame)
            # WARNING add at the beggining of the seq an imaginary star, necessary to stick to gene prediction
            if add_to_inf == 0:
                pos_start = [frame - 1] + pos_start

            if pos_start:
                # Search the correct stop position
                # distance between the start and the end of the MetaG
                distance = length_seq - pos_start[-1]
                # -3 is -1 to get index of the last nt because here we start counting at 0 and -2 because we need the first nt of the last codon
                stop = length_seq - 3 - (distance % 3)
                yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=stop + 2, complet=False, border=True)
        return

    # research of potential orf between 0 position and the first stop -->  the program doesn't manage circular DNA
    if pos_stop[0] + 3 > threshold:  # check if the length from 0 to the first stop is bigger than the trhresold
        pos_start = []
        twoStopLen = pos_stop[0] + 3
        twoStopLen = twoStopLen - twoStopLen % 3
        add_to_inf = to_fit_len_max(twoStopLen)  # +3 is to get the len between 0 to the first stop
        pos_start = codon_finder(starts, seq['data'], frame=frame, inf=add_to_inf, sup=(
            pos_stop[0] + 3 - threshold))  # research of starts from 0 to (first stop - thresold)

        if add_to_inf == 0:
            pos_start = [frame - 1] + pos_start
        if pos_start:
            yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=pos_stop[0] + 2, complet=True, border=True)

    for i in range(len(pos_stop) - 1):  # search between stop i and stop i+1. Stop before the end.

        # First +3 to include stop codin in final length as in the ncbi orfinder.
        if (pos_stop[i + 1] + 3 - (pos_stop[i] + 3)) < threshold:
            continue

        else:
            pos_start = []
            add_to_inf = to_fit_len_max(twoStopLen=pos_stop[i + 1] - pos_stop[i])
            pos_start = codon_finder(
                starts, seq['data'], inf=pos_stop[i] + 3 + add_to_inf, sup=(pos_stop[i + 1] + 3 - threshold))

            if pos_start:
                yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=pos_stop[i + 1] + 2, complet=True, border=False)

    # Search of orf from the last stop codon to the end of the seq
    if length_seq - pos_stop[-1] + 3 > threshold:
        pos_start = []
        twoStopLen = length_seq - pos_stop[-1] + 3
        twoStopLen = twoStopLen - twoStopLen % 3
        add_to_inf = to_fit_len_max(twoStopLen)  # %3 to get a multiple of 3
        # research of starts from 0 to (first stop - thresold)
        pos_start = codon_finder(
            starts, seq['data'], inf=pos_stop[-1] + 3 + add_to_inf, sup=(length_seq - threshold))
        if pos_start:
            # looking for the appropriate stop:
            # It is a bit tricky : The gene is at the border of the MetaG
            # We are looking for the position of the first nt
            # of the last full codon in the frame of the gene
            # distance between the start and the end of the MetaG
            distance = length_seq - pos_start[-1]
            # -3 is -1 to get index of the last nt because here we start counting at 0 and -2 because we need the first nt of the last codon
            stop = length_seq - 3 - (distance % 3)

            yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=stop + 2, complet=False, border=True)

    return


def to_fit_len_max(twoStopLen):
    """
    Determine the inferior limit for searching the start codons between two stop codons.

    Take the length between the two stop codon
    And give back the inferior borne where to begin the searc of start codon
    IF YOU WANT TO HAVE THE FULL ORF WITHOUT A MAX SIZE LIMITATION YOU CAN RETURN 0 ONLY
    """
    if twoStopLen <= obj.Gene.length_max:  # if the space between the two stop <= to len max then no need to add anything
        return 0
    return twoStopLen - obj.Gene.length_max


def complement_reverse(seq):
    """Complement reverse sequence."""
    compl_rev_seq = {}
    compl_rev_seq = reverse(seq)
    compl_rev_seq = complement(compl_rev_seq)
    compl_rev_seq['rev'] = not seq['rev']
    return compl_rev_seq


def reverse(seq):
    """Reverse sequence."""
    rev_seq = {}
    rev_seq['data'] = seq['data'][::-1]
    return rev_seq


def complement(seq):
    """Complement sequence."""
    result = ''
    complement = {}
    for i in seq['data']:
        if i == 'A':
            result = result + 'T'
        elif i == 'T':
            result = result + 'A'
        elif i == 'C':
            result = result + 'G'
        elif i == 'G':
            result = result + 'C'
    complement['data'] = result
    return complement


def getGeneticCode(ncbi_genetic_code_id):
    """Return the genetic code corresponding to a given ncbi_genetic_code_id.

    This function is written by Theo Falgarone.

    Args:
        ncbi_genetic_code_id: NCBI identifier for the genetic code

    Returns:
        table: genetic code table

    """
    if ncbi_genetic_code_id == 1:
        base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
        start = "---M---------------M---------------M----------------------------"
        aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    elif ncbi_genetic_code_id == 11:
        base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
        start = "---M---------------M------------MMMM---------------M------------"
        aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    else:
        raise ValueError(f'The given ncbi id {ncbi_genetic_code_id} is not managed yet.')
        return
    table = {}
    codonstart = []
    codonstop = []
    for i in range(0, len(aas)):
        codon = (base1[i] + base2[i] + base3[i])
        aa = aas[i]
        table[codon] = aa
        if aa == '*':
            codonstop.append(codon)
        elif start[i] == 'M':
            codonstart.append(codon)
    table['stop'] = codonstop
    table['start'] = codonstart
    return table


def translate(nucl_seq, codonTable=11):
    """Return the translated proteic sequence based on the standard genetic code."""
    proteic_seq = {}
    codons = getWords(nucl_seq, 3)

    code = getGeneticCode(codonTable)

    nucl_seq["data"] = nucl_seq["data"].upper()
    # proteic_seq["description"]="translated proteic sequence from "+nucl_seq["description"]
    proteic_seq["data"] = ""
    for codon in codons:
        if codon not in code:
            proteic_seq["data"] = proteic_seq["data"] + "?"
        else:
            proteic_seq["data"] = proteic_seq["data"] + code[codon]

    return proteic_seq


def getWords(seq, wlength):
    """Return a list of non-overlapping words of a given length in a seq."""
    seq["data"] = seq["data"].upper()
    words = []
    for i in range(0, len(seq["data"]) - wlength + 1, wlength):
        words.append(seq["data"][i:i + wlength])

    return words
