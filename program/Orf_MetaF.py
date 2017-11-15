# coding: utf-8
from itertools import chain
import Object_MetaF as obj
from subprocess import call
import Function_MetaF as fct2
# import re
# import csv
# from operator import attrgetter, itemgetter
# import find_orf as orf
# from math import log


def rescue_lonely_gene(dico_orf, dico_gff, scaffold, tmp_adjorf_faa):
    """
    Manage all the step of rescuing lonely gene
    Argument :
        dico_orf: has the contig sequence file open and the current line. allow to retrieve the sequence
        dico_gff: gff file open and curent line. used by get_gff_ends() to get stop postion of the predicted gene in order
        to not take into acount the related orf.
        scaffold: name of the current scaffold
        tmp_adjorf_faa: name of the file where adj orf will be write
    """
    fl = open(tmp_adjorf_faa, 'w')
    # Extraction of the sequence of the contig in order to find the orf inside.
    seq, dico_orf["line"] = fct2.get_fast_fasta(dico_orf['fl'], dico_orf['line'], scaffold)
    obj.Orf.seq = seq
    # retrieve the end position of the predicted gene to then skip the related ORF:
    gff_ends, highestGeneNumber = get_gff_ends(dico_gff, scaffold)

    # TEST IF THERE ARE LONELY GENE IN STRAND plu
    if obj.TA_gene.lonely['+']:
        # print 'THERE IS LONELY GENE ON STRAND +'
        generator = findORF(scaffold, seq, rev=1)
        # read the generator: check if the orf is a predicted gene
        # then check if the orf is adjacent with a TA_gene
        # If so then the faa ORF file is written...
        orf_manager(generator, '+', obj.TA_gene.genes_strand['+'], gff_ends, fl)  # orf manager will do everything !!
    # TEST IF THERE ARE LONELY GENES IN STRAND minus
    if obj.TA_gene.lonely['-']:
        # print 'THERE IS LONELY GENE ON STRAND -'
        generator = findORF(scaffold, seq, rev=-1)
        orf_manager(generator, '-', obj.TA_gene.genes_strand['-'], gff_ends, fl)
    # else:
        # raw_input(scaffold + " THERE IS NO LONELY GENE ON STRAND -")
    fl.close()
    # print obj.Orf.adj_orf
    # If the program have found adj orf then several step are applied
    if obj.Orf.adj_orf:
        # launch hmmsearch
        table_hmm = fct2.HMM_launcher(tmp_adjorf_faa, add_to_name='adj_orf')
        # parsing result and store it in obj.Orf.hmm_orf
        adjOrf_HMM(table_hmm)
        # taking into account the domain in the possible start of the ORF
        # give gene number to hmm orf that fit with genes from gff file
        adjust_orf_attribut(obj.Orf.hmm_orf, highestGeneNumber)

        # Adding adjOrf to list of TA gene and finding with who there are adjacent
        hmm_orf_get_adj(obj.Orf.hmm_orf)


def hmm_orf_get_adj(dico_hmmorf):
    # Add the hmm orf into the main lists !
    # obj.TA_gene.genes.extend(dico_hmmorf['+'] + dico_hmmorf['-'])
    # obj.TA_gene.genes_plus.extend(dico_obj['+'])
    # obj.TA_gene.genes_minus.extend(dico_obj['-'])

    for strand in dico_hmmorf:
        for g in dico_hmmorf[strand]:
            print g.gene_number
        # Adding adjOrf to the main list of TA gene
        # obj.TA_gene.genes_strand[strand].extend(dico_hmmorf[strand])
        # obj.TA_gene.genes.extend(dico_hmmorf[strand])
        for hmmorf in dico_hmmorf[strand]:
            for gene in obj.TA_gene.genes_strand[strand]:
                if hmmorf.gene_number == gene.gene_number:  # MAGIC METHOD
                    # print hmmorf
                    # print gene
                    continue
                if gene.is_pre_adj_to(hmmorf):
                    add_post_pre_attr(strand, g_post=hmmorf, g_prev=gene)
                    # obj.TA_gene.linked.add(gene)
                    # obj.TA_gene.linked.add(hmmorf)
                    # if strand == '+':
                    #     gene.post.append(hmmorf)
                    #     hmmorf.prev.append(gene)
                    # else:
                    #     hmmorf.post.append(gene)
                    #     gene.prev.append(hmmorf)
                elif hmmorf.is_pre_adj_to(gene):
                    add_post_pre_attr(strand, g_post=gene, g_prev=hmmorf)
                    # obj.TA_gene.linked.add(gene)
                    # obj.TA_gene.linked.add(hmmorf)
                    # if strand == '+':
                    #     hmmorf.post.append(gene)
                    #     gene.prev.append(hmmorf)
                    # else:
                    #     gene.post.append(hmmorf)
                    #     hmmorf.prev.append(gene)
            obj.TA_gene.genes_strand[strand].append(hmmorf)
            obj.TA_gene.genes.append(hmmorf)


def add_post_pre_attr(strand, g_post, g_prev):
    # add to the set linked
    obj.TA_gene.linked.add(g_post)
    obj.TA_gene.linked.add(g_prev)
    if strand == '+':
        g_prev.post.append(g_post)
        g_post.prev.append(g_prev)
    else:  # in minus strand the post rev are inverted
        g_post.post.append(g_prev)
        g_prev.prev.append(g_post)


def adjust_orf_attribut(dico_obj, highestGeneNumber):
    """
    adjust_possible_start_to_domain
    give post prev attr
    give gene number to hmm orf that fit with genes from gff file
    """
    for strand in dico_obj:
        for o in dico_obj[strand]:
            # print 'CT BORDER', o.domain_Ct_border
            # print 'before', o.possible_start
            # starts_zip = zip(o.possible_start, o.possible_start_orf)
            o.possible_start = [s for s in o.possible_start if s <= o.domain_Ct_border]
            # print 'after', o.possible_start
            o.distanceMin = obj.Gene.distanceMin - abs(o.possible_start[-1] - o.possible_start[0])
            o.post = []
            o.prev = []
            highestGeneNumber += 1
            o.gene_number = highestGeneNumber


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
    for i in xrange(inf + frame - 1, sup, 3):  # -1 to get frame from 0 to 2 and not from 1 to 3
        if seq[i:i + 3] in liste:
            position.append(i)
    return position


def orf_manager(generator, strand, lonelyGenes, gffEnds, fl):
    """
    lonely gene are a list of lonely gene from the same strand
    check to see if orf is a predicted gene
    """
    # print 'STRAND ', strand
    # print 'GFF ENDS\n', gffEnds[strand]
    for o in generator:
        if o.real_end() in gffEnds[strand]:
            # print 'Is a predicted gene..', o
            gffEnds[strand].remove(o.real_end())
            continue
        # resize if too long.. Not possible to have something to small
        # because orfinder has already a minimum size limit tha fit the threshold
        # set also the g.distanceMin personalized of the orf to get the tandem gene.
        # EDIT now orfinder size max and min correctly...
        # o.manage_size()  # A enlever quand tu le sens attention au distanceMin ! !
        o.distanceMin = obj.Gene.distanceMin - (o.possible_start[-1] - o.possible_start[0])
        # .distanceMin = Gene.distanceMin - abs(o.possible_start[-1] - o.possible_start[0])
        if o.is_adj(lonelyGenes):
            # print 'IS adj', o
            obj.Orf.adj_orf_index += 1
            obj.Orf.adj_orf[obj.Orf.adj_orf_index] = o
            # print obj.Gene.distanceMin, '-(', (o.possible_start[-1], '-', o.possible_start[0]), ') = ', o.distanceMin
            # print 'Index', obj.Orf.adj_orf_index, 'distance min', o.distanceMin
            # try:
            #     obj.Orf.adj_orf.setdefault(strand, []).append(o)
            # except AttributeError:
            #     obj.Orf.adj_orf = {'+': [], '-': []}
            #     obj.Orf.adj_orf[strand].append(o)
            o.write_faa(fl, index=obj.Orf.adj_orf_index)  # write the file to tmp faa file
    # ==============================================
        # #REMOVE##
        # if o.strand == '+' and o.real_end() > 24000:
        #     print 'OLALALA'
        #     print o, 'frame :', o.frame, 'real end', o.real_end()
        ####
    if gffEnds[strand]:
        print o.scaffold
        print strand
        print 'len du saffold', len(o.seq['data'])
        print "THERE ARE SOME PEDICTED GENE THAt DIDNT FOUND THEIR ORF :-("
        # TODO LOG CHANGE THAT
        with open(obj.Gene.output_way + 'Predicted_not_found.err', 'a') as fl:
            fl.write("scaffold " + o.scaffold + ' strand :' + strand + '\n' + str(gffEnds[strand]))
            for g in gffEnds[strand]:
                fl.write(str(g) + '\n')
        # print gffEnds[strand]
    # if strand == "+":
    #     with open(obj.Gene.output_way + 'sca_len.txt', 'a') as fl:
    #         fl.write(o.scaffold + '   ' + str(len(o.seq['data'])) + '\n')
    #         fl.write(o.seq['data'][len(o.seq['data']) - 20:len(o.seq['data'])])


def adjOrf_HMM(table_hmm):
    """

    """

    # tablehmm = obj.Gene.output_way + '/table_hmm_adjorf.txt'
    # hmm_db = obj.Gene.hmmdb
    #
    # bash_commande = "hmmsearch -E 0.5  --domtblout {}/output/MP0313/test/adjOrf_table_hmm.txt  {}/dependence/ALL_plus_MET_curatted.hmm {}/data/MP0313/MP0313.549.faa".format(home, home, home)
    # print bash_commande
    #
    # call(bash_commande, shell=True)
    # call(obj.Orf.hmm_adjorf_launcher, shell=True)

    with open(table_hmm, 'r') as fl:
        for line in fl:
            if line[0] != "#":
                # print line
                gene_number, domain = fct2.hmmtable_parser(line)
                hmm_orf = obj.Orf.adj_orf[gene_number]
                # hmm_orf.gene_number = obj.Gene.highestGeneNumber

                try:
                    hmm_orf.domain.append(domain)
                    hmm_orf.domain_Ct_border = max(hmm_orf.domain_Ct_border, domain.ali_from * 3)
                except:
                    hmm_orf.domain_Ct_border = domain.ali_from * 3
                    hmm_orf.domain = [domain]
                    obj.Orf.hmm_orf.setdefault(hmm_orf.strand, []).append(hmm_orf)

                # Append to hmm_orf dico of class Attribute


def get_gff_ends(dico_gff, scaffold):
    '''
    dico_gff is composed of the line and the fl of the gff file
    Work in the same way as get_fast_fasta
    Yield the gff stop according the frame
    Then return a dictionary composed of the generator by frame
    '''
    gff_ends = {'+': [], '-': []}
    line = dico_gff["line"]
    fl = dico_gff["csv"]
    while line and line[0] != scaffold:
        line = next(fl)
    for line in fl:
        if line[0] != scaffold:
            break
        if line[2] not in ['CDS', 'gene']:
            line = next(fl)
            continue
        # this number is used when we write the gff into file to have a number that doesn't overlap with other gene in round2
        highestGeneNumber = int(line[8].split(";")[0].split('|')[1])

        if line[6] == '+':
            gff_ends['+'].append(int(line[4]))
        else:
            gff_ends['-'].append(int(line[3]))
    dico_gff["line"] = line
    return gff_ends, highestGeneNumber


def findORF(scaffold, seq, rev):
    """This function determine the orf of a nucleic sequence
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
    # print 'IN FIND ORF...'
    if (rev == -1 and seq["rev"] is False) or (rev == 1 and seq['rev'] is True):
        # print 'Pre rev?', seq["rev"], 'rev=', rev
        seq = complement_reverse(seq)
        # print 'Post rev? ', seq["rev"], 'rev=', rev
    # generator_dict = {}
    # for frame in [1, 2, 3]:
    #     generator_dict[frame * rev] = orf_by_frame(seq, threshold, starts, stops, frame=frame, rev=rev)
    # return generator_dict
    generator1 = orf_by_frame(seq, threshold, starts, stops, frame=1, rev=rev)
    generator2 = orf_by_frame(seq, threshold, starts, stops, frame=2, rev=rev)
    generator3 = orf_by_frame(seq, threshold, starts, stops, frame=3, rev=rev)

    return chain(generator1, generator2, generator3)


def orf_by_frame(seq, threshold, starts, stops, frame, rev):
    """This function find the orf in a sequence for a given frame

    First the function determine all the stop codon of the sequence in frame by calling codon_finder().
    Then it will find all the start codon flan by two stop codon by calling codon_finder
    And finally will build the dictionaries with the iformation relative to each orf by calling the remplissage() function

    Args:
        seq : sequence where the orf will be determine
        frame : give the frame of the search
        threshold : the minimum size of orf
        dico_orf: store start and stop codon of the geneic code 11
        rev : 1 or -1 to know if the sequence has been completed reverted in order to fill correctly the frame.
        EDIT: reverse info is stored in seq dictionary

    Returns:
        The list of the Orf of the frame as a list of dictionaries

        """
    length_seq = len(seq['data'])
    pos_stop = codon_finder(stops, seq['data'], frame=frame)
    if not pos_stop:  # Check the case if there is no stop codon in te sequence
        # TODO LOG
        print "NO STOP IN THE SEQUENCE"
        if length_seq > threshold:
            twoStopLen = length_seq
            twoStopLen = twoStopLen - twoStopLen % 3
            add_to_inf = to_fit_len_max(twoStopLen)

            pos_start = []
            add_to_inf = to_fit_len_max(twoStopLen=length_seq)  # to have a len that is a multiple of 3
            pos_start = codon_finder(starts, seq['data'], inf=add_to_inf, frame=frame)  # research of starts in the whole sequence !
            # ATTENTION ajout du debut de la sequence aux start... start imaginaire donc, mais necessaire pour coller au gene predit ! :-(
            if add_to_inf == 0:
                pos_start = [frame - 1] + pos_start

            if pos_start:
                # Search the correct stop position
                distance = length_seq - pos_start[-1]  # distance between the start and the end of the MetaG
                stop = length_seq - 3 - (distance % 3)  # -3 is -1 to get index of the last nt because here we start counting at 0 and -2 because we need the first nt of the last codon
                yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=stop + 2, complet=False, border=True)
        return

    # research of potential orf between 0 position and the first stop -->  the program doesn't manage circular DNA
    if pos_stop[0] + 3 > threshold:  # check if the length from 0 to the first stop is bigger than the trhresold
        pos_start = []
        twoStopLen = pos_stop[0] + 3
        twoStopLen = twoStopLen - twoStopLen % 3
        add_to_inf = to_fit_len_max(twoStopLen)  # +3 is to get the len between 0 to the first stop
        pos_start = codon_finder(starts, seq['data'], frame=frame, inf=add_to_inf, sup=(pos_stop[0] + 3 - threshold))  # research of starts from 0 to (first stop - thresold)
        # ATTENTION ajout du debut de la sequence aux start... start imaginaire donc, mais necessaire pour coller au gene predit ! :-(
        if add_to_inf == 0:
            pos_start = [frame - 1] + pos_start
        if pos_start:
            yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=pos_stop[0] + 2, complet=True, border=True)

    for i in range(len(pos_stop) - 1):  # recherche entre le stop i et i+1 donc on s'arrete a un stop avant la fin

        if (pos_stop[i + 1] + 3 - (pos_stop[i] + 3)) < threshold:  # premier +3 pour inclure le codon stop dans la longueur final comme dans ncbi orf finder...
            continue

        else:
            pos_start = []
            add_to_inf = to_fit_len_max(twoStopLen=pos_stop[i + 1] - pos_stop[i])
            pos_start = codon_finder(starts, seq['data'], inf=pos_stop[i] + 3 + add_to_inf, sup=(pos_stop[i + 1] + 3 - threshold))

            if pos_start:
                yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=pos_stop[i + 1] + 2, complet=True, border=False)

    # Search of orf from the last stop codon to the end of the seq
    if length_seq - pos_stop[-1] + 3 > threshold:
        pos_start = []
        twoStopLen = length_seq - pos_stop[-1] + 3
        twoStopLen = twoStopLen - twoStopLen % 3
        add_to_inf = to_fit_len_max(twoStopLen)  # %3 to get a multiple of 3
        pos_start = codon_finder(starts, seq['data'], inf=pos_stop[-1] + 3 + add_to_inf, sup=(length_seq - threshold))  # research of starts from 0 to (first stop - thresold)
        if pos_start:
            # looking for the appropriate stop:
            # It is here a bit tricky : The gene is at the border of the MetaG
            # We are looking for the position of the first nt
            # of the last full codon in the frame of the gene
            distance = length_seq - pos_start[-1]  # distance between the start and the end of the MetaG
            stop = length_seq - 3 - (distance % 3)  # -3 is -1 to get index of the last nt because here we start counting at 0 and -2 because we need the first nt of the last codon

            yield obj.Orf(frame=frame * rev, possible_start_orf=pos_start, end_orf=stop + 2, complet=False, border=True)

    return


def to_fit_len_max(twoStopLen):
    """
    Determine the borne inf for searching the start codons between two stop codons:
    Take the length between the two stop codon
    And give back the inferior borne where to begin the searc of start codon
    IF YOU WANT TO HAVE THE FULL ORF WITHOUT A MAX SIZE LIMITATION YOU CAN RETURN 0 ONLY
    """
    if twoStopLen <= obj.Gene.length_max:  # if the space between the two stop <= to len max then no need to add anything
        return 0
    # print 'twoStopLen {} - obj.Gene.length_max {}'.format(twoStopLen, obj.Gene.length_max)
    # print 'MODULO', (twoStopLen - obj.Gene.length_max) % 3
    return twoStopLen - obj.Gene.length_max


def complement_reverse(seq):
    compl_rev_seq = {}
    compl_rev_seq = reverse(seq)
    compl_rev_seq = complement(compl_rev_seq)
    # print 'In comp rev fct'
    compl_rev_seq['rev'] = not seq['rev']
    # print 'REV is ...', compl_rev_seq['rev']
    return compl_rev_seq


def reverse(seq):
    rev_seq = {}
    rev_seq['data'] = seq['data'][::-1]
    return rev_seq


def complement(seq):
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


def getGeneticCode(NCBI_ID):
    """Returns the genetic code corresponding to a given NCBI_ID

    This function is written by Theo Falgarone.

    Args:
        NCBI_ID: NCBI identifier for the genetic code

    Returns:
        table: genetic code table

    """
    if NCBI_ID == 1:
        base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
        start = "---M---------------M---------------M----------------------------"
        aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    elif NCBI_ID == 11:
        base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
        base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
        base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
        start = "---M---------------M------------MMMM---------------M------------"
        aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    else:
        print('ERROR : votre NCBI_ID ne correspond pas, veuillez entrer 1 ou 11')
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
    """Returns the translated proteic sequence based on the standard genetic code"""

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
    """Returns a list of non-overlapping words of a given length in a seq"""

    seq["data"] = seq["data"].upper()
    words = []
    for i in range(0, len(seq["data"]) - wlength + 1, wlength):
        words.append(seq["data"][i:i + wlength])

    return words
