# coding: utf-8
import Object_MetaF as obj
import csv
# import Orf_MetaF as orf2
import re
from operator import attrgetter, itemgetter
from subprocess import call


def HMM_launcher(faa_file, add_to_name=''):
    hmm_db = obj.Gene.hmmdb
    table_hmm = obj.Gene.output_way + '/output_HMM_table_' + add_to_name + '.txt'

    bash_commande = "hmmsearch -E 0.5 --domtblout {} {} {} > /dev/null".format(table_hmm, hmm_db, faa_file)

    call(bash_commande, shell=True)

    return table_hmm


def get_fast_fasta(fl, line, scaffold):
    # line est la dernière ligne lu dans le fichier et doit correspondre à un entête >
    # la liste de scaffold est trié donc normalement y'a pas de problem.. Et les scaffold cherché devrait être dans le bonne ordre
    seq = {'description': scaffold, 'data': '', 'rev': False}
    # print scaffold + '=' + line[:len(scaffold) + 1]
    while line and line[:len(scaffold) + 1] != '>' + scaffold:
        line = next(fl)
    for line in fl:
        if line[0] == '>':
            break
        seq['data'] += line.rstrip()
    if seq['data']:
        seq['data'] = seq["data"].upper()
        return seq, line
    else:
        print ("seq['data'] is empty.. get_fast_fasta has failed to retreive the sequence {}".format(scaffold))


def get_hmm_genes(scaffold, table_hmm, gff_file):
    """
    Take all gene in table hmm and build an object off class ListGene() with them.
    scaffold = contig name
    table_hmm = file name
    gff_lines = list of gff line of the contig

    OUTPUT:
    genes : objet class ListGenes. Contain every gene of the given contig that are in table hmm.
    Info of the gff line are parsed and info about domain form hmm are also parsed.
    """

    obj.Gene.scaffold = scaffold  # attribut of the classe gene like every objet (orf and TA_gene) will have this attribut, then no need to give it that each time
    domains_dict = {}  # keys : gene numbers, value: liste of the domain objet !!
    fl = open(table_hmm, 'r')
    for l in fl:
        # print l[:len(scaffold)]
        if l[:len(scaffold)] == scaffold:
            #  domain is an OBJECT of class Domain. It gather info about the domain found.
            nb_gene, domain = hmmtable_parser(l)

            domains_dict.setdefault(nb_gene, []).append(domain)
    fl.close()

    gff_fl = csv.reader(open(gff_file, 'r'), delimiter='\t')
    gff_line = next(gff_fl)
    for n in sorted(domains_dict):
        """
        Search the gene info of the gene in the gff file
        n is the keys of domains dict and then take the value all the gene number with a hmm hit
        """
        gene, gff_line, gff_fl = get_gff_info(gff_fl, gff_line, scaffold, n)

        gene.domain = domains_dict[n]

        gene.domain_Ct_border = max((d.ali_from * 3 for d in domains_dict[n]))
        obj.TA_gene.genes.append(gene)
        # if gene.strand == '+':
        #     obj.TA_gene.genes_plus.append(gene)
        # else:
        #     obj.TA_gene.genes_minus.append(gene)
        obj.TA_gene.genes_strand[gene.strand].append(gene)
    # return obj.TA_gene.genes


def get_gff_info(gff_handler, gff_line, scaffold, gene_number):
    """
    gene_number need to come sorted
    return: the file handler of the gff file
    """
    for l in gff_handler:
        """
        it is a bit obscure but l is a list. To know the number of the gene we need to access to the last (8)
        tab delimiter area [8] and then split by ';' then select the first item [0]
        and split by | to finally get the gene number by selecting item 1.
        Example of gff line:
        ICM0007MP0313_1000001	GPMF	CDS	1	1101	.	+	0	ID=ICM0007MP0313_1000001|1;partial=10;s...
        CP000569.1	Genbank	CDS	2019629	2020228	.	-	0	ID=CP000569.1|1797;Parent=gene1855;Dbx.....
        """

        if l[2] not in ['CDS', 'ORF', 'gene'] or l[0] != scaffold:
            continue
        number = int(l[8].split(";")[0].split('|')[1])

        if number == gene_number:
            gene = obj.TA_gene()
            gene.feature = l[2]
            gene.start = int(l[3])
            gene.end = int(l[4])
            gene.strand = l[6]
            gene.gene_number = number

            s = re.search(ur'locus_tag=([^;\n]+)', l[8])
            if s:
                gene.locus_tag = s.group(1)

            return gene, l, gff_handler


def hmmtable_parser(line):
    """
    Regex : parse a line of the hmmtable and return a object from Domain class with all the info of the line
    """
    pattern = re.compile(r"""
    (?P<scaffold>[^|]+)\|(?P<gene_number>\d{1,4}) #scaffold name and gene_number separated by a |
    \s+-\s+
    (?P<gene_len>\d+) #length of the gene sequence in residu
    \s+
    (?P<domain>[\w.\(\)_+-]+) #domain name
    \s+
    (?P<domain_acc>[\w.+-]+) #domain acc pour pfam
    \s+\d+\s+
    (?P<evalue>[\w.-]+) # E-value of the overall sequence/profile comparison (including all domains).
    \s+
    (?P<score>[\d.]+) # Bit score of the overall sequence/profile comparison (
    \s+[\d.]+\s+
    (?P<domain_number>\d+) # This domain’s number
    \s+
    (?P<total_domain>\d+) # The total number of domains reported in the sequence, ndom.
    \s+[\w.-]+\s+[\w.+-]+\s+[\w.-]+\s+[\w.-]+\s+\d+\s+\d+\s+
    (?P<ali_from>\d+) # The start of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
    \s+
    (?P<ali_to>\d+) # The end of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
    \s+
    (?P<env_from>\d+) # The start of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
                      # The envelope defines a subsequence for which their is substantial probability mass supporting a homologous domain,
                      # whether or not a single discrete alignment can be identified.
                      # The envelope may extend beyond the endpoints of the MEA alignment, and in fact often does, for weakly scoring domains.
    \s+
    (?P<env_to>\d+) # The end of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
    \s
    """, re.VERBOSE)
    """
    regex in one line:
    (?P<scaffold>[^|]+)\|(?P<gene_number>\d{1,4})\s+-\s+(?P<gene_len>\d+)\s+(?P<domain>[\w.-]+)\s+(?P<domain_acc>[\w.+-]+)\s+\d+\s+(?P<evalue>[\w.-]+)\s+(?P<score>[\d.]+)\s+[\d.]+\s+(?P<domain_number>\d+)\s+(?P<total_domain>\d+)\s+[\w.+-]+\s+([\w.+-]+)\s+[\w.+-]+\s+[\w.+-]+\s+\d+\s+\d+\s+(?P<ali_from>\d+)\s+(?P<ali_to>\d+)\s+(?P<env_from>\d+)\s+(?P<env_to>\d+)\s

    link : https://regex101.com/r/kKHzsB/2
    """

    result = pattern.match(line)

    domain = obj.Domain(
        domain_number=result.group("domain_number"), ali_from=result.group("ali_from"), ali_to=result.group("ali_to"), env_from=result.group("env_from"),
        env_to=result.group("env_to"), domain_name=result.group("domain"), domain_acc=result.group("domain_acc"), e_value=result.group("evalue"), score=result.group("score"), line=line)
    gene_number = int(result.group("gene_number"))
    return gene_number, domain


def get_start_po(dico):
    """
    Process the search of start position over the main liste of gene obj.TA_gene.genes !

    """
    obj.TA_gene.genes = sorted(obj.TA_gene.genes, key=attrgetter('gene_number'))
    non_valides = []
    for g in obj.TA_gene.genes:
        start_po = g.give_start_po(dico, min_max_intervall=True)
        # print 'gene number', g.gene_number
        # print 'start po : ', start_po
        # print 'START PO \t', start_po
        # print 'START PO', start_po
        # g.possible_starts = start_po # nouvelle attribut tous les position des start possibles

        if not start_po:
            # print 'invalide', g
            non_valides.append(g)
            continue

        g.distanceMin = obj.Gene.distanceMin - abs(start_po[-1] - start_po[0])

        # Uncomment this part if instead of having the position of start inside the intial gene,
        # it is desire to have the position inthe contig

        if g.strand == '+':
            # start_po = [x + g.start for x in start_po]
            g.start += start_po[0]

        else:
            # start_po = [g.end - x for x in start_po]
            g.end -= start_po[0]

        g.possible_start = start_po  # nouvelle attribut tous les position des start possibles
        # print 'POSSIBLE START', g.possible_start
    # ATTENTION need to discard the gene in the main liste and also in minus and plus
    for g in non_valides:
        # print 'INVALIDE',g.gene_number
        obj.TA_gene.genes.remove(g)
        obj.TA_gene.genes_strand[g.strand].remove(g)
        # if g.strand == '+':
        #     obj.TA_gene.genes_plus.remove(g)
        # else:
        #     obj.TA_gene.genes_minus.remove(g)


def get_list_scaffold(table_hmm):
    # input : Name of the table hmm file
    # output : List of all the scaffold of the table

    set_of_scaffold = set()
    fl_hmm = open(table_hmm, 'r')
    for line in fl_hmm:
        if line[0] == "#":
            continue
        set_of_scaffold.add(line[:line.index('|')])
    fl_hmm.close()
    return list(set_of_scaffold)


def get_adj():
    """
    Treat independently list plus and minus
    Sort the list by the end position
    then take the first one and look for the next gene (the gene has the its end after the end of the first one).
    call tandem gene function to check if the gene fit the distance threshold
    """
    # obj.TA_gene.genes_plus = sorted(obj.TA_gene.genes_plus, key=attrgetter('end'))
    # obj.TA_gene.genes_minus = sorted(obj.TA_gene.genes_minus, key=attrgetter('start'))
    obj.TA_gene.genes_strand['+'] = sorted(obj.TA_gene.genes_strand['+'], key=attrgetter('end'))
    obj.TA_gene.genes_strand['-'] = sorted(obj.TA_gene.genes_strand['-'], key=attrgetter('start'))
    #
    # TODO Ameliorer mettre un break quand les gene sont de toute facon trop loin.. bu actually harder than excepted because
    # Passe si le end du gene suivant est plus loin que la taille max + la
    # distance max dans ces cas là ya pas de possbilité pour qu'il soit adj!

    # Strand plus !
    # adj_by_strand(obj.TA_gene.genes_plus)
    # adj_by_strand(obj.TA_gene.genes_minus)
    adj_by_strand(obj.TA_gene.genes_strand['+'])
    adj_by_strand(obj.TA_gene.genes_strand['-'])


def adj_by_strand(liste):
    """
    liste: list of hmm gene with homogenous strand
    Check if the gene is in tandem with another and if so store the gene inside a set obj.TA_gene.linked
    In parallel it clean up the list obj.TA_gene.genes
    by removing the genes that forme a tandem. Then TA_gene.genes has only the lonely_gene
    """
    for gi, gene in enumerate(liste):
        # print obj.TA_gene.genes_plus[gi].gene_number,  obj.TA_gene.genes_plus[gi].len_val
        for gpost in liste[gi + 1:]:
            if gpost.end - gene.end + 1 > obj.Gene.length_max + obj.Gene.distanceMax:
                """
                if the distance between gene.end and gpost.end is superior to lenmax + distmax
                Then the two gene won't be in pair and the next postgene either because they are sorted by their start
                So we can break the gpost for loop and check the next gene
                """
                # print 'distance between gene number {} .end and gpost {}.end is superior to lenmax + distmax'.format(gene.gene_number, gpost.gene_number)
                break
            # it is a simple test that ckeck if the two gene are adjacent
            if gene.is_pre_adj_to(gpost):
                #  store the information of prev and post according the strand
                if gene.strand == '+':
                    gene.post.append(gpost)
                    gpost.prev.append(gene)
                else:
                    gpost.post.append(gene)
                    gene.prev.append(gpost)
                obj.TA_gene.linked.add(gene)  # add the gene because it has a link in the set linked of class TA_gene
                obj.TA_gene.linked.add(gpost)  # add the gene because it has a link in the set linked of class TA_gene


def only_lonely_gene():
    """
    Sine Round1 and Round2 are merged in a single run this fct is hasbeen
    Function only useful in round1 when we focus only on rescuing lonely gene
    Make space by deleting all the obj that are grouped
    Normally at the end genes_minus and genes_plus lists contain only lonely gene.

    """
    for link_gene in obj.TA_gene.linked:
        obj.TA_gene.genes.remove(link_gene)
        try:
            obj.TA_gene.genes_strand['+'].remove(link_gene)
        except ValueError:
            obj.TA_gene.genes_strand['-'].remove(link_gene)
    obj.TA_gene.linked.clear()


def create_lonely_gene_list():
    """
    From Gene.genes_minus and plus liste we build a lonely_plus and minus liste that gather lonely genes
    so gene that are not in the set linked
    If there is no lonely gene then the attribute lonely stays None
    """
    if len(obj.TA_gene.genes) > len(obj.TA_gene.linked):  # that mean there are lonely gene
        obj.TA_gene.lonely = {'+': [], '-': []}
        for gene in obj.TA_gene.genes:
            if gene not in obj.TA_gene.linked:
                obj.TA_gene.lonely[gene.strand].append(gene)


# def rescue_lonely_gene(dico_orf, dico_gff, scaffold, tmp_adjorf_faa):
#     fl = open(tmp_adjorf_faa, 'w')
#     # Extraction of the sequence of the contig in order to find the orf inside.
#     seq, dico_orf["line"] = get_fast_fasta(dico_orf['fl'], dico_orf['line'], scaffold)
#     obj.Orf.seq = seq
#     # retrieve the end position of the predicted gene to then skip the related ORF:
#     gff_ends = orf2.get_gff_ends(dico_gff, scaffold)
#
#     # TEST IF THERE ARE LONELY GENE IN STRAND plu
#     if obj.TA_gene.lonely['+']:
#         print 'THERE IS LONELY GENE ON STRAND +'
#         generator = orf2.findORF(scaffold, seq, rev=1)
#         # read the generator: check if the orf is a predicted gene
#         # then check if the orf is adjacent with a TA_gene
#         # If so then the faa fa gff ORF file are written...
#         orf2.orf_manager(generator, '+', obj.TA_gene.genes_strand['+'], gff_ends, fl)  # orf manager will do everything !!
#     # else:
#         # raw_input(scaffold + " THERE IS NO LONELY GENE ON STRAND +")
#     # TEST IF THERE ARE LONELY GENES IN STRAND minus
#     if obj.TA_gene.lonely['-']:
#         print 'THERE IS LONELY GENE ON STRAND -'
#         generator = orf2.findORF(scaffold, seq, rev=-1)
#         orf2.orf_manager(generator, '-', obj.TA_gene.genes_strand['-'], gff_ends, fl)
#     # else:
#         # raw_input(scaffold + " THERE IS NO LONELY GENE ON STRAND -")
#     fl.close()
#     # print obj.Orf.adj_orf
#     if obj.Orf.adj_orf:
#         orf2.adjOrf_HMM()
#         adjust_possible_start_to_domain(obj.Orf.hmm_orf)
#         hmm_orf_get_adj(obj.Orf.hmm_orf)
#     # for k in obj.Orf.hmm_orf:
#     #     for o in obj.Orf.hmm_orf[k]:
#     #         print o
#     #         print o.distanceMin
#     #         print 'gene number', o.gene_number
#     #         print 'prev', [g.gene_number for g in o.prev]
#     #         print 'post', [g.gene_number for g in o.post]
#     #         print 'Domain border', o.domain_Ct_border
#     #         print o.possible_start


# def adjust_possible_start_to_domain(dico_obj):
#     for strand in dico_obj:
#         for o in dico_obj[strand]:
#             # print 'CT BORDER', o.domain_Ct_border
#             # print 'before', o.possible_start
#             # starts_zip = zip(o.possible_start, o.possible_start_orf)
#             o.possible_start = [s for s in o.possible_start if s <= o.domain_Ct_border]
#             # print 'after', o.possible_start
#             o.distanceMin = obj.Gene.distanceMin - abs(o.possible_start[-1] - o.possible_start[0])
#             o.post = []
#             o.prev = []

def check_size(gene_strand):
    """
    Remove genes with a length that does not fit the length thresholds.
    """
    for strand in gene_strand:
        invalide = []
        for gene in gene_strand[strand]:
            # print gene.gene_number, '.....*'
            if obj.Gene.length_min <= len(gene) <= obj.Gene.length_max:
                gene.possible_start = [0]
            else:
                invalide.append(gene)
        for g in invalide:
            obj.TA_gene.genes.remove(g)
            obj.TA_gene.genes_strand[strand].remove(g)


def transformation(dico, multiple):
    maximum = max(dico.itervalues())
    minimum = min(dico.itervalues())
    for k in dico:
        dico[k] = ((dico[k] - minimum) / (maximum - minimum)) * multiple
    return dico


def give_proba_dict(inf, sup, dico, k, N_tot):
    result = {}
    for nb in xrange(inf, sup + 1):
        somme = 0
        for i in xrange(nb - k, nb + k + 1):
            # print i
            somme += dico.get(nb + i, 0)
        result[nb] = (somme / float(N_tot) / (2 * k + 1))

    return result


# def hmm_orf_get_adj(dico_hmmorf):
#     # Add the hmm orf into the main lists !
#     # obj.TA_gene.genes.extend(dico_hmmorf['+'] + dico_hmmorf['-'])
#     # obj.TA_gene.genes_plus.extend(dico_obj['+'])
#     # obj.TA_gene.genes_minus.extend(dico_obj['-'])
#     for strand in dico_hmmorf:
#         obj.TA_gene.genes_strand[strand].extend(dico_hmmorf[strand])
#         obj.TA_gene.genes.extend(dico_hmmorf[strand])
#         for hmmorf in dico_hmmorf[strand]:
#             for gene in obj.TA_gene.genes_strand[strand]:
#                 if hmmorf.gene_number == gene.gene_number:  # MAGIC METHOD
#                     # print hmmorf
#                     # print gene
#                     continue
#                 if gene.is_pre_adj_to(hmmorf):
#                     obj.TA_gene.linked.add(gene)
#                     obj.TA_gene.linked.add(hmmorf)
#                     if strand == '+':
#                         gene.post.append(hmmorf)
#                         hmmorf.prev.append(gene)
#                     else:
#                         hmmorf.post.append(gene)
#                         gene.prev.append(hmmorf)
#                 elif hmmorf.is_pre_adj_to(gene):
#                     obj.TA_gene.linked.add(gene)
#                     obj.TA_gene.linked.add(hmmorf)
#                     if strand == '+':
#                         hmmorf.post.append(gene)
#                         gene.prev.append(hmmorf)
#                     else:
#                         gene.post.append(hmmorf)
#                         hmmorf.prev.append(gene)


def score_TA_list(genes_strand):
    for strand in genes_strand:
        for gene in genes_strand[strand]:
            for post in gene.post:
                # print 'Gene :\n', gene
                # print 'Gene Post :\n', post
                score_pair(gene, post)  # GIVE THE SCOREEE
                # fct2.write_human_result(g, g_post, fl_hu_res)


def score_pair(pre, post):  # post is a gene located upstream of pre !
    """
    Give the score of the pair
    return False if not possible to pair them
    TODO : le post et prev gene are post and prev according their position in the system and not in the contig
    strand + : ---===self===>---==post==>--
    strand - : ---<===post==---<===self==---
    """

    # get the distancebetween the two gene
    initial_dist = post.real_start() - pre.real_end()
    if pre.strand == "-":
        initial_dist *= -1
    compatible_starts = [s for s in post.possible_start if
                         obj.Gene.distanceMin < initial_dist + s < obj.Gene.distanceMax]

    # print "MIN {} --> MAX {}".format(obj.Gene.distanceMin, obj.Gene.distanceMax)
    # print "possible start of g"
    # print ('initial distance and strand ', initial_dist, pre.strand)
    # print dir(post)
    # print ('post g possible start', post.possible_start)
    # print ('compatible start', compatible_starts)
    # print ('possible start of self ! ', pre.possible_start)
    score_post = get_score(post, starts=compatible_starts, distance=initial_dist)
    score = get_score(pre, pre.possible_start)
    # print pre.feature
    # print pre.dict_score
    pre.dict_score[post.gene_number] = score
    post.dict_score[pre.gene_number] = score_post
    # for s in score:
    #     print s


def get_score(gene, starts, distance=None):
    """
    Return a double dico with first start as a key and another dico
    with score of domain and length and when pre_end isn't Noneth distance score
    """
    score = []
    gene.domain = sorted(gene.domain, key=attrgetter('score'), reverse=True)
    domains = iter(gene.domain)
    d = domains.next()
    # print 'starts process by get_score()', starts
    for start in starts:
        while d.ali_from * 3 < start:
            try:
                d = domains.next()
            except StopIteration:
                break

        length = len(gene) - start
        # print '===GENNNNNNE==='
        # print gene
        # print '===LEN==='
        # print 'gene.end {} - gene.start {} + start {} + 1  GIVE THE LEN of {} nt, {} aa '.format(gene.end, gene.start, start, length, length / 3)

        proba_len = obj.Gene.length_proba[length / 3]

        # print 'proba len ', proba_len
        dico_score = {"start": start, "domain": d.score_transformed(), "len_score": proba_len, 'length': length, "sum": d.score_transformed() + proba_len}  # Ajout de distance:None? ? ?
        if distance is not None:
            # print('distance', distance + start)
            dico_score['dist_score'] = obj.Gene.distance_proba[distance + start]
            dico_score['distance'] = distance + start
            dico_score['sum'] += dico_score['dist_score']
        score.append(dico_score)
    score = sorted(score, key=itemgetter('sum'), reverse=True)
    return score


def from_file_to_dict(file_name):
    dico = {}
    with open(file_name) as csvfile:
        reader = csv.reader(csvfile)
        N_tot = 0
        for rows in reader:
            dico[int(rows[0])] = int(rows[1])
            N_tot += int(rows[1])
    return dico, N_tot


def delete_files(listeFiles):
    from os import remove
    for outfile in listeFiles:
        try:
            remove(outfile)
        except OSError, e:  # if failed, report it back to the user ##
            print ("Error: %s - %s." % (e.filename, e.strerror))


def contig_stat_manager(writer_stat, scaffold, initial_nb_lonely, rescue):
    contig_stat = {}
    contig_stat['gene with TA domain'] = len(obj.TA_gene.genes)
    contig_stat['linked gene'] = len(obj.TA_gene.linked)
    contig_stat['lonely gene'] = contig_stat['gene with TA domain'] - contig_stat['linked gene']
    if rescue:
        contig_stat['lonely gene rescue'] = initial_nb_lonely - contig_stat['lonely gene']
        contig_stat['adjacent orf'] = obj.Orf.adj_orf_index
        contig_stat['orf with TA domain'] = len(obj.Orf.hmm_orf)

    for k in contig_stat: #add the value of the row in obj.Gene.metaG_stat to make the total at the end
        obj.Gene.metaG_stat[k] += contig_stat[k]  # contig is not there yet because it is not numerical value
    contig_stat['contig'] = scaffold
    writer_stat.writerow(contig_stat)


def write_result(set_linked, dict_output, scaffold):
    i = 0
    contig_header = "\n" + '==' * 2 + scaffold + '==' * 2 + '\n'
    if dict_output['result_H']:
        dict_output['result_H'].write(contig_header)
    if dict_output['result_S']:
        dict_output['result_S'].write(contig_header)
    for gene in sorted(set_linked, key=attrgetter('start')):
        for g_post in gene.post:
            i += 1  # to give a number to each TA pair
            g_score = gene.dict_score[g_post.gene_number]
            post_score = g_post.dict_score[gene.gene_number]
            if dict_output['result_H']:
                write_human_result(gene, g_post, dict_output['result_H'], i, g_score, post_score)
            if dict_output['result_S']:
                write_short_result(gene, g_post, dict_output['result_S'], i, g_score, post_score)
            if dict_output['result_T']:
                pass


def write_short_result(g, post, fl, i, g_score, post_score):
    if hasattr(g, 'locus_tag'):
        tag_g = g.locus_tag
    else:
        tag_g = '_' + str(g.gene_number)

    if hasattr(post, 'locus_tag'):
        tag_p = post.locus_tag
    else:
        tag_p = '_' + str(post.gene_number)

    fl.write("{}. Genes {} & {}\tstrand {}\tscore {}\n".format(
        i, tag_g, tag_p, g.strand, g_score[0]['sum'] + post_score[0]['sum']))


def write_human_result(g, post, fl, i, g_score, post_score):
    fl.write("\nPRE GENE\n" + write_line(g, g_score))
    fl.write("\nPOST GENE\n" + write_line(post, post_score) + '\n')
    fl.write("DISTANCE {} ({})\n".format(post_score[0]['distance'], post_score[0]['dist_score']))
    # fl.write(visualisation_genes(g, post, post_score[0]['distance']))


def write_line(g, score):
    line = "Gene {}\tfrom {} to {}\t{}aa ({})\tstart {}\t{}".format(
        g.gene_number, g.real_start(), g.real_end(), score[0]["length"] / 3, score[0]["len_score"], score[0]["start"], g.feature)
    domain_va = map(str, g.valid_domain(score[0]['start']))
    domain_va = '/'.join(domain_va)
    line += ("\tdomain: {}\n".format(domain_va))
    return line


def visualisation_genes(pre, post, distance):
    pre_str = visual_str(len(pre))
    post_str = visual_str(len(post))
    dist_str = str(distance) + 'nt'
    # print pre_str
    # print post_str
    # print 'disatance ', distance
    sign = (distance + 1) / abs(distance + 1)
    visual_dist = int(distance / 30.0 + 0.98 * sign)
    # print visual_dist
    final_str = pre_str + '\n'
    position_g2 = (len(pre_str) + visual_dist)
    final_str += position_g2 * ' ' + post_str + '\n'
    positions = [position_g2, len(pre_str)]
    final_str += (min(positions) - 1) * ' ' + '/' + abs(positions[0] - positions[1]) * ' ' + '\\\n'
    final_str += (-1 + min(positions) + abs(positions[0] - positions[1]) - len(dist_str) / 2) * ' ' + dist_str
    # print final_str
    return final_str


def visual_str(size):
    # min and max size of the gene in characteres on the viual representation
    # vlen = ((length - gene.length_min)) * ((visual_max - visual_min) / (gene.length_max - gene.length_min)) + visual_min
    vlen = (size / 50) + 1

    string = vlen * '=' + str(size / 3) + 'aa' + vlen * '=' + '>'
    return string
