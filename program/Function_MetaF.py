#!/usr/bin/env python3
import Object_MetaF as obj
import csv
import re
from subprocess import call
from operator import attrgetter
import os
import logging


def HMM_launcher(faa_file, add_to_name=''):
    hmm_db = obj.Gene.hmmdb
    table_hmm = obj.Gene.outdir + '/output_HMM_table_' + add_to_name + '.txt'

    bash_commande = "hmmsearch -E 0.5 --domtblout {} {} {} > /dev/null".format(
        table_hmm, hmm_db, faa_file)

    call(bash_commande, shell=True)

    if not os.path.isfile(table_hmm):
        logging.warning(f'bash command failed: {bash_commande}')
        raise FileNotFoundError(f'hmmsearch command failed to create the file {table_hmm}')

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
        print("seq['data'] is empty.. get_fast_fasta has failed to retreive the sequence {}".format(scaffold))


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

    # attribut of the classe gene like every objet (orf and TA_gene) will have this attribut, then no need to give it that each time
    obj.Gene.scaffold = scaffold
    domains_dict = {}  # keys : gene numbers, value: liste of the domain objet !!
    fl = open(table_hmm, 'r')
    for line in fl:
        # print l[:len(scaffold)]
        if line[:len(scaffold)] == scaffold:
            #  domain is an OBJECT of class Domain. It gather info about the domain found.
            nb_gene, domain = hmmtable_parser(line)

            domains_dict.setdefault(nb_gene, []).append(domain)
    fl.close()

    gff_fl = csv.reader(open(gff_file, 'r'), delimiter='\t')
    gff_line = ""  # next(gff_fl)
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
    for line in gff_handler:
        """
        it is a bit obscure but l is a list. To know the number of the gene we need to access to the last (8)
        tab delimiter area [8] and then split by ';' then select the first item [0]
        and split by | to finally get the gene number by selecting item 1.
        Example of gff line:
        ICM0007MP0313_1000001	GPMF	CDS	1	1101	.	+	0	ID=ICM0007MP0313_1000001|1;partial=10;s...
        CP000569.1	Genbank	CDS	2019629	2020228	.	-	0	ID=CP000569.1|1797;Parent=gene1855;Dbx.....
        """

        if line[2] not in ['CDS', 'ORF', 'gene'] or line[0] != scaffold:
            continue
        number = int(line[8].split(";")[0].split('|')[1])

        if number == gene_number:
            gene = obj.TA_gene()
            gene.feature = line[2]
            gene.start = int(line[3])
            gene.end = int(line[4])
            gene.strand = line[6]
            gene.gene_number = number

            # s = re.search(ur'protein_id=([^;\n]+)', l[8])
            s = re.search('protein_id=([^;\n]+)', line[8])
            if s:
                gene.protein_id = s.group(1)

            return gene, line, gff_handler


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
    if not result:
        raise Exception(
            'The parser of HMM table has failed... you may check the HMM table output ---> hmmtable line with a problem', line)
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

        # WARNING absolute value of distance min should not be greater than the length of the gene !!
        # because it would allow overlap of more than the length of the gene
        # and then could give a post gene before a pre gene so a non sense

        if abs(g.distanceMin) >= len(g):
            g.distanceMin = -len(g) + 1  # if so distance min is len -1

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
                # add the gene because it has a link in the set linked of class TA_gene
                obj.TA_gene.linked.add(gene)
                # add the gene because it has a link in the set linked of class TA_gene
                obj.TA_gene.linked.add(gpost)


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
                # WARNING absolute value of distance min should not be greater than the length of the gene !!
                # because it would allow overlap of more than the length of the gene
                # and then could give a post gene before a pre gene so a non sense
                if abs(gene.distanceMin) >= len(gene):
                    # print 'OOOOOOh'
                    # print gene.distanceMin
                    gene.distanceMin = -len(gene) + 1  # if so distance min is len -1
                    # print gene
                    # print gene.distanceMin
                    # raw_input()
            else:
                invalide.append(gene)
        for g in invalide:
            obj.TA_gene.genes.remove(g)
            obj.TA_gene.genes_strand[strand].remove(g)


def delete_files(listeFiles):
    from os import remove
    for outfile in listeFiles:
        try:
            remove(outfile)
        except OSError as e:  # if failed, report it back to the user ##
            print("Error: %s - %s." % (e.filename, e.strerror))


# def contig_stat_manager(writer_stat, scaffold, initial_nb_lonely, rescue, total_stat):
#     contig_stat = {}
#     contig_stat['gene with TA domain'] = len(obj.TA_gene.genes)
#     contig_stat['linked gene'] = len(obj.TA_gene.linked)
#     contig_stat['lonely gene'] = contig_stat['gene with TA domain'] - contig_stat['linked gene']
#     if rescue:
#         contig_stat['lonely gene rescue'] = initial_nb_lonely - contig_stat['lonely gene']
#         contig_stat['adjacent orf'] = obj.Orf.adj_orf_index
#         contig_stat['orf with TA domain'] = len(obj.Orf.hmm_orf)
#
#     for k in contig_stat: #add the value of the row in obj.Gene.metaG_stat to make the total at the end
#         total_stat[k] += contig_stat[k]  # contig is not there yet because it is not numerical value
#     total_stat['contig'] = scaffold
#     writer_stat.writerow(contig_stat)
#
#
# def write_result(set_linked, dict_output, scaffold):
#     i = 0
#     contig_header = "\n" + '==' * 2 + scaffold + '==' * 2 + '\n'
#     if dict_output['result_H']:
#         dict_output['result_H'].write(contig_header)
#     if dict_output['result_S']:
#         dict_output['result_S'].write(contig_header)
#     for gene in sorted(set_linked, key=attrgetter('start')):
#         for g_post in gene.post:
#             i += 1  # to give a number to each TA pair
#             g_score = gene.dict_score[g_post.gene_number]
#             post_score = g_post.dict_score[gene.gene_number]
#             if dict_output['result_H']:
#                 write_human_result(gene, g_post, dict_output['result_H'], i, g_score, post_score)
#             if dict_output['result_S']:
#                 write_short_result(gene, g_post, dict_output['result_S'], i, g_score, post_score)
#             if dict_output['result_T']:
#                 pass
#
#
# def write_short_result(g, post, fl, i, g_score, post_score):
#     if hasattr(g, 'protein_id'):
#         tag_g = g.protein_id
#     else:
#         tag_g = '_' + str(g.gene_number)
#
#     if hasattr(post, 'protein_id'):
#         tag_p = post.protein_id
#     else:
#         tag_p = '_' + str(post.gene_number)
#
#     fl.write("{}. Genes {} & {}\tstrand {}\tscore {}\n".format(
#         i, tag_g, tag_p, g.strand, g_score[0]['sum'] + post_score[0]['sum']))
#
#
# def write_human_result(g, post, fl, i, g_score, post_score):
#     fl.write("\nPRE GENE\n" + write_line(g, g_score))
#     fl.write("\nPOST GENE\n" + write_line(post, post_score) + '\n')
#     fl.write("DISTANCE {} ({})\n".format(post_score[0]['distance'], post_score[0]['dist_score']))
#     # fl.write(visualisation_genes(g, post, post_score[0]['distance']))
#
#
# def write_line(g, score):
#     line = "Gene {}\tfrom {} to {}\t{}aa ({})\tstart {}\t{}".format(
#         g.gene_number, g.real_start(), g.real_end(), score[0]["length"] / 3, score[0]["len_score"], score[0]["start"], g.feature)
#     domain_va = map(str, g.valid_domain(score[0]['start']))
#     domain_va = '/'.join(domain_va)
#     line += ("\tdomain: {}\n".format(domain_va))
#     return line
#
#
# def visualisation_genes(pre, post, distance):
#     pre_str = visual_str(len(pre))
#     post_str = visual_str(len(post))
#     dist_str = str(distance) + 'nt'
#     # print pre_str
#     # print post_str
#     # print 'disatance ', distance
#     sign = (distance + 1) / abs(distance + 1)
#     visual_dist = int(distance / 30.0 + 0.98 * sign)
#     # print visual_dist
#     final_str = pre_str + '\n'
#     position_g2 = (len(pre_str) + visual_dist)
#     final_str += position_g2 * ' ' + post_str + '\n'
#     positions = [position_g2, len(pre_str)]
#     final_str += (min(positions) - 1) * ' ' + '/' + abs(positions[0] - positions[1]) * ' ' + '\\\n'
#     final_str += (-1 + min(positions) + abs(positions[0] - positions[1]) - len(dist_str) / 2) * ' ' + dist_str
#     # print final_str
#     return final_str
#
#
# def visual_str(size):
#     # min and max size of the gene in characteres on the viual representation
#     # vlen = ((length - gene.length_min)) * ((visual_max - visual_min) / (gene.length_max - gene.length_min)) + visual_min
#     vlen = (size / 50) + 1
#
#     string = vlen * '=' + str(size / 3) + 'aa' + vlen * '=' + '>'
#     return string
