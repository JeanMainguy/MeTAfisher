# coding: utf-8
import csv
import Object_MetaF as obj
from operator import attrgetter, itemgetter

def score_manager(mini, maxi, csv_file, k, mltp):

    dico, Ntot = from_file_to_dict(csv_file)

    proba = give_proba_dict(mini, maxi, dico, k, Ntot)
    proba_tranformed = transformation(proba, mltp)

    return proba_tranformed

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

def from_file_to_dict(file_name):
    dico = {}
    with open(file_name) as csvfile:
        reader = csv.reader(csvfile)
        N_tot = 0
        for rows in reader:
            dico[int(rows[0])] = int(rows[1])
            N_tot += int(rows[1])
    return dico, N_tot


def score_TA_list(genes_strand, bonus_start):
    # for strand in genes_strand:
    #     for gene in genes_strand[strand]:
    #         for post in gene.post:
    #             # print 'Gene :\n', gene
    #             # print 'Gene Post :\n', post
    #             score_pair(gene, post, bonus_start)  # GIVE THE SCOREEE
    #             # fct2.write_human_result(g, g_post, fl_hu_res)
    for gene in genes_strand:
        for post in gene.post:
            score_pair(gene, post, bonus_start)  # GIVE THE SCOREEE

    # for gene in genes_strand:
    #     print '=='*20
    #     print 'GENE', gene
    #     print gene.dict_score
    #     print '=='*20


def score_pair(pre, post, bonus_start):  # post is a gene located upstream of pre !
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
    # # print dir(post)
    # print ('post g possible start', post.possible_start)
    # print ('compatible start', compatible_starts)
    # print ('possible start of self ! ', pre.possible_start)
    score_post = get_score(post, compatible_starts, bonus_start, distance=initial_dist)
    # print 'score post de ',post.gene_number,   score_post
    score = get_score(pre, pre.possible_start, bonus_start)
    # print pre.feature
    # print pre.dict_score
    pre.dict_score[post.gene_number] = score
    post.dict_score[pre.gene_number] = score_post
    # for s in score:
    #     print s


def get_score(gene, starts, bonus_start, distance=None):
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
        # Give a bonus to the initial start !!
        dico_score['bonus_init_start'] = bonus_start if start == 0 else 0
        dico_score['sum'] += dico_score["bonus_init_start"]

        score.append(dico_score)
    score = sorted(score, key=itemgetter('sum'), reverse=True)
    return score
