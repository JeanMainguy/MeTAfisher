# coding: utf-8
import csv
import Object_MetaF as obj
from operator import attrgetter, itemgetter


def score_manager(inf, sup, csv_file, k):
    dico = from_file_to_dict(csv_file, inf, sup)
    Ntot = sum(dico.values())
    assert k >= 0
    proba = give_proba_dict(dico, k, Ntot)
    # assert 1.01 > sum(proba.values()) > 0.99
    return proba  # proba_tranformed


def give_proba_dict(dico, k, N_tot):
    k = 30
    total = 0
    result = {}
    for nb in range(min(dico)-k, max(dico)+k + 1):
        somme = 0
        list_to_sum = [dico.get(i, 0) for i in range(nb - k, nb + k + 1)]
        mean_value = sum(list_to_sum)  # / float(len(list_to_sum))
        result[nb] = (mean_value / float(N_tot))
    return result


def from_file_to_dict(file_name, inf=None, sup=None):
    # anything will always be greater than None and lower than False
    # Not limit by default
    dico = {}

    with open(file_name) as csvfile:
        reader = csv.reader(csvfile)
        dico = {int(rows[0]): int(rows[1]) for rows in reader}
        if inf is None:
            inf = min(dico)
        if sup is None:
            sup = max(dico)
        dico = {k: v for k, v in dico.items() if inf <= k <= sup}
    return dico


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


def conflaction_proba(proba1, proba2):

    proba_final = (proba1 * proba2) / (proba1 * proba2 + (1-proba1)*(1-proba2))

    return proba_final


def confl(*proba):
    p_m = 1
    p_inv = 1
    for p in proba:
        assert 0 <= p <= 1, 'probability value is not between 0 and 1'
        p_m *= p
        p_inv *= (1 - p)

    conflation = p_m / (p_m + p_inv)
    return conflation


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

    score_post = get_score(post, compatible_starts, bonus_start, distance=initial_dist)
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

        dico_score = {"start": start, "domain": d.score, "len_score": proba_len,
                      'length': length, "score": proba_len}  # Ajout de distance:None? ? ?
        if distance is not None:
            # print('distance', distance + start)
            dico_score['dist_score'] = obj.Gene.distance_proba[distance + start]
            dico_score['distance'] = distance + start
            dico_score['score'] = conflaction_proba(dico_score['dist_score'], proba_len)

        # Give a bonus to the initial start !!
        dico_score['bonus_init_start'] = bonus_start if start == 0 else 0

        score.append(dico_score)
    score = sorted(score, key=itemgetter('score'), reverse=True)
    return score
