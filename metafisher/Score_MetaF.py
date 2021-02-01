# coding: utf-8
import csv
import Object_MetaF as obj
from operator import attrgetter, itemgetter
import json


def decoder(filename):
    with open(filename, 'r') as file:
        return json.load(file)


def score_manager(inf, sup, csv_file, k):
    dico = from_file_to_dict(csv_file, inf, sup)
    Ntot = sum(dico.values())
    assert k >= 0
    proba = give_proba_dict(dico, k, Ntot)
    # assert 1.01 > sum(proba.values()) > 0.99
    return proba  # proba_tranformed


def give_proba_dict(dico, k, N_tot):
    k = 30
    result = {}
    for nb in range(min(dico)-k, max(dico)+k + 1):
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


def score_TA_list(genes, score_dict):
    for gene in genes:
        for post in gene.post:
            score_pair(gene, post, score_dict)  # GIVE THE SCOREEE


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


def domain_association_prct(do1, do2, association, domain_types):
    try:
        number_of_association = float(association[do1][do2])
    except KeyError:
        number_of_association = 0.0
    nb_do1_found = domain_types[do1]['AT'] + domain_types[do1]['T']
    nb_do2_found = domain_types[do2]['AT'] + domain_types[do2]['T']

    if nb_do1_found < nb_do2_found:
        perct_association = number_of_association / nb_do1_found
    else:
        perct_association = number_of_association / nb_do2_found

    return perct_association


def score_do_association(post, pre, dict_domain_association, dict_domain_gene_type):

    asso = []
    for pre_do in pre.domain:

        for post_do in post.domain:
            perct_asso = domain_association_prct(
                pre_do.domain_name, post_do.domain_name, dict_domain_association, dict_domain_gene_type)
            asso.append(perct_asso)
    return max(asso)


def score_pair(pre, post, score_dict):  # post is a gene located upstream of pre !
    """
    Give the score of the pair
    return False if not possible to pair them
    TODO : le post et prev gene are post and prev according their position in the system and not in the contig
    strand + : ---===self===>---==post==>--
    strand - : ---<===post==---<===self==---
    """

    # get the distance between the two genes
    initial_dist = post.real_start() - pre.real_end()
    if pre.strand == "-":
        initial_dist *= -1
    # possible_starts attribute is 0 when we have TA gene and resize flag is false
    #
    compatible_starts = [s for s in post.possible_starts if
                         obj.Gene.distanceMin < initial_dist + s < obj.Gene.distanceMax]

    score_post = get_score(post, compatible_starts, distance=initial_dist,
                           length_proba=score_dict['length_proba'], distance_proba=score_dict['distance_proba'])
    score = get_score(pre, pre.possible_starts,
                      length_proba=score_dict['length_proba'], distance_proba=score_dict['distance_proba'])

    score_do_asso = score_do_association(post, pre, score_dict['domain_association'],
                                         score_dict['domain_gene_type'])

    for s in score_post:  # don't manage alternative start.. !!
        s['domain_association_score'] = score_do_asso

    pre.dict_score[post.gene_number] = score
    post.dict_score[pre.gene_number] = score_post


def get_score(gene, starts, length_proba, distance_proba, distance=None):
    """
    Return a double dict with first start as a key and another dico
    with score of domain and length and when pre_end isn't Noneth distance score
    """
    scores = []
    gene.domain = sorted(gene.domain, key=attrgetter('score'), reverse=True)
    # print([d.score for d in gene.domain])
    domains = iter(gene.domain)
    d = next(domains)

    for start in starts:
        while d.ali_from * 3 < start:  # to not take into account domain that would be before of the current start investiagted
            try:
                d = next(domains)
            except StopIteration:
                break

        length = len(gene) - start
        assert length % 3 == 0

        proba_len = length_proba[int(length / 3)]

        dico_score = {"start": start, "domain": d.score, "len_score": proba_len,
                      'length': length, "score": proba_len}
        if distance is not None:
            dico_score['dist_score'] = distance_proba[distance + start]
            dico_score['distance'] = distance + start
            dico_score['score'] = conflaction_proba(dico_score['dist_score'], proba_len)

        scores.append(dico_score)
    scores = sorted(scores, key=itemgetter('score'), reverse=True)
    return scores
