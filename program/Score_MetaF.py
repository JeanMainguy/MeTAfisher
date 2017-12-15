# coding: utf-8
import csv

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
