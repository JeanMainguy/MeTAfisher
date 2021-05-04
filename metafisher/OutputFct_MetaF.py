#!/usr/bin/env python3

"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 27 nov. 2020
License     : MIT
Maintainer  : jean.mainguy@outlook.fr


Program to retrieve toxin antitoxin (TA) systems in genomes or metagenomes.
"""

import csv
import Function_MetaF as fct
from operator import attrgetter
import Score_MetaF as score
import os


def output_manager(output_way, metaG_name, thresholds, dict_output, info_contig_stat, rescue=False, resize=False):
    headinfo, complement = output_headinfo_creation(metaG_name, thresholds, rescue, resize)
    # Output H, T and S initialization
    # dico output is update no need to give it back
    output_file_creation(output_way, metaG_name, dict_output, headinfo, complement)
    # Output T as a csv file initialization
    dict_output["result_TA_genes"] = result_ta_genes_config(dict_output["result_TA_genes"])
    dict_output["result_TA_pairs"] = result_ta_pairs_config(dict_output["result_TA_pairs"])
    # Stat file initialisation
    writer_stat, total_stat, fl_stat = stat_file_creation(
        output_way, metaG_name, info_contig_stat, headinfo, complement, rescue)
    return writer_stat, total_stat, fl_stat


def output_headinfo_creation(metaG_name, thresholds, rescue, resize):
    complement = ''
    if rescue:
        complement += '_rescue'
    if resize:
        complement += '_resize'

    headinfo = f'## Name of the sequence analysed: {metaG_name}\n'
    headinfo += f"\n## Rescue lonely gene : {rescue}\n"
    headinfo += f"## Resize gene : {resize}\n"
    headinfo += f"## Distance threshold from {thresholds['distanceMin']}nt to {thresholds['distanceMax']}nt\n"
    headinfo += f"## Length threshold from {int(thresholds['lenMin']/3)}aa to {int(thresholds['lenMax']/3)}aa\n"
    return headinfo, complement


def result_ta_genes_config(fl):
    if not fl:
        return False  # fl is False
    # Chnage the regular file writer into a csv file writer with header
    header = ["contig", "gene_id", "start", "end",
              "length", "strand", "feature", "possible_partners", 'TA_domains', 'TA_families', 'TADB_hits']
    writer_table = csv.DictWriter(fl, fieldnames=header, delimiter='\t')
    writer_table.writeheader()
    return writer_table


def result_ta_pairs_config(fl):
    if not fl:
        return False  # fl is False
    # Chnage the regular file writer into a csv file writer with header
    header = ["contig",  "strand",
              "gene1_id", "gene1_start", "gene1_end", "gene1_length", "gene1_length_score",
              "gene2_id", "gene2_start", "gene2_end", "gene2_length", "gene2_length_score",
              "distance", "distance_score", "TA_association_score", "system_score",
              'gene1_TA_domains', "gene1_TA_families", 'gene1_TADB_hits', 'gene2_TA_domains', "gene2_TA_families", 'gene2_TADB_hits', 'shared_family']
    writer_table = csv.DictWriter(fl, fieldnames=header, delimiter='\t')
    writer_table.writeheader()
    return writer_table


def output_file_creation(output_way, metaG_name, dict_output, headinfo, complement):
    # the fl of each kind of output are stored in a dictionnary:
    # key : name of the output | value : fl or False if not wanted
    is_output = False  # by default it is False and then if one of the output is True, it will become True

    for out_name in dict_output:
        if out_name:  # if the flag is not False
            extension = 'txt'
            if out_name.startswith('result_TA_'):
                extension = 'tsv'
            elif out_name == 'result_GFF':
                extension = 'gff'
                # out_name = 'TA_Genes'

            file_out = os.path.join(output_way, f"{metaG_name}_{out_name}{complement}.{extension}")
            flout = open(file_out, "w")
            # flout.write("## {}\n".format(out_name))
            # flout.write(headinfo)
            dict_output[out_name] = flout
            is_output = True
    dict_output['is_output'] = is_output


def stat_file_creation(output_way, metaG_name, info_contig_stat, headinfo, complement, rescue):

    if info_contig_stat:
        header = ['contig', 'gene with TA domain', 'lonely gene', 'linked gene']

        if rescue:
            header += ['adjacent orf', 'orf with TA domain', 'lonely gene rescue']

        fl_stat = open('{}/{}_contig_stat{}.tsv'.format(output_way, metaG_name, complement), 'w')
        fl_stat.write(headinfo)
        writer_stat = csv.DictWriter(fl_stat, fieldnames=header, delimiter='\t')
        writer_stat.writeheader()
        total_stat = dict.fromkeys(header, 0)
        total_stat['contig'] = metaG_name
        return writer_stat, total_stat, fl_stat
    return False, False, False


def contig_stat_manager(writer_stat, contig, initial_nb_lonely, rescue, total_stat, genes, adj_orfs=[]):
    contig_stat = {}
    contig_stat['contig'] = contig
    contig_stat['gene with TA domain'] = len(genes)
    contig_stat['linked gene'] = len(list(fct.get_linked_genes(genes)))
    contig_stat['lonely gene'] = contig_stat['gene with TA domain'] - contig_stat['linked gene']
    if rescue:
        contig_stat['lonely gene rescue'] = initial_nb_lonely - contig_stat['lonely gene']
        contig_stat['adjacent orf'] = len(adj_orfs)
        contig_stat['orf with TA domain'] = len([gene for gene in genes if gene.feature == 'ORF'])

    for k in contig_stat:  # add the value of the row in obj.Gene.metaG_stat to make the total at the end
        total_stat[k] += contig_stat[k]  # contig is not there yet because it is not numerical value
    total_stat['contig'] = contig
    writer_stat.writerow(contig_stat)


def write_result(set_linked, dict_output, contig):
    i = 0
    contig_header = "\n" + '==' * 2 + contig + '==' * 2 + '\n'
    # if dict_output['result_T']:
    #     dict_output['result_T'].write(contig_header)
    if dict_output['result_H']:
        dict_output['result_H'].write(contig_header)
    if dict_output['result_S']:
        dict_output['result_S'].write(contig_header)

    for gene in sorted(set_linked, key=attrgetter('start')):
        for g_post in gene.post:
            i += 1  # to give a number to each TA pair
            if dict_output['result_H']:
                write_human_result(gene, g_post, dict_output['result_H'], i)
            if dict_output['result_S']:
                write_short_result(gene, g_post, dict_output['result_S'], i)
            if dict_output['result_TA_pairs']:
                write_table_pairs_result(gene, g_post, dict_output['result_TA_pairs'])
        if dict_output['result_TA_genes']:
            write_table_genes_result(gene, dict_output['result_TA_genes'])

        if dict_output['result_GFF']:
            write_gff(gene, dict_output['result_GFF'])


def simplify_families_field(families_str):
    families_str = families_str.replace('|', ';')
    families = [family.strip() for family in families_str.split(';')]
    families_count = {family: families.count(family) for family in set(families)}
    simplified_family_str = ';'.join([family for family, count in sorted(
        families_count.items(), key=lambda item: str(item[1])+item[0], reverse=True)])
    return simplified_family_str


def write_table_pairs_result(gene_prev, gene_post, out_tsv_fl):
    line = {}

    line["contig"] = gene_prev.contig
    line["strand"] = gene_prev.strand
    assert gene_prev.strand == gene_post.strand and gene_prev.contig == gene_post.contig

    for i, gene in [(1, gene_prev), (2, gene_post)]:
        line[f"gene{i}_id"] = give_id(gene)
        line[f"gene{i}_start"] = gene.start
        line[f"gene{i}_end"] = gene.end
        line[f"gene{i}_length"] = len(gene)
        line[f"gene{i}_TA_domains"] = ';'.join(
            [d.domain_info['acc'] for d in gene.domain if d.source == "hmmsearch"])
        line[f"gene{i}_TA_families"] = ';'.join(
            [d.domain_info['family'].replace('|', ';') for d in gene.domain if d.source == "hmmsearch"])
        line[f"gene{i}_TADB_hits"] = ';'.join(
            [d.name for d in gene.domain if d.source == "diamond"])

    gene_prev_score = gene_prev.dict_score[gene_post.gene_number]
    gene_post_score = gene_post.dict_score[gene_prev.gene_number]

    line[f"gene1_TA_families"] = simplify_families_field(line[f"gene1_TA_families"])
    line[f"gene2_TA_families"] = simplify_families_field(line[f"gene2_TA_families"])

    line['shared_family'] = ';'.join(sorted(set(line[f"gene1_TA_families"].split(
        ';')) & set(line[f"gene2_TA_families"].split(';'))))
    # print(line[f"gene1_TA_families"])
    # print(line[f"gene2_TA_families"])
    # print('SHARED', line['shared_family'])
    npc = 3
    domain_asso = round(gene_post_score[0]['domain_association_score'], npc)
    dist_score = round(gene_post_score[0]['dist_score'], npc)
    system_score = score.confl(gene_prev_score[0]["len_score"], gene_post_score[0]
                               ["len_score"], gene_post_score[0]['dist_score'], gene_post_score[0]['domain_association_score'])

    line[f"gene1_length_score"] = round(gene_prev_score[0]["len_score"], npc)
    line[f"gene2_length_score"] = round(gene_post_score[0]["len_score"], npc)
    line["distance"] = gene_post_score[0]['distance']
    line["distance_score"] = dist_score
    line["TA_association_score"] = domain_asso
    line["system_score"] = system_score

    out_tsv_fl.writerow(line)


def write_gff(gene, fl):
    neighbors = "possible_partners="+','.join([give_id(n) for n in gene.prev+gene.post])
    list_attrib = ["ID="+give_id(gene), 'best_'+gene.domain[0].writeGffLike(), neighbors]
    attributes = '{}|{};'.format(gene.contig, gene.gene_number)
    attributes += ';'.join(list_attrib)
    list_gff = [gene.contig, 'metaF', gene.feature, str(
        gene.start), str(gene.end), ".", gene.strand, ".", attributes+'\n']
    fl.write('\t'.join(list_gff))


def write_table_genes_result(gene, tsvfl):
    line = {}
    line["contig"] = gene.contig
    # line["gene_number"] = gene.gene_number
    line["gene_id"] = give_id(gene)
    line["start"] = gene.start
    line["end"] = gene.end
    line["length"] = len(gene)
    # line["length_score"] = score[0]["length"]
    line["strand"] = gene.strand  # + gene.frame
    line["feature"] = gene.feature
    line["possible_partners"] = ';'.join([give_id(g) for g in gene.prev + gene.post])
    line["TA_domains"] = ';'.join([d.domain_info['acc']
                                   for d in gene.domain if d.source == "hmmsearch"])
    line["TADB_hits"] = ';'.join([d.name for d in gene.domain if d.source == "diamond"])

    line[f"TA_families"] = ';'.join(
        [d.domain_info['family'].replace('|', ';').strip() for d in gene.domain if d.source == "hmmsearch"])

    line[f"TA_families"] = simplify_families_field(line[f"TA_families"])

    tsvfl.writerow(line)


def write_adj_gene(gene, neighbours, position):
    npc = 3
    info = ''
    for n in neighbours:
        gene_score = gene.dict_score[n.gene_number]
        n_score = n.dict_score[gene.gene_number]
        distance = gene_score[0]["distance"] if 'distance' in gene_score[0] else n_score[0]["distance"]
        dist_score = gene_score[0]["dist_score"] if 'dist_score' in gene_score[0] else n_score[0]["dist_score"]

        info += '{} {} distance {} (score {}) system score {}|'.format(position, give_id(
            n), distance, round(dist_score, npc), round(score.confl(n_score[0]['score'], gene_score[0]['score']), npc))

    return info[: -1]


def write_short_result(g, post, fl, i):
    g_score = g.dict_score[post.gene_number]
    post_score = post.dict_score[g.gene_number]
    # domain_asso = post_score[0]['domain_association_score']
    # system_score = score.confl(g_score[0]['score'], post_score[0]['score'], domain_asso)

    tag_g = give_id(g)
    tag_p = give_id(post)

    system_score = score.confl(g_score[0]['score'], post_score[0]
                               ['score'], post_score[0]['domain_association_score'])

    fl.write("{}. Genes {} & {}\tstrand {}\tscore {}\n".format(
        i, tag_g, tag_p, g.strand, system_score))


def write_human_result(g, post, fl, i):
    npc = 3
    g_score = g.dict_score[post.gene_number]
    post_score = post.dict_score[g.gene_number]
    domain_asso = round(post_score[0]['domain_association_score'], npc)
    system_score = score.confl(g_score[0]['score'], post_score[0]
                               ['score'], post_score[0]['domain_association_score'])
    dist_score = round(post_score[0]['dist_score'], npc)

    fl.write("\nPRE GENE\n" + write_line(g, g_score))
    fl.write("\nPOST GENE\n" + write_line(post, post_score) + '\n')
    fl.write("DISTANCE {} ({})\tDomain-Domain ASSOCIATION: {}\tSYSTEM score: {}\n".format(
        post_score[0]['distance'], dist_score, domain_asso, round(system_score, npc)))
    fl.write(visualisation_genes(g, post, post_score[0]['distance'])+'\n')


def write_line(g, score):
    line = f"{give_id(g)}\tfrom {g.real_start()} to {g.real_end()}"  # ({g.feature}_{g.gene_number})
    line += f"\t{int(score[0]['length']/3)}aa ({score[0]['len_score']:.3})"
    # line += f"\tstart {score[0]['start']}"

    domain_va = g.valid_domain(score[0]['start'])
    domain_va = write_domain_lines(domain_va)
    line += "\nTA HIT:\n{}".format(domain_va)
    return line


def write_domain_lines(domains):
    domain_str_list = []
    do_str = ''
    for d in domains:
        type_tot = sum([n for n in d.domain_info['type_prct'].values()])
        type_prct_type_list = ['{}:{}%'.format(t, int(round(float(n)*100/type_tot)))
                               for t, n in d.domain_info['type_prct'].items()]
        type_prct_type = ' | '.join(sorted(type_prct_type_list))
        # print type_prct_type
        domain_str_list.append([d.domain_info["acc"],
                                d.domain_info['family'], type_prct_type])

    len_max_col0 = max([len(line[0]) for line in domain_str_list])
    len_max_col1 = max([len(line[1]) for line in domain_str_list])

    for line in domain_str_list:

        line[0] += ' '*(len_max_col0 - len(line[0]))
        line[1] += ' '*(len_max_col1 - len(line[1]))
        do_str += '\t'.join(line) + '\n'

    return do_str


def visualisation_genes(pre, post, distance):
    pre_str = visual_str(len(pre))
    post_str = visual_str(len(post))
    dist_str = str(distance) + 'nt'

    sign = (distance + 1) / abs(distance + 1)
    visual_dist = int(distance / 30.0 + 0.98 * sign)

    final_str = pre_str + '\n'
    position_g2 = (len(pre_str) + visual_dist)
    final_str += position_g2 * ' ' + post_str + '\n'
    positions = [position_g2, len(pre_str)]
    final_str += (min(positions) - 1) * ' ' + '/' + abs(positions[0] - positions[1]) * ' ' + '\\\n'
    final_str += int(min(positions) +
                     (abs(positions[0] - positions[1]) - len(dist_str))/2) * ' ' + dist_str

    return final_str


def visual_str(size):
    # min and max size of the gene in characteres on the viual representation
    # vlen = ((length - gene.length_min)) * ((visual_max - visual_min) / (gene.length_max - gene.length_min)) + visual_min
    vlen = int((size / 50) + 1)

    string = vlen * '=' + str(size / 3) + 'aa' + vlen * '=' + '>'
    return string


def give_id(g):
    if hasattr(g, 'gene_id'):
        return g.gene_id  # + ':' + str(g.gene_number)
    else:
        return g.feature+str(g.gene_number)
