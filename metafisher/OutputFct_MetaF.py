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
import logging


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
              "length", "strand", "feature", "possible_partners",
              'lonely_gene', 'TA_domains', 'TA_families', 'TADB_hits',
              "toxin_score", "antitoxin_score", "type"]
    writer_table = csv.DictWriter(fl, fieldnames=header, delimiter='\t')
    writer_table.writeheader()
    return writer_table


def result_ta_pairs_config(fl):
    if not fl:
        return False  # fl is False
    # Chnage the regular file writer into a csv file writer with header
    # header = ["contig",  "strand",
    #           "gene1_id", "gene1_start", "gene1_end", "gene1_length", "gene1_length_score",
    #           "gene2_id", "gene2_start", "gene2_end", "gene2_length", "gene2_length_score",
    #           "distance", "distance_score", "TA_association_score", "system_score",
    #           'gene1_TA_domains', "gene1_TA_families", 'gene1_TADB_hits',
    #           'gene1_best_hit', 'gene1_best_hit_evalue', 'gene1_best_hit_bitscore',
    #           'gene2_TA_domains', "gene2_TA_families", 'gene2_TADB_hits',
    #           'gene2_best_hit', 'gene2_best_hit_evalue', 'gene2_best_hit_bitscore',
    #           'shared_family', 'top_family',
    #           "gene1_toxin_score", "gene1_antitoxin_score", "gene2_toxin_score",
    #           "gene2_antitoxin_score", "gene1_type", "gene2_type", "valid_system"]

    header = ["contig",
              "gene1_id",
              "gene2_id",
              "system_score",
              'gene1_TA_domains',
              'gene1_TADB_hits',
              'gene2_TA_domains',
              'gene2_TADB_hits',
              'top_family',
              "gene1_type",
              "gene2_type",
              "valid_system"]

    writer_table = csv.DictWriter(fl, fieldnames=header, delimiter='\t')
    writer_table.writeheader()
    return writer_table


def output_file_creation(output_way, metaG_name, dict_output, headinfo, complement):
    # the fl of each kind of output are stored in a dictionnary:
    # key : name of the output | value : fl or False if not wanted
    is_output = False  # by default it is False and then if one of the output is True, it will become True

    for out_name, output_flag in dict_output.items():
        if not output_flag:
            continue
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


def write_result(ta_genes, dict_output, contig):
    i = 0
    contig_header = "\n" + '==' * 2 + contig + '==' * 2 + '\n'
    # if dict_output['result_T']:
    #     dict_output['result_T'].write(contig_header)
    if dict_output['result_H']:
        dict_output['result_H'].write(contig_header)
    if dict_output['result_S']:
        dict_output['result_S'].write(contig_header)

    for gene in sorted(ta_genes, key=attrgetter('start')):
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

        # if dict_output['result_GFF']:
            # write_gff(gene, dict_output['result_GFF'])


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
    # line["strand"] = gene_prev.strand
    assert gene_prev.strand == gene_post.strand and gene_prev.contig == gene_post.contig

    for i, gene in [(1, gene_prev), (2, gene_post)]:
        best_hit = gene.get_best_hit()

        line[f"gene{i}_id"] = give_id(gene)
        # line[f"gene{i}_start"] = gene.start
        # line[f"gene{i}_end"] = gene.end
        # line[f"gene{i}_length"] = len(gene)
        line[f"gene{i}_TA_domains"] = ';'.join(
            [d.domain_info['acc'] for d in gene.domains if d.source == "hmmsearch"])
        # line[f"gene{i}_TA_families"] = ';'.join(gene.ta_families)
        line[f"gene{i}_TADB_hits"] = ';'.join(
            [d.name for d in gene.domains if d.source == "diamond"])
        # line[f"gene{i}_best_hit"] = best_hit.domain_info['acc']
        # line[f"gene{i}_best_hit_evalue"] = best_hit.e_value
        # line[f"gene{i}_best_hit_bitscore"] = best_hit.score

        # line[f"gene{i}_toxin_score"] = gene.toxin_score
        # line[f"gene{i}_antitoxin_score"] = gene.antitoxin_score
        line[f"gene{i}_type"] = gene.type

    gene_prev_score = gene_prev.dict_score[gene_post.gene_number]
    gene_post_score = gene_post.dict_score[gene_prev.gene_number]

    # line['shared_family'] = ';'.join(
    #     sorted(set(gene_post.ta_families) & set(gene_prev.ta_families)))

    best_shared_families = get_best_shared_families(gene_prev, gene_post)
    line['top_family'] = ';'.join(sorted(best_shared_families))
    if len(best_shared_families) > 1:
        logging.info(f'More than one best family {line["top_family"]}')

    # print(line[f"gene1_TA_families"])
    # print(line[f"gene2_TA_families"])
    # print('SHARED', line['shared_family'])
    npc = 3
    domain_asso = round(gene_post_score[0]['domain_association_score'], npc)
    dist_score = round(gene_post_score[0]['dist_score'], npc)
    system_score = score.confl(gene_prev_score[0]["len_score"], gene_post_score[0]
                               ["len_score"], gene_post_score[0]['dist_score'], gene_post_score[0]['domain_association_score'])

    # line[f"gene1_length_score"] = round(gene_prev_score[0]["len_score"], npc)
    # line[f"gene2_length_score"] = round(gene_post_score[0]["len_score"], npc)
    # line["distance"] = gene_post_score[0]['distance']
    # line["distance_score"] = dist_score
    # line["TA_association_score"] = domain_asso
    line["system_score"] = system_score
    line['valid_system'] = is_a_valid_system(gene_prev.type, gene_post.type, system_score)
    out_tsv_fl.writerow(line)


def is_a_valid_system(type1, type2, system_score):
    if sorted([type1, type2]) != ["antitoxin", 'toxin']:
        return False
    if system_score < 0.8:
        return False
    return True


def get_best_shared_families(gene1, gene2):
    shared_families = sorted(set(gene1.ta_families) & set(gene2.ta_families))
    # print(set(shared_families))
    # print(gene1, gene2)
    # for do in gene1.domains:
    #     print(do)

    shared_families2sum_score = {
        family: gene1.family2bitscore[family]+gene2.family2bitscore[family] for family in shared_families}
    best_sum_score = max(shared_families2sum_score.values())

    best_families = [family for family, score in shared_families2sum_score.items()
                     if score == best_sum_score]
    return best_families


def write_gff(gene, fl):
    neighbors = "possible_partners="+','.join([give_id(n) for n in gene.prev+gene.post])
    list_attrib = ["ID="+give_id(gene), 'best_'+gene.domains[0].writeGffLike(), neighbors]
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
    line["lonely_gene"] = False if gene.prev or gene.post else True
    line["TA_domains"] = ';'.join([d.domain_info['acc']
                                   for d in gene.domains if d.source == "hmmsearch"])
    line["TADB_hits"] = ';'.join([d.name for d in gene.domains if d.source == "diamond"])

    line[f"TA_families"] = ';'.join(gene.ta_families)

    line["toxin_score"] = gene.toxin_score
    line["antitoxin_score"] = gene.antitoxin_score
    line["type"] = gene.type

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
    line += "\n\nTA HITs:\n{}".format(domain_va)
    return line


def write_domain_lines(domains):
    do_str = ''
    headers = ['accession', 'family', 'TA type', 'Evalue', 'bit score']
    dash_sep = ['-'*len(h) for h in headers]
    domain_str_list = [headers, dash_sep]

    for d in domains:
        type_tot = sum([n for n in d.domain_info['type_prct'].values()])
        type_prct_type_list = ['{}:{}%'.format(t, int(round(float(n)*100/type_tot)))
                               for t, n in d.domain_info['type_prct'].items()]
        type_prct_type = ' | '.join(sorted(type_prct_type_list))
        # print type_prct_type
        domain_str_list.append([d.domain_info["acc"],
                                d.domain_info['family'], type_prct_type, f"{d.e_value}", f"{d.score}"])

    for col_index in range(len(headers)):

        len_max_col = max([len(line[col_index]) for line in domain_str_list])

        for line in domain_str_list:

            line[col_index] += ' '*(len_max_col - len(line[col_index]))

    do_str = '\n'.join(['\t'.join(line) for line in domain_str_list])+'\n'

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
