# coding: utf-8
import csv
import Object_MetaF as obj
from operator import attrgetter
import Score_MetaF as score


def output_manager(output_way, metaG_name, thresholds, dict_output, info_contig_stat, rescue, resize):
    headinfo, complement = output_headinfo_creation(metaG_name, thresholds, rescue, resize)
    # Output H, T and S initialization
    # dico output is update no need to give it back
    output_file_creation(output_way, metaG_name, dict_output, headinfo, complement)
    # Output T as a csv file initialization
    dict_output["result_T"] = result_table_config(dict_output["result_T"])
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
    headinfo = '## Name of the sequence analysed: ' + metaG_name
    headinfo += "\n## Rescue lonely gene : {}\n".format(rescue)
    headinfo += "## Resize gene : {}\n".format(resize)
    headinfo += "## Distance threshold from {}nt to {}nt\n".format(
        thresholds['distanceMin'], thresholds['distanceMax'])
    headinfo += "## Length threshold from {}aa to {}aa\n".format(
        thresholds['lenMin'], thresholds['lenMax'])
    return headinfo, complement


def result_table_config(fl):
    if not fl:
        return False  # fl is False
    # Chnage the regular file writer into a csv file writer with header
    header = ["Contig", "Gene number", "Gene id", "start", "end", "length",
              "length_score", "strand", "feature", "domain", "Neighbor gene"]
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
            if out_name == 'result_T':
                extension = 'csv'
            elif out_name == 'result_GFF':
                extension = 'gff'
                # out_name = 'TA_Genes'

            file_out = '{}/{}_{}{}.{}'.format(output_way, metaG_name,
                                              out_name, complement, extension)
            flout = open(file_out, "w")
            flout.write("## {}\n".format(out_name))
            flout.write(headinfo)
            dict_output[out_name] = flout
            is_output = True
    dict_output['is_output'] = is_output


def stat_file_creation(output_way, metaG_name, info_contig_stat, headinfo, complement, rescue):

    if info_contig_stat:
        header = ['contig', 'gene with TA domain', 'lonely gene', 'linked gene']

        if rescue:
            header += ['adjacent orf', 'orf with TA domain', 'lonely gene rescue']

        fl_stat = open('{}/{}_contig_stat{}.csv'.format(output_way, metaG_name, complement), 'w')
        fl_stat.write(headinfo)
        writer_stat = csv.DictWriter(fl_stat, fieldnames=header, delimiter='\t')
        writer_stat.writeheader()
        total_stat = dict.fromkeys(header, 0)
        total_stat['contig'] = metaG_name
        return writer_stat, total_stat, fl_stat
    return False, False, False


def contig_stat_manager(writer_stat, scaffold, initial_nb_lonely, rescue, total_stat):
    contig_stat = {}
    contig_stat['contig'] = scaffold
    contig_stat['gene with TA domain'] = len(obj.TA_gene.genes)
    contig_stat['linked gene'] = len(obj.TA_gene.linked)
    contig_stat['lonely gene'] = contig_stat['gene with TA domain'] - contig_stat['linked gene']
    if rescue:
        contig_stat['lonely gene rescue'] = initial_nb_lonely - contig_stat['lonely gene']
        contig_stat['adjacent orf'] = obj.Orf.adj_orf_index
        contig_stat['orf with TA domain'] = len(obj.Orf.hmm_orf)

    for k in contig_stat:  # add the value of the row in obj.Gene.metaG_stat to make the total at the end
        total_stat[k] += contig_stat[k]  # contig is not there yet because it is not numerical value
    total_stat['contig'] = scaffold
    writer_stat.writerow(contig_stat)


def write_result(set_linked, dict_output, scaffold):
    i = 0
    contig_header = "\n" + '==' * 2 + scaffold + '==' * 2 + '\n'
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
        if dict_output['result_T']:
            write_table_result(gene, dict_output['result_T'])

        if dict_output['result_GFF']:
            write_gff(gene, dict_output['result_GFF'])


def write_gff(gene, fl):
    neighbors = "possible_partners="+','.join([give_id(n) for n in gene.prev+gene.post])
    list_attrib = ["protein_id="+give_id(gene), 'best_'+gene.domain[0].writeGffLike(), neighbors]
    attributes = '{}|{};'.format(gene.scaffold, gene.gene_number)
    attributes += ';'.join(list_attrib)
    list_gff = [gene.scaffold, 'metaF', gene.feature, str(
        gene.start), str(gene.end), ".", gene.strand, ".", attributes+'\n']
    fl.write('\t'.join(list_gff))


def write_table_result(gene, csvfl):
    line = {}
    line["Contig"] = gene.scaffold
    line["Gene number"] = gene.gene_number
    line["Gene id"] = give_id(gene)
    line["start"] = gene.start
    line["end"] = gene.end
    line["length"] = len(gene)
    # line["length_score"] = score[0]["length"]
    line["strand"] = gene.strand  # + gene.frame
    line["feature"] = gene.feature
    line["domain"] = writeDomain(gene)
    line["Neighbor gene"] = write_adj_gene(gene, gene.prev, 'Down:')
    line["Neighbor gene"] += write_adj_gene(gene, gene.post, 'Up:  ')
    csvfl.writerow(line)


def writeDomain(gene):
    """
    return only best domain
    the attribut domain is a sorted list according the score. Then the first domain is the best domain
    """

    # print gene.valid_domain(gene.start)
    return gene.domain[0]


def write_adj_gene(gene, neighbours, position):
    npc = 3
    info = ''
    for n in neighbours:
        # print gene.dict_score
        gene_score = gene.dict_score[n.gene_number]
        n_score = n.dict_score[gene.gene_number]
        # print n_score
        # print gene_score
        # print gene_score[0]['distance']
        distance = gene_score[0]["distance"] if 'distance' in gene_score[0] else n_score[0]["distance"]
        dist_score = gene_score[0]["dist_score"] if 'dist_score' in gene_score[0] else n_score[0]["dist_score"]

        info += '{}\t{}\t distance {} (score {})\tsystem score {}|'.format(position, give_id(
            n), distance, round(dist_score, npc), round(score.confl(n_score[0]['score'], gene_score[0]['score']), npc))

    return info[:-1]


def write_short_result(g, post, fl, i):
    g_score = g.dict_score[post.gene_number]
    post_score = post.dict_score[g.gene_number]
    # domain_asso = post_score[0]['domain_association_score']
    # system_score = score.confl(g_score[0]['score'], post_score[0]['score'], domain_asso)

    tag_g = give_id(g)
    tag_p = give_id(post)
    fl.write("{}. Genes {} & {}\tstrand {}\tscore {}\n".format(
        i, tag_g, tag_p, g.strand, score.confl(g_score[0]['score'], post_score[0]['score'])))


def write_human_result(g, post, fl, i):
    npc = 3
    g_score = g.dict_score[post.gene_number]
    post_score = post.dict_score[g.gene_number]
    domain_asso = round(post_score[0]['domain_association_score'], npc)
    system_score = score.confl(g_score[0]['score'], post_score[0]['score'], domain_asso)
    dist_score = round(post_score[0]['dist_score'], npc)

    fl.write("\n\nPRE GENE\n" + write_line(g, g_score))
    fl.write("\nPOST GENE\n" + write_line(post, post_score) + '\n')
    fl.write("DISTANCE {} ({})\tDomain-Domain ASSOCIATION: {}\tSYSTEM score: {}\n".format(
        post_score[0]['distance'], dist_score, domain_asso, round(system_score, npc)))
    fl.write(visualisation_genes(g, post, post_score[0]['distance']))


def write_line(g, score):
    npc = 3  # number post coma
    line = "Gene {}\tfrom {} to {}\t{}aa ({})\tstart {}\t{}".format(
        g.gene_number, g.real_start(), g.real_end(), int(score[0]["length"] / 3), round(score[0]["len_score"], npc), score[0]["start"], g.feature)

    domain_va = g.valid_domain(score[0]['start'])
    domain_va = write_domain_lines(domain_va)
    line += "\nDOMAIN:\n{}".format(domain_va)
    return line


def write_domain_lines(domains):
    domain_str_list = []
    do_str = ''
    for d in domains:
        type_tot = sum([n for n in d.dict_info['type_prct'].values()])
        type_prct_type_list = ['{}:{}%'.format(t, int(round(float(n)*100/type_tot)))
                               for t, n in d.dict_info['type_prct'].items()]
        type_prct_type = ' | '.join(sorted(type_prct_type_list))
        # print type_prct_type
        domain_str_list.append([d.dict_info["acc"],
                                d.dict_info['family'], type_prct_type])

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
    final_str += int(-1 + min(positions) +
                     (abs(positions[0] - positions[1]) - len(dist_str))/2) * ' ' + dist_str
    # print final_str
    return final_str


def visual_str(size):
    # min and max size of the gene in characteres on the viual representation
    # vlen = ((length - gene.length_min)) * ((visual_max - visual_min) / (gene.length_max - gene.length_min)) + visual_min
    vlen = int((size / 50) + 1)

    string = vlen * '=' + str(size / 3) + 'aa' + vlen * '=' + '>'
    return string


def give_id(g):
    if hasattr(g, 'protein_id'):
        return g.protein_id  # + ':' + str(g.gene_number)
    else:
        return g.feature+str(g.gene_number)
