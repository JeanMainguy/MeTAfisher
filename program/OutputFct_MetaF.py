# coding: utf-8
import csv

def output_manager(output_way, metaG_name, thresholds, dict_output, info_contig_stat, rescue, resize):
    headinfo, complement = output_headinfo_creation(metaG_name, thresholds, rescue, resize)
    output_file_creation(output_way, metaG_name, dict_output, headinfo, complement) # dico output is update no need to give it back
    writer_stat, total_stat, fl_stat = stat_file_creation(output_way, metaG_name, info_contig_stat, headinfo, complement, rescue)
    return writer_stat, total_stat, fl_stat

def output_headinfo_creation(metaG_name, thresholds, rescue, resize):
    complement = ''
    if rescue:
        complement += '_rescue'
    if resize:
        complement += '_resize'
    headinfo = '## Name of the sequence analysed: ' + metaG_name
    headinfo = "## Rescue lonely gene : {}\n".format(rescue)
    headinfo += "## Resize gene : {}\n".format(resize)
    headinfo += "## Distance threshold from {}nt to {}nt\n".format(thresholds['distanceMin'], thresholds['distanceMax'])
    headinfo += "## Length threshold from {}aa to {}aa\n".format(thresholds['lenMin'], thresholds['lenMax'])
    return headinfo, complement

def output_file_creation(output_way, metaG_name, dict_output, headinfo, complement):
    # the fl of each kind of output are stored in a dictionnary:
    # key : name of the output | value : fl or False if not wanted
    is_output = False # by default it is False and then if one of the output is True, it will become True

    for out_name in dict_output:

        if out_name:  # if the flag is not False
            file_out = '{}/{}_{}{}.txt'.format(output_way, metaG_name, out_name, complement)
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
            header += ['adjacent orf', 'rescue flag', 'orf with TA domain', 'lonely gene rescue']

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
    contig_stat['gene with TA domain'] = len(obj.TA_gene.genes)
    contig_stat['linked gene'] = len(obj.TA_gene.linked)
    contig_stat['lonely gene'] = contig_stat['gene with TA domain'] - contig_stat['linked gene']
    if rescue:
        contig_stat['lonely gene rescue'] = initial_nb_lonely - contig_stat['lonely gene']
        contig_stat['adjacent orf'] = obj.Orf.adj_orf_index
        contig_stat['orf with TA domain'] = len(obj.Orf.hmm_orf)

    for k in contig_stat: #add the value of the row in obj.Gene.metaG_stat to make the total at the end
        total_stat[k] += contig_stat[k]  # contig is not there yet because it is not numerical value
    total_stat['contig'] = scaffold
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
