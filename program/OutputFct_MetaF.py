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
