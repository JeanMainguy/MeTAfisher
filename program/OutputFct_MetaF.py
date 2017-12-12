
def output_manager(output_way, metaG_name, rescue, resize, thresholds, rescue, resize)

def output_header_creation(metaG_name, rescue,resize, thresholds, rescue, resize):
    complement = ''
    if rescue:
        complement += '_rescue'
    if resize:
        complement += '_resize'
    header = '## Name of the sequence analysed: ' + metaG_name
    header = "## Rescue lonely gene : {}\n".format(rescue)
    header += "## Resize gene : {}\n".format(resize)
    header += "## Distance threshold from {}nt to {}nt\n".format(thresholds['distanceMin'], thresholds['distanceMax'])
    header += "## Length threshold from {}aa to {}aa\n".format(thresholds['lenMin'], thresholds['leneMax'])
    return header, complement

def output_file_creation(output_human, output_short, output_table, header):
    # the fl of each kind of output are stored in a dictionnary:
    # key : name of the output | value : fl or False if not wanted
    dict_output = {'result_H': output_human, "result_S": output_human, 'result_T': output_table, 'is_output': False}

    for out_name in dict_output:
        if out_name:  # if the flag is not False
            file_out = '{}/{}_{}{}.txt'.format(output_way, metaG_name, out_name, complement)
            flout = open(file_out, "w")
            flout.write("## {}\n".format(out_name))
            flout.write(header)
            dict_output[out_name] = flout
            dict_output['is_output'] = True

def stat_file_creation(info_contig_stat, header):
    if info_contig_stat:

        header = ['contig', 'gene with TA domain', 'lonely gene', 'linked gene']
        if rescue:
            header += ['adjacent orf', 'rescue flag', 'orf with TA domain', 'lonely gene rescue']
        fl_stat = open('{}/{}_contig_stat{}.csv'.format(output_way, metaG_name, complement), 'w')
        fl_stat.write("#Rescue lonely gene : {}\n".format(rescue))
        fl_stat.write("#Resize gene : {}\n".format(resize))
        writer_stat = csv.DictWriter(fl_stat, fieldnames=header, delimiter='\t')
        writer_stat.writeheader()
        obj.Gene.metaG_stat = dict.fromkeys(header, 0)
        obj.Gene.metaG_stat['contig'] = metaG_name
