#!/usr/bin/env python3

"""
Module      : Orf
Description : The main entry point for the program.
Copyright   : (c) Jean Mainguy, 27 nov. 2020
License     : MIT
Maintainer  : jean.mainguy@outlook.fr


Tool to retrieve toxin antitoxin (TA) systems in genomes or metagenomes.
"""

import Object_MetaF as obj
import csv
import re
from subprocess import call
from operator import attrgetter
import os
import logging
from collections import defaultdict, namedtuple
import gzip


def run_command(command, outfile):
    cmd_name = command.split()[0]
    logging.info(f"{cmd_name} command: {command}")
    exit_code = os.system(command)

    if exit_code != 0:
        raise ValueError(f"{cmd_name} has failed. exit_code={exit_code}.")
    # call(hmmsearch_command, shell=True)
    if not os.path.isfile(outfile):
        # logging.warning(f'hmmsearch has failed.')
        raise FileNotFoundError(f'{cmd_name} command has failed to create the file {outfile}')

    return exit_code


def hmmsearch(faa_file, hmm_db, outfile, force=True):

    if os.path.isfile(outfile) and not force:
        logging.info(f'Hmmresult file {outfile} exists already. We use it.')
        return

    else:
        logging.info('Running hmmsearch to identify TA genes.')

    if faa_file.endswith('.gz'):
        logging.info(
            f'Protein sequence file {faa_file} ends with gz. It will be unzipped first to be used in hmmsearch.')

        zipped_faa_file = faa_file
        faa_file = os.path.join(os.path.dirname(outfile), os.path.basename(faa_file)[:-3])
        zcat_command = f"zcat {zipped_faa_file} > {faa_file}"
        run_command(zcat_command, faa_file)

    hmmsearch_command = f"hmmsearch -E 0.5 --domtblout {outfile} {hmm_db} {faa_file} > {outfile}.LOG"
    run_command(hmmsearch_command, outfile)


def diamond_blastp(faa_file, diamond_db, outfile):
    fields = 'qseqid sseqid qstart qend sstart send pident qcovhsp evalue bitscore stitle'
    hmmsearch_command = f"diamond blastp -q {faa_file} -d {diamond_db} -o {outfile} --outfmt 6 {fields} > /dev/null"
    run_command(hmmsearch_command, outfile)


def get_fast_fasta(fl, line, contig):
    seq = {'id': contig, 'data': '', 'rev': False}

    while line and line[:len(contig) + 1] != f'>{contig}':
        line = next(fl)

    for line in fl:
        if line[0] == '>':
            break
        seq['data'] += line.rstrip().upper()

    if seq['data']:
        return seq, line
    else:
        raise ValueError(
            f"seq['data'] is empty.. get_fast_fasta has failed to retreive the sequence of {contig}")


def get_ta_genes_from_diamond(diamond_result,  gene_to_hits=defaultdict(list), min_coverage=95, min_pident=95):

    fields = 'qseqid sseqid qstart qend sstart send pident qcovhsp evalue bitscore stitle'.split()

    DiamondLineParser = namedtuple('DiamondLineParser', fields)

    with open(diamond_result) as fl:
        for line in fl:

            diamondhit = DiamondLineParser._make(line.strip().split('\t'))

            if float(diamondhit.qcovhsp) < min_coverage and float(diamondhit.pident) < min_pident:
                continue

            gene_id = diamondhit.qseqid

            domain = obj.TaHit(
                ali_from=diamondhit.qstart,
                ali_to=diamondhit.qend,
                name=diamondhit.sseqid,
                accession=' '.join(diamondhit.stitle.split(' ')[1:]),
                e_value=diamondhit.evalue,
                score=diamondhit.bitscore,
                line=line,
                source='diamond')

            gene_to_hits[gene_id].append(domain)

    return gene_to_hits


def hmm_result_parser(hmm_result):

    domtblout_headers = ("target_name", "target_accession",
                         "tlen", "query_name", "query_accession",
                         "qlen", "seq_e_value", "seq_score", "seq_bias",
                         "domain_count", "domain_of", "domain_c_Evalue", "domain_i_Evalue",
                         "domain_score", "domain_bias", "hmm_coord_from", "hmm_coord_to",
                         "ali_coord_from", "ali_coord_to", "env_coord_from",
                         "env_coord_to", "acc", "description_of_target")

    HmmLineParser = namedtuple('HmmLineParser', domtblout_headers)

    with open(hmm_result) as fl:
        for line in fl:
            if line.startswith('#'):
                continue
            yield HmmLineParser._make(line.split()[:23])


def get_ta_genes_from_hmmsearch(hmm_result, gene_to_hits=None):

    if gene_to_hits is None:
        gene_to_hits = defaultdict(list)

    for hmmhit in hmm_result_parser(hmm_result):

        gene_id = hmmhit.target_name

        domain = obj.TaHit(
            ali_from=hmmhit.ali_coord_from,
            ali_to=hmmhit.ali_coord_to,
            name=hmmhit.query_name,
            accession=hmmhit.query_accession,
            e_value=hmmhit.seq_e_value,
            score=hmmhit.seq_score,
            line=hmmhit,
            source='hmmsearch')

        gene_to_hits[gene_id].append(domain)

    return gene_to_hits


def get_hmm_genes(contig, table_hmm, gff_file):
    """ """
    genes = []

    # attribut of the classe gene like every objet (orf and TA_gene) will have this attribut, then no need to give it that each time
    # obj.Gene.scaffold = contig
    domains_dict = {}  # keys : gene numbers, value: liste of the domain objet !!
    with open(table_hmm) as fl:
        for line in fl:
            if line[:len(contig)] == contig:
                #  domain is an OBJECT of class TaHit. It gather info about the domain found.
                gene_number, domain = hmmtable_parser(line)

                domains_dict.setdefault(gene_number, []).append(domain)
    gff_headers = ("seqname", "_3", "feature", "start", "end", "_2", "strand", "_1", "attribute")
    with open(gff_file, 'r') as gff_fl:
        gff_tsv_dict_reader = csv.DictReader(gff_fl, delimiter='\t', fieldnames=gff_headers)

        for line_dict in gff_tsv_dict_reader:
            if line_dict["feature"] not in ['CDS', 'ORF', 'gene'] or line_dict["seqname"] != contig:
                continue

            gene_number = int(line_dict['attribute'].split(";")[0].split('|')[1])

            if gene_number in domains_dict:
                gene = build_ta_gene_from_gff_line(line_dict)
                gene.domain = domains_dict[gene_number]

                gene.domain_Ct_border = max((d.ali_from * 3 for d in domains_dict[gene_number]))

                genes.append(gene)

    return genes


def build_ta_gene_from_gff_line(gff_line):
    gene = obj.TA_gene()

    gene.contig = gff_line["seqname"]
    gene.feature = gff_line["feature"]
    gene.start = int(gff_line['start'])
    gene.end = int(gff_line['end'])
    gene.strand = gff_line['strand']
    # gene.gene_number = cds_number
    # gene_id_search = re.search('ID=([^;\n]+)', gff_line['attribute'])
    # gene.gene_id = protein_id_search.group(1)

    protein_id_search = re.search('protein_id=([^;\n]+)', gff_line['attribute'])
    if protein_id_search:
        gene.protein_id = protein_id_search.group(1)

    return gene


def hmmtable_parser(line):
    """
    Parse a line of the hmmtable and return a object from TaHit class with all the info of the line.
    """
    pattern = re.compile(r"""
    (?P<scaffold>[^|]+)\|(?P<gene_number>\d{1,4}) #scaffold name and gene_number separated by a |
    \s+-\s+
    (?P<gene_len>\d+) # length of the gene sequence in residu
    \s+
    (?P<domain>[\w.\(\)_+-]+) # domain name
    \s+
    (?P<domain_acc>[\w.+-]+) # domain acc pour pfam
    \s+\d+\s+
    (?P<evalue>[\w.-]+) # E-value of the overall sequence/profile comparison (including all domains).
    \s+
    (?P<score>[\d.]+) # Bit score of the overall sequence/profile comparison (
    \s+[\d.]+\s+
    (?P<domain_number>\d+) # This domain’s number
    \s+
    (?P<total_domain>\d+) # The total number of domains reported in the sequence, ndom.
    \s+[\w.-]+\s+[\w.+-]+\s+[\w.-]+\s+[\w.-]+\s+\d+\s+\d+\s+
    (?P<ali_from>\d+) # The start of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
    \s+
    (?P<ali_to>\d+) # The end of the MEA alignment of this domain with respect to the sequence, numbered 1..L for a sequence of L residues.
    \s+
    (?P<env_from>\d+) # The start of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
                      # The envelope defines a subsequence for which their is substantial probability mass supporting a homologous domain,
                      # whether or not a single discrete alignment can be identified.
                      # The envelope may extend beyond the endpoints of the MEA alignment, and in fact often does, for weakly scoring domains.
    \s+
    (?P<env_to>\d+) # The end of the domain envelope on the sequence, numbered 1..L for a sequence of L residues.
    \s
    """, re.VERBOSE)
    """
    regex in one line:
    (?P<scaffold>[^|]+)\|(?P<gene_number>\d{1,4})\s+-\s+(?P<gene_len>\d+)\s
    +(?P<domain>[\w.-]+)\s+(?P<domain_acc>[\w.+-]+)\s+\d+\s+(?P<evalue>[\w.-]
    +)\s+(?P<score>[\d.]+)\s+[\d.]+\s+(?P<domain_number>\d+)\s+(?P<total_domain>\d+)\
    s+[\w.+-]+\s+([\w.+-]+)\s+[\w.+-]+\s+[\w.+-]+\s+\d+\s+\d+\s+(?P<ali_from>\d+)\s+(?P
    <ali_to>\d+)\s+(?P<env_from>\d+)\s+(?P<env_to>\d+)\s

    link : https://regex101.com/r/kKHzsB/2
    """

    result = pattern.match(line)
    if not result:
        raise Exception(
            'The parser of HMM table has failed... you may check the HMM table output ---> hmmtable line with a problem', line)

    domain = obj.TaHit(
        domain_number=result.group("domain_number"),
        ali_from=result.group("ali_from"), ali_to=result.group("ali_to"), env_from=result.group("env_from"),
        env_to=result.group("env_to"), name=result.group("domain"),
        accession=result.group("domain_acc"), e_value=result.group("evalue"), score=result.group("score"), line=line)
    gene_number = int(result.group("gene_number"))
    return gene_number, domain


def get_start_po(fna_seq_dict, genes):
    """
    Process the search of start position over the list of genes.
    """

    non_valides = []
    for gene in sorted(genes, key=attrgetter('gene_number')):
        start_po = gene.give_start_po(fna_seq_dict, min_max_intervall=True)

        if not start_po:
            # print 'invalide', gene
            non_valides.append(gene)
            continue

        gene.distanceMin = obj.Gene.distanceMin - abs(start_po[-1] - start_po[0])

        # WARNING absolute value of distance min should not be greater than the length of the gene !!
        # because it would allow overlap of more than the length of the gene
        # and then could give a post gene before a pre gene so a non sense

        if abs(gene.distanceMin) >= len(gene):
            gene.distanceMin = -len(gene) + 1  # if so distance min is len -1

        # Uncomment this part if instead of having the position of start inside the intial gene,
        # it is desire to have the position inthe contig

        if gene.strand == '+':
            # start_po = [x + gene.start for x in start_po]
            gene.start += start_po[0]

        else:
            # start_po = [gene.end - x for x in start_po]
            gene.end -= start_po[0]

        gene.possible_starts = start_po

    for gene in non_valides:
        genes.remove(gene)


def get_contig_from_hmmsearch_result(table_hmm):
    # input : Name of the table hmm file
    # output : List of all the scaffold of the table

    set_of_scaffold = set()
    fl_hmm = open(table_hmm, 'r')
    for line in fl_hmm:
        if line[0] == "#":
            continue
        set_of_scaffold.add(line[:line.index('|')])
    fl_hmm.close()
    return list(set_of_scaffold)


def compute_gene_adjacency(genes):
    """
    Treat independently list plus and minus
    Sort the list by the end position
    then take the first one and look for the next gene (the gene has the its end after the end of the first one).
    call tandem gene function to check if the gene fit the distance threshold
    """

    genes_strand_plus = sorted(
        (gene for gene in genes if gene.strand == '+'), key=attrgetter('end'))
    genes_strand_minus = sorted(
        (gene for gene in genes if gene.strand == '-'), key=attrgetter('start'))

    # adj_by_strand(obj.TA_gene.genes_plus)
    # adj_by_strand(obj.TA_gene.genes_minus)
    linked_genes = adj_by_strand(genes_strand_plus)
    linked_genes |= adj_by_strand(genes_strand_minus)

    return linked_genes


def adj_by_strand(genes):
    """
    liste: list of hmm gene with homogenous strand
    Check if the gene is in tandem with another and if so store the gene inside a set obj.TA_gene.linked
    In parallel it clean up the list obj.TA_gene.genes
    by removing the genes that forme a tandem. Then TA_gene.genes has only the lonely_gene
    """
    linked_genes = set()

    for gi, gene in enumerate(genes):
        # print obj.TA_gene.genes_plus[gi].gene_number,  obj.TA_gene.genes_plus[gi].len_val
        for gpost in genes[gi + 1:]:
            if gpost.end - gene.end + 1 > obj.Gene.length_max + obj.Gene.distanceMax:
                """
                if the distance between gene.end and gpost.end is superior to lenmax + distmax
                Then the two gene won't be in pair and the next postgene either because they are sorted by their start
                So we can break the gpost for loop and check the next gene
                """
                break
            # it is a simple test that ckeck if the two gene are adjacent
            if gene.is_pre_adj_to(gpost):
                #  store the information of prev and post according the strand
                if gene.strand == '+':
                    gene.post.append(gpost)
                    gpost.prev.append(gene)
                else:
                    gpost.post.append(gene)
                    gene.prev.append(gpost)
                # add the gene because it has a link in the set linked of class TA_gene
                linked_genes.add(gene)
                # add the gene because it has a link in the set linked of class TA_gene
                linked_genes.add(gpost)
    return linked_genes


def check_size(genes):
    """
    Remove genes with a length that does not fit the length thresholds.
    """
    invalide_genes = []
    for gene in genes:
        if obj.Gene.length_min <= len(gene) <= obj.Gene.length_max:
            gene.possible_starts = [0]

            # WARNING absolute value of distance min should not be greater than the length of the gene !!
            # because it would allow overlap of more than the length of the gene
            # and then could give a post gene before a pre gene so a non sense
            if abs(gene.distanceMin) >= len(gene):
                logging.warning(f'gene distance min >= length of gene')
                gene.distanceMin = -len(gene) + 1  # if so distance min is len -1

        else:
            invalide_genes.append(gene)

    for gene in invalide_genes:
        genes.remove(gene)


def get_linked_genes(genes):
    return (gene for gene in genes if gene.post or gene.prev)


def get_lonely_genes(genes):
    """
    Retrieve gene that are not linked to any other genes.
    """
    return (gene for gene in genes if not gene.post and not gene.prev)


def delete_files(listeFiles):
    for outfile in listeFiles:
        try:
            os.remove(outfile)
        except OSError as e:  # if failed, report it back to the user ##
            print(f"Error: {e.filename} - {e.strerror}.")


def annotate_ta_hits(domains, info_domains, dict_domain_gene_type):
    for domains_grp in domains:
        for domain in domains_grp:
            domain.annotate_ta_hit(info_domains, dict_domain_gene_type)


def check_protein_ids(gene_to_domains):
    """
    Some ncbi protein files have protein id as follow: lcl|<seqname>_prot_<prot_id>_<prot_rank_in_seqname>
    ie:lcl|NC_010645.1_prot_WP_012415714.1_3
    Let's identify these ids and keep only the protein id.
    """
    gene_to_domains_cleaned = {}
    pattern = re.compile("lcl\|[A-Z]{2}_\d+\.\d+_prot_([A-Z]{2}_\d+\.\d+)_\d+")
    for gene_id, domains in gene_to_domains.items():
        if pattern.match(gene_id):
            id = pattern.match(gene_id).group(1)
        else:
            id = gene_id
        gene_to_domains_cleaned[id] = domains

    return gene_to_domains_cleaned


def get_genes_by_contigs(gene_to_domains, gff_file):

    gene_id_retrieved = []
    genes = []

    current_seqname = ""
    gene_count_by_contig = 0

    gff_headers = ("seqname", "_3", "feature", "start", "end", "_2", "strand", "_1", "attribute")
    gene_id_pattern = re.compile('ID=([^;\n]+)')
    protein_id_pattern = re.compile('protein_id=([^;\n]+)')

    gene_to_domains = check_protein_ids(gene_to_domains)
    proper_open = gzip.open if gff_file.endswith('.gz') else open

    with proper_open(gff_file, 'rt') as gff_fl:
        gff_tsv_dict_reader = csv.DictReader(gff_fl, delimiter='\t', fieldnames=gff_headers)

        for line_dict in gff_tsv_dict_reader:
            seqname = line_dict["seqname"]
            if line_dict["seqname"] != current_seqname:
                gene_count_by_contig = 0
                if genes:
                    yield (current_seqname, genes)
                    genes = []
                current_seqname = line_dict["seqname"]

            if line_dict["feature"] not in ['CDS', 'ORF']:
                continue
            gene_count_by_contig += 1

            id_search = protein_id_pattern.search(line_dict['attribute'])
            if not id_search:
                id_search = gene_id_pattern.search(line_dict['attribute'])

            try:
                gene_id = id_search.group(1)
            except AttributeError as err:
                logging.critical(
                    f'Attribute of gff line has no protein_id or ID: {line_dict["attribute"]}')
                raise err

            # check for lcl protein id from ncbi faa file ie >lcl|NC_010645.1_prot_WP_012415716.1_5

            lcl_prot_id = f"lcl|{seqname}_prot_{gene_id}_{gene_count_by_contig}"

            if lcl_prot_id in gene_to_domains:
                gene_id = lcl_prot_id
            if gene_id in gene_to_domains:
                gene_id_retrieved.append(gene_id)

                gene = build_ta_gene_from_gff_line(line_dict)
                gene.gene_number = gene_count_by_contig
                gene.domain = gene_to_domains[gene_id]
                gene.domain_Ct_border = max((d.ali_from * 3 for d in gene_to_domains[gene_id]))
                genes.append(gene)

        if genes:
            yield (line_dict["seqname"], genes)

    if len(gene_id_retrieved) < len(gene_to_domains):
        missing_prot = set(gene_to_domains) - set(gene_id_retrieved)

        raise ValueError(
            f"Only {len(gene_id_retrieved)}/{len(gene_to_domains)} proteins identified by hmmsearch or diamond have been found in the gff file ({gff_file}). Missing proteins: {missing_prot}")
