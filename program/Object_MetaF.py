# coding: utf-8
import re
import csv
from operator import attrgetter, itemgetter
# import find_orf as orf
import Orf_MetaF as orf2
from math import log
import Function_MetaF as fct2


class Gene:
    # Mettre des stats.. ici en mémoire qui enregistre tout
    # Nb de gene avec un hit de HMM
    # nb d'orf donné à
    #  Class attribut

    # Metafisher Parameters :
    # --> treshold size of gene
    length_min = None  # 26 * 3
    length_max = None  # 500 * 3
    # --> threshold distance for pair
    distanceMin = None  # -100
    distanceMax = None  # 300

    # Allowance of integrity loss of domain (5%) arbitrary number
    allowance = 0.05  # IS not USE anyway...
    scaffold = None  # contig name
    highestGeneNumber = 0
    metaG_stat = {}
    domain_dict = {}

    def __init__(self):
        self.start = None
        self.end = None
        self.possible_start = []
        self.strand = None
        self.feature = None
        self.dict_score = {}
        self.locus_tag = None

    def __str__(self):
        presentation = 'length: {}\tfrom {} to {}\tstrand: {}\tfeature: {}\n'.format(
            self.length(), self.start, self.end, self.strand, self.feature)
        return presentation

    # def __eq__(self, other):
    #     if self.gene_number == other.gene_number:
    #         return True
    #     return False

    def length(self):
        return self.end - self.start + 1

    def __len__(self):  # MAGIC method that the same as length(self) but in a cooler way
        return self.end - self.start + 1

    def real_end(self):
        if self.strand == '+':
            return self.end
        else:
            return self.start

    def real_start(self):
        if self.strand == '+':
            return self.start
        else:
            return self.end

    def distance(self, gene):
        return max(gene.real_end() - self.real_start(), self.real_end() - self.real_start())

    def is_pre_adj_to(self, gene):
        """
        Check if self is in tandem with gene in parameter.
        end of gene is located a front of end of self.------===self==>--===gene(post)===>------
        however for strand minus the gene post is inversed : ------<===self(post)==--<===gene===------
        the if statement is there to take care of that
        """

        if gene.distanceMin < gene.start - self.end < Gene.distanceMax:  # If the distance between the two gene are fitting the thresholds?
            return True
        # elif gene.start - self.end > Gene.dist_max:
            # return False  # false le gene start trop loin
        return False  # Plus simple il me semble qu'avant

    def valid_domain(self, start):
        # Test if the domain is integre avec le start donné en argument !!
        # retourne la liste des domaines qui ne sont pas integre.

        start = start / 3  # transformation en start aa
        valid = []
        for do in self.domain:
            # print 'Domain',do.ali_from,'-->',do.ali_to
            if start <= do.ali_from + (do.ali_to - do.ali_from) * Gene.allowance:
                # print 'domain discard'
                valid.append(do)
        return valid


class TA_gene(Gene):
    genes = []
    genes_strand = {'+': [], '-': []}
    linked = set()
    lonely = None  # {'+': [], '-': []}

    def __init__(self):

        # info of the gff line
        self.gene_number = None
        self.feature = None

        # info of hmm lines
        self.domain = []
        self.best_domain = []
        # self.rasta_domain = []
        # border of the last domain in C terminal! In order to use alternative
        # start that doesn't affect at least one domain.
        self.domain_Ct_border = 0

        # Validity features
        # self.len_val = None
        self.post = []
        self.prev = []
        # N'est pas hyper utile pour le moment....
        self.statue = None  # statue dit si le gene valide ou non.. il est pas valide dans le cas ou le resize ne marche pas
        # score:
        self.dict_score = {}

    def __str__(self):
        presentation = 'nb: {}\tlength: {}\tfrom {} to {}\tstrand: {}\tfeature: {}\n'.format(
            self.gene_number, self.length(), self.start, self.end, self.strand, self.feature)
        presentation += 'nb_domain: {}\tbest: {}\t Possible post_partner: {} pre partner {}'.format(
            len(self.domain), len(self.best_domain), [p.gene_number for p in self.post], [p.gene_number for p in self.prev])
        return presentation

    def give_start_po(self, dico, min_max_intervall=False):
        """
        Give back the list of all the position of the start codon in the dna sequence of the gene.
        The position are relatif to the first start
        If the gene is partial and don't start with a start codon the function add the position 0 as a start

        Adding a feature to this fct :
        Possibility to find start codons only in the correct interval using the min and max length thrshold
        and taking also into account that it is necessary to have at least one domain safe in the sequence
        using attribut domain_Ct_border
        """

        if self.feature == 'ORF':
            fl_name = 'fl_orf'
            line_name = 'line_orf'
        else:
            fl_name = 'fl'
            line_name = 'line'
        start_codon = dico["codon_start"]

        scaffold_gnb = self.scaffold + '|' + str(self.gene_number)  # we have to build this because of fast_fasta function
        seq, dico[line_name] = fct2.get_fast_fasta(dico[fl_name], dico[line_name], scaffold_gnb)  # use of get fast fasta as before.
        # print seq['description']
        # print seq["data"]
        """Setting min and max parameter : where to start and end the start dectection
        in order to have ibly start position fitting the length threshold and giving at least one domain intact"""
        if min_max_intervall:
            len_seq = len(seq['data'])

            # INF
            if len_seq <= Gene.length_max:
                inf = 0
            else:
                inf = len_seq - Gene.length_max

            if len_seq - Gene.length_min > 0:
                sup = min((len_seq - Gene.length_min), self.domain_Ct_border)
                # print 'self.domain_Ct_border', self.domain_Ct_border
            else:
                print(len_seq, len_seq - Gene.length_min)
                print('the length of the gene is smaller than the minimal length threshold ')
                return []  # if the length of the gene is smaller than the minimal length threshold then no start are compatible and then the fct return an empty fct
        else:
            inf = 0
            sup = len(seq['data'])
        # print self.gene_number, 'a une sequence de', len(seq['data']), 'sup and inf', sup, inf

        start_po = orf2.codon_finder(start_codon, seq["data"][inf:sup + 1])  # give evry start position of the sequence

        if inf == 0 and (len(start_po) == 0 or start_po[0] != 0):
            start_po.insert(0, 0)

        return start_po
        # if len(start_po) == 0:
        #     if start_po[0] == 0 or inf != 0:
        #         """
        #         if it starts already at 0 or if the search has not started  at 0 then,
        #         no need to add 0. as a potential start
        #         """
        #         return start_po
        #     else:
        #         """ even if 0 is not a start codon it is still a start position
        #         (e.i. if the gene is at the border of the contig) """
        #         return start_po.insert(0, 0)
        # else:
        #     return []


class Orf(Gene):
    # dict with seq, header of the fasta
    seq = {'data': 'ADCBD'}  # scaffold name
    adj_orf = {}
    adj_orf_index = 0
    hmm_orf = {}

    def __init__(self, frame, possible_start_orf, end_orf, complet=True, border=False):
        """
        end_orf and possible_start_orf in the context of the orf
        meaning position on the personal strand of the gene and the sequence start with 0
        """
        # print "possible_start OF ORF", possible_start_orf
        # possible start in the context of the orf
        self.end_orf = end_orf
        self.start_orf = possible_start_orf[0]  # should be start
        self.frame = frame
        self.possible_start_orf = possible_start_orf
        self.possible_start = [s - self.start_orf for s in possible_start_orf]
        # print 'self.possible_start_orf ', self.possible_start_orf
        # print 'self.possible_start ', self.possible_start
        # Conversion to get start and stop found with official orf finder
        # and to match with predicted gene
        if self.frame < 0:
            len_contig = len(Orf.seq['data'])
            # self.possible_start = [len_contig - x for x in possible_start_orf]

            self.start = len_contig - end_orf  # ncbi start

            self.end = len_contig - possible_start_orf[0]  # ncbi_end
            self.strand = '-'
            # self.gff_start = self.ncbi_end # USE real start instead
        else:
            # self.possible_start = [x + 1 for x in possible_start_orf]
            self.end = end_orf + 1
            self.start = possible_start_orf[0] + 1
            self.strand = '+'

        self.complet = complet
        self.border = border
        self.feature = 'ORF'
        self.dict_score = {}

    def is_adj(self, liste_genes):
        """
        liste_genes : list of gene objet
        Check if self is adjacent to at least one of the gene in liste_genes
        The two gene need to have a personalized distanceMin attribute to be fully accurate on the adjacency
        """
        for gene in liste_genes:
            if gene.is_pre_adj_to(self) or self.is_pre_adj_to(gene):
                return True
        return False

    def manage_size(self):
        """
        TO DELETE !!
        possible start is not start relativ the the realstart  of the gene and not according the contig..
        So fct HASBEEN !
        """
        if len(self) > Gene.length_max:
            print "MANAGE SiZe METHOD"
            print self
            raw_input("size > max")
            del self.possible_start[0]
            while self.end - self.possible_start[0] > Gene.length_max:
                print len(self), self.possible_start[0]
                del self.possible_start[0]  # delete the start
            if self.strand == '+':
                self.start = self.possible_start[0]
            else:
                self.end = self.possible_start[0]
        self.distanceMin = Gene.distanceMin - abs(self.possible_start[-1] - self.possible_start[0])

    def write_faa(self, fl, index):
        """
        Depend of the seq in memory that why it is not a method of Gene
        """
        # print 'IN write FAA'
        if (Orf.seq["rev"] is not True and self.strand == '-') or (Orf.seq["rev"] is True and self.strand == '+'):
            # print 'Need to be reverse comp... seq is rev?', Orf.seq["rev"], ' strand :', self.strand
            Orf.seq = orf2.complement_reverse(Orf.seq)
            # print 'Need to be reverse comp... seq is rev?', Orf.seq["rev"], ' strand :', self.strand
        # print 'seq is rev? ', Orf.seq["rev"], 'self is ', self.strand
        fl.write('>{}|{}\n'.format(self.scaffold, index))
        fl.write(re.sub("(.{60})", "\\1\n", self.protein(), 0, re.DOTALL) + '\n')

    def data_nt(self):
        return self.seq['data'][self.start_orf:self.end_orf + 1]  # !!!### Verifier !

    def protein(self):
        seq = {'data': self.data_nt()[3:]}
        return 'M' + orf2.translate(seq)['data']  # ATTENTION every prot seq start with a M...


class Domain:
    # Threshold overlap !!
    threshold = 0.1  # 10% of overlaping after that  the domains are considered as overlaping

    def __init__(self, domain_number, ali_from, ali_to, env_from, env_to, domain_name, domain_acc, e_value, score, line):
        self.domain_number = int(domain_number)
        self.ali_from = int(ali_from)
        self.ali_to = int(ali_to)
        self.env_from = int(env_from)
        self.env_from = int(env_from)
        self.domain_name = domain_name
        self.domain_acc = domain_acc
        self.e_value = float(e_value)
        self.score = float(score)
        self.line = line

    def __str__(self):
        # return self.name + '\nfrom %d to %d \n' % (self.ali_from * 3, self.ali_to * 3)
        try:
            dico_do = Gene.domain_dict[self.domain_name]
        except KeyError:
            return 'Domain not found in the csv domain file.. {} ({}) '.format(self.domain_name, log(self.score))
        return "name:{} score:{} type:{} family:{}".format(dico_do["acc"], log(self.score), dico_do['type'], dico_do['family'])

    def score_transformed(self):
        # log transformation
        return log(self.score)

    def overlap(self, d):
        if self.ali_from <= d.ali_from <= self.ali_to:  # overlap
            if self.ali_to - d.ali_from > (self.ali_to - self.ali_from) * Domain.threshold:
                return True
            else:
                return False
        elif self.ali_from <= d.ali_to <= self.ali_to:  # overlap
            if d.ali_to - self.ali_from > (self.ali_to - self.ali_from) * Domain.threshold:
                return True
            else:
                return False
        elif d.ali_from <= self.ali_from and self.ali_to < d.ali_to:  # overlap est self anglobé dans d
            return True
        else:  # overlap pas
            return False
