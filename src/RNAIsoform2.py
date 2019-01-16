from collections import defaultdict, Counter, namedtuple, OrderedDict
from copy import copy
import datetime
from operator import itemgetter, methodcaller, mul
from itertools import chain, combinations, combinations_with_replacement, groupby, permutations, product
from math import exp
import glob
import numpy as np
from random import shuffle, randint
import primer3
import os
import sys
import time
import tempfile
import subprocess
from subprocess import check_output, CalledProcessError, Popen
from primer3.bindings import designPrimers, setP3Globals
import networkx as nx
from pyfaidx import Fasta

sys.path.append("/raid1/projects/CGDB/models/database_merging/scripts")
from RNAIsoform import RNAIsoform
import designParams

import pdb


class RNAIsoform2(RNAIsoform):
    def __init__(self, CGDB_name, chromosome, strand):
        RNAIsoform.__init__(self, None, None, CGDB_name, chromosome, strand)
        self.is_target = None
        self.amplicons = defaultdict(list)
        self.primer_regions = {"fwd":set(), "rev":set()}  # The sets are of tuples (start, stop), (start, stop), ...

    def setAsTarget(self):
        self.is_target = True

    def setAsNontarget(self):
        self.is_target = False


    def areSequential(self, genomic_pos_1, genomic_pos_2):
        return (abs(self.coords_genome_to_mRNA[genomic_pos_1] - self.coords_genome_to_mRNA[genomic_pos_2]) == 1)


    def hasPositions(self, genomic_positions):
        return all(map(lambda p: p in self.coords_genome_to_mRNA, genomic_positions))


    def areAllSequential(self, genomic_positions):
        return len(genomic_positions) == abs(self.coords_genome_to_mRNA[max(genomic_positions)] - self.coords_genome_to_mRNA[min(genomic_positions)]) + 1


    def getRevPrimerStartingAt(self, primer_genomic_5p_pos, primer_len):
        '''NOTE: method returns the reverse complement of what the actual primer sequence would be.
        The returned primer_genomic_positions are ordered to correspond to the primer sequence 5'->3'.'''
        primer_seq = None
        primer_genomic_positions = None

        primer_mRNA_5p_pos = self.coords_genome_to_mRNA[primer_genomic_5p_pos] 
        if (primer_mRNA_5p_pos < self.length and primer_mRNA_5p_pos >= primer_len-1):
            primer_genomic_positions = tuple(map(lambda mRNA_pos: self.coords_mRNA_to_genome[mRNA_pos], range(primer_mRNA_5p_pos, primer_mRNA_5p_pos-primer_len, -1)))
            assert (len(primer_genomic_positions) == primer_len)
            assert ((self.strand=='+' and primer_genomic_positions[0] > primer_genomic_positions[-1]) or
                (self.strand=='-' and primer_genomic_positions[0] < primer_genomic_positions[-1]))
            primer_seq = self.nuc_seq[primer_mRNA_5p_pos-primer_len+1:primer_mRNA_5p_pos+1]
            assert (len(primer_seq) == primer_len)

        return (primer_seq, primer_genomic_positions)

        
    def getRevPrimerEndingAt(self, primer_genomic_3p_pos, primer_len):
        '''NOTE: method returns the reverse complement of what the actual primer sequence would be.'''
        primer_seq = None
        primer_genomic_positions = None

        primer_mRNA_3p_pos = self.coords_genome_to_mRNA[primer_genomic_3p_pos] 
        if (primer_mRNA_3p_pos <= self.length - primer_len):
            primer_genomic_positions = [primer_genomic_3p_pos]
            for mRNA_pos in range(primer_mRNA_3p_pos+1, primer_mRNA_3p_pos + primer_len):
                genomic_pos = self.coords_mRNA_to_genome[mRNA_pos]
                primer_genomic_positions.append( genomic_pos )
            assert (len(primer_genomic_positions) == primer_len)
            primer_genomic_positions.sort()
            primer_genomic_positions = tuple(primer_genomic_positions)
            primer_seq = self.nuc_seq[primer_mRNA_3p_pos:primer_mRNA_3p_pos+primer_len]
            assert (len(primer_seq) == primer_len)

        return (primer_seq, primer_genomic_positions)

        
    def getFwdPrimerStartingAt(self, primer_genomic_5p_pos, primer_len):
        '''The returned primer_genomic_positions are ordered to correspond to the primer sequence 5'->3'.'''

        primer_seq = None
        primer_genomic_positions = None

        primer_mRNA_5p_pos = self.coords_genome_to_mRNA[primer_genomic_5p_pos]
        if (primer_mRNA_5p_pos + primer_len <= self.length):
            primer_genomic_positions = tuple(map(lambda mRNA_pos: self.coords_mRNA_to_genome[mRNA_pos], range(primer_mRNA_5p_pos, primer_mRNA_5p_pos + primer_len)))
            assert (len(primer_genomic_positions) == primer_len)
            assert ((self.strand=='+' and primer_genomic_positions[0] < primer_genomic_positions[-1]) or
                (self.strand=='-' and primer_genomic_positions[0] > primer_genomic_positions[-1]))
            primer_seq = self.nuc_seq[primer_mRNA_5p_pos:primer_mRNA_5p_pos + primer_len]
            assert (len(primer_seq) == primer_len)
                
        return (primer_seq, primer_genomic_positions)

        
    def getFwdPrimerEndingAt(self, primer_genomic_3p_pos, primer_len):
        primer_seq = None
        primer_genomic_positions = None

        primer_mRNA_3p_pos = self.coords_genome_to_mRNA[primer_genomic_3p_pos]
        if (primer_mRNA_3p_pos - primer_len + 1 >= 0):
            primer_genomic_positions = [primer_genomic_3p_pos]
            for mRNA_pos in range(primer_mRNA_3p_pos - primer_len + 1, primer_mRNA_3p_pos):
                genomic_pos = self.coords_mRNA_to_genome[mRNA_pos]
                primer_genomic_positions.append( genomic_pos )
            assert (len(primer_genomic_positions) == primer_len)
            primer_genomic_positions.sort()
            primer_genomic_positions = tuple(primer_genomic_positions)
            primer_seq = self.nuc_seq[primer_mRNA_3p_pos-primer_len+1:primer_mRNA_3p_pos+1]
            assert (len(primer_seq) == primer_len)
                
        return (primer_seq, primer_genomic_positions)


#    def addAmpliconAnnotation(self, genomic_start, genomic_stop, experiment_number, amplicon_number, ppr_ID):
#        # Get all exonic genomic positions
#        mrna_pos1 = self.coords_genome_to_mRNA[genomic_start]
#        mrna_pos2 = self.coords_genome_to_mRNA[genomic_stop]
#        assert (mrna_pos1 < mrna_pos2)
#        amplicon_genomic_positions = set(self.coords_mRNA_to_genome[mrna_pos1:mrna_pos2+1])
#        print >> sys.stderr, "%s\tExperiment/amplicon = %d/%d,\tamplicon_length=%d,\tppr ID=%s" % \
##          (self.ID, experiment_number, amplicon_number, len(amplicon_genomic_positions), ppr_ID)
#
#        for exon in self.ordered_exons:
#            exon_positions = set(range(exon_start, exon_stop+1))
#            common_positions = exon_positions & amplicon_genomic_positions
#            if (len(common_positions) > 0):
#                start_pos = min(common_positions)
#                stop_pos = max(common_positions)
#                self.amplicons[(experiment_number,amplicon_number)].append( (start_pos, stop_pos, exon.getAsGFF3()) )


    def getAsGFF3(self, parent_gff3_ID):
        mRNA_tuple, exon_tuples = RNAIsoform.getAsGFF3(self, parent_gff3_ID)

        primer_region_tuples = []
        for start, stop, primer_pair_region_ID, fwd_primer_region_ID in self.primer_regions["fwd"]:
            parent_exons = list(filter(lambda x: x[1]<=start and x[2]>=stop, exon_tuples))
            assert (len(parent_exons) == 1)
            primer_region_tuples.append( (self.chromosome, start, stop, self.strand, fwd_primer_region_ID, parent_exons[0]) ) # primer_pair_region_ID, 

        for start, stop, primer_pair_region_ID, rev_primer_region_ID in self.primer_regions["rev"]:
            parent_exons = list(filter(lambda x: x[1]<=start and x[2]>=stop, exon_tuples))
            assert (len(parent_exons) == 1)
            primer_region_tuples.append( (self.chromosome, start, stop, self.strand, rev_primer_region_ID, parent_exons[0]) ) # primer_pair_region_ID, 

        amplicon_tuples = []
        for score_val, amplicon_parts in self.amplicons.items():
            for part_counter, (start, stop, exon_tuple) in enumerate(amplicon_parts,1):
                amplicon_tuples.append( (self.chromosome, start, stop, score_val, self.strand, "amp%d.%d" % (score_val,part_counter), exon_tuple) )

        return (mRNA_tuple, exon_tuples, primer_region_tuples, amplicon_tuples)
