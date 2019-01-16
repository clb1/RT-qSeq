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

import pdb


class GenomicPositionData(object):
    def __init__(self, genomic_position):
        self.genomic_position = genomic_position
        self.isoforms = set()

        # Each is a list of [(primer, primer genomic positions, set of isoforms)]
        self.fwd_primer_groups = []
        self.rev_primer_groups = []


    def hasNoData(self):
        return (len(self.fwd_primer_groups)==0 and len(self.rev_primer_groups)==0)

        
    def getGenomicPosition(self):
        return self.genomic_position


    def getFwdPrimerGroups(self):
        return self.fwd_primer_groups 


    def getRevPrimerGroups(self):
        return self.rev_primer_groups


    def addIsoform(self, new_isoform):
        assert (new_isoform not in self.isoforms)
        self.isoforms.add(new_isoform)


    def getIsoforms(self):
        return self.isoforms
    

    def setPrimerGroups(self):
        '''Record the different fwd and rev primers whose 5' end is the genomic position of self.
        Associate isoforms with specific primers.'''
        min_primer_len = designParams.min_primer_len
        
        fwd_primer_groups = defaultdict(list)
        rev_primer_groups = defaultdict(list)

        fwd_primer_genomic_position_paths = defaultdict(list)
        rev_primer_genomic_position_paths = defaultdict(list)
        
        for isoform in self.isoforms:
            fwd_primer_seq, fwd_primer_genome_positions = isoform.getFwdPrimerStartingAt(self.genomic_position, min_primer_len)
            if (fwd_primer_seq != None):
                assert (len(fwd_primer_genome_positions) == min_primer_len and len(fwd_primer_seq) == min_primer_len)
                fwd_primer_groups[fwd_primer_seq].append(isoform)
                fwd_primer_genomic_position_paths[fwd_primer_seq].append( (isoform, fwd_primer_genome_positions) )
                #fwd_primer_stop_per_isoform[fwd_primer_seq].append( (isoform, fwd_primer_genome_positions[-1]) )

            rev_primer_seq, rev_primer_genome_positions = isoform.getRevPrimerStartingAt(self.genomic_position, min_primer_len)
            if (rev_primer_seq != None):
                assert (len(rev_primer_genome_positions) == min_primer_len and len(rev_primer_seq) == min_primer_len)
                rev_primer_groups[rev_primer_seq].append(isoform)
                rev_primer_genomic_position_paths[rev_primer_seq].append( (isoform,rev_primer_genome_positions) )

        for primer_seq, isoforms_w_primer in fwd_primer_groups.items():
            isoforms_tuple = sorted(isoforms_w_primer, key=methodcaller("getCGDBName"))    # Instantiate an ordering based on the isoforms' names.
            self.fwd_primer_groups.append( (primer_seq, tuple(isoforms_tuple), fwd_primer_genomic_position_paths[primer_seq]) )

        for primer_seq, isoforms_w_primer in rev_primer_groups.items():
            isoforms_tuple = sorted(isoforms_w_primer, key=methodcaller("getCGDBName"))    # Instantiate an ordering based on the isoforms' names.
            self.rev_primer_groups.append( (primer_seq, tuple(isoforms_tuple), rev_primer_genomic_position_paths[primer_seq]) )


    def getPrimerGroups(self, rev_fwd):
        assert (rev_fwd in ["rev", "fwd"])
        return self.rev_primer_groups if (rev_fwd=="rev") else self.fwd_primer_groups
