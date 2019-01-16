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


class SetOfPrimers(object):
    def __init__(self, chromosome, OS_ID, PPR_combo, primer_pairs_descriptor, primers):
        self.chromosome = chromosome
        self.OS_ID = OS_ID
        self.PPR_combo = PPR_combo
        self.primer_pairs_descriptor = primer_pairs_descriptor

        # List, where each of local/global isoform group is a tuple of isoform IDs
        self.target_isoform_groups_by_ID = None
        # List of isoform groups, where each isoform group is a tuple of isoform IDs
        self.nontarget_isoform_groups_by_ID = None
        self.other_isoform_groups_by_ID = None
        
        self.fwd_primers = {}
        self.rev_primers = {}
        for erID, primer in primers:
            try:
                if (erID[0] == 'F'):
                    if (erID not in self.fwd_primers):
                        self.fwd_primers[erID] = primer
                    else:
                        assert (self.fwd_primers[erID] == primer)
                else:
                    if (erID not in self.rev_primers):
                        self.rev_primers[erID] = primer
                    else:
                        assert (self.rev_primers[erID] == primer)
            except AssertionError as ae:
                pdb.set_trace()
                print(erID, primer, file=sys.stderr, flush=True)
                

    def setIsoformGroupsFromCategorized(self, target_isoform_groups_by_ID, nontarget_isoform_groups_by_ID, other_isoform_groups_by_ID):
        assert (len(target_isoform_groups_by_ID) > 0), "No target isoform groups"
        self.target_isoform_groups_by_ID = target_isoform_groups_by_ID
        self.nontarget_isoform_groups_by_ID = nontarget_isoform_groups_by_ID
        self.other_isoform_groups_by_ID = other_isoform_groups_by_ID


    def getGlobalTargetIsoformIDs(self):
        assert (self.target_isoform_groups_by_ID != None)
        return set(chain.from_iterable(self.target_isoform_groups_by_ID))


    def getAllFwdPrimers(self):
        return self.fwd_primers.values()

    
    def getAllRevPrimers(self):
        return self.rev_primers.values()


    def getAllPrimers(self):
        return list(self.fwd_primers.values()) + list(self.rev_primers.values())


    def getAllAnnotPrimers(self):
        return list(self.fwd_primers.items()) + list(self.rev_primers.items())


    def isSameAs(self, other_sop):
        my_primers = set(self.getAllAnnotPrimers())
        other_primers = set(other_sop.getAllAnnotPrimers())
        return my_primers == other_primers
