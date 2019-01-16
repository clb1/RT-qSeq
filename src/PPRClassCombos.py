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


class PPRClassCombos(object):
    def __init__(self, olap_set, transcriptome_ref, transcriptome_vsearch_udb, antitargets_fasta_name): # genome_ref transcriptomic_MFE,
        self.olap_set = olap_set
        olap_set.checkDirectionality()

        self.oligo_thermo = OligoThermodynamics()
        self.transcriptome_ref = transcriptome_ref
        self.transcriptome_vsearch_udb = transcriptome_vsearch_udb
        self.antitargets_fasta_name = antitargets_fasta_name

        self.context_strand = olap_set.getStrand()
        self.target_isoform_groups = olap_set.getTargetIsoformGroups("primerable")
        self.target_isoform_groups_IDs = set(map(lambda x: tuple(map(lambda y: y.getCGDBName(), x)), self.target_isoform_groups))
        self.target_isoforms = set(chain.from_iterable(self.target_isoform_groups))
        self.target_isoforms_IDs = set(chain.from_iterable(self.target_isoform_groups_IDs))

        self.all_primer_pair_regions = olap_set.getEquivPrimerRegionPairs(only_with_default_primer_pair=True)
        self.ppr_classes = []
        self.pprs_per_class = []

        self.complete_TIG_IDs_w_nonlocal_members = self.olap_set.getGlobalTIGsOfLocalTIGs()

        # The number of *isoforms* that are not in a target isoform group but that are primed.
        # Among these isoforms can be target isoforms that are in a non-target isoform group
        # (and so, not distinguishable from non-target isoforms, but primed nonetheless).
        #
        # Redundant isoforms counts, accumulated later, refer to the number of extra times that target
        # isoforms are primed as part of a target isoform group (ie redundantly informative)
        self.ppr_classes_num_offtarget = []

        assert (len(self.target_isoform_groups)>0), "No primerable target isoform groups"

        self.transcriptomic_MFE = None # transcriptomic_MFE

        self.lowest_num_offtarget = None # These four are based on knowing the pair of equivalence regions and primer start positions
        self.lowest_num_redundant = None
        self.lowest_num_unwanted = None # = self.lowest_num_offtarget + self.lowest_num_redundant
        self.lowest_expected_num_unwanted_from_primers = None # Based on global interactions among primers in a complete solution
        self.smallest_Npart = None      # The number of separate RT-qSeq experiments needed to prime all of the TIGs

        print >> sys.stderr, "Finding conflicting pprs...",
        self.conflicting_pprs = olap_set.getConflictingPrimerPairRegions(include_ppr_IDs=True)
        print >> sys.stderr, "done. Size = %d" % len(self.conflicting_pprs)
        self.all_primer_pair_regions_lookup = {}

        for ppr in self.all_primer_pair_regions:
            fwd_and_rev_region_IDs_tup = ppr.getFwdRevIDs()
            self.all_primer_pair_regions_lookup[ fwd_and_rev_region_IDs_tup ] = ppr

        # Place the PPR classes in a ranked order based on 1) having a unique TIG, 2) fewest number of offtargets IGs and then 3) larger TIGs
        L = {}
        P = defaultdict(list)
        self.ig_occur_counter = Counter()

        for ppr in filter(lambda x: x.couldAmplifyAnyTargetGroups(self.target_isoform_groups), self.all_primer_pair_regions):
            igs = ppr.getIsoformGroupsDualPrimed()
            igs_in_set = set(igs)
            num_tigs = sum(map(lambda x: x in self.target_isoform_groups, igs))
            assert (num_tigs > 0), "Adding PPR that doesn't have a target isoform group"
            num_offtarget_isoforms = len(list(chain.from_iterable(igs_in_set - self.target_isoform_groups)))
            L[(num_offtarget_isoforms, -num_tigs, igs)] = igs_in_set
            P[igs].append(ppr)
            self.ig_occur_counter.update(igs_in_set)

        # Get the TIGs that only occur in one PPR class and place those classes first
        uniquely_occur_TIGs = set(map(itemgetter(0), filter(lambda t: t[0] in self.target_isoform_groups and t[1]==1, self.ig_occur_counter.most_common())))
        for key_tup in sorted(L.keys(), key=itemgetter(0,1)):
            if (not L[key_tup].isdisjoint(uniquely_occur_TIGs)):
                self.ppr_classes.append( L[key_tup] )
                self.pprs_per_class.append( P[key_tup[2]] )
                self.ppr_classes_num_offtarget.append( key_tup[0] )

        # Now the other PPR classes
        for key_tup in sorted(L.keys(), key=itemgetter(0,1)):
            if (L[key_tup].isdisjoint(uniquely_occur_TIGs)):
                self.ppr_classes.append( L[key_tup] )
                self.pprs_per_class.append( P[key_tup[2]] )
                self.ppr_classes_num_offtarget.append( key_tup[0] )

        assert (len(self.ppr_classes) > 0)

        self.expire_time = None
        self.time_expired = False

        self.tempdir = None
        self.local_transcriptome_fasta = None


    def setTemporaryDirectory(self, tempdir):
        self.tempdir = tempdir


    def setLocalTranscriptomeFasta(self):
        assert (self.local_transcriptome_fasta == None)
        assert (self.tempdir != None)
        
        self.local_transcriptome_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=self.tempdir, delete=True)
        isoforms_seqs = self.olap_set.getAllOverlappingIsoformsSeqs()
        for isoform_ID, isoform_seq in isoforms_seqs.items():
            self.local_transcriptome_fasta.write(">%s\n%s\n" % (isoform_ID, isoform_seq))
        self.local_transcriptome_fasta.flush()


    def unsetLocalTranscriptomeFasta(self):
        assert (self.local_transcriptome_fasta != None)
        self.local_transcriptome_fasta.close()
        self.local_transcriptome_fasta = None


    def getAntitargetsFastaFileName(self):
        return self.antitargets_fasta_name


    def setExpireTime(self, begin_time, minutes_from_now):
        self.expire_time = time.time() + minutes_from_now * 60
        expire_datetime = begin_time + datetime.timedelta(minutes = minutes_from_now)
        expire_datetime_str = expire_datetime.strftime("%Y-%m-%d %H:%M:%S")
        print("Expire time set to %s" % expire_datetime_str, file=sys.stderr, flush=True)


    def setPPRsMandatoryInSolution(self):
        '''Each PPR Class that contains the only instance of a target isoform group must necessarily be part of a complete solution.
        This method identifies such cases and initializes the solution with them, which has the added effect of reducing the search space
        for the remainder of the solution.'''
        S = np.zeros((len(self.ppr_classes),), dtype=np.uint16)
        num_offtarget = 0
        should_have_unique_TIG = True # To confirm that all PPR classes with a unique TIG occur before those that do not
        last_mandatory_col = -1
        counter_TIG_isoforms_in_mandatory = Counter()
        TIGs_in_mandatory = set()
        for c, igs_in_set in enumerate(self.ppr_classes):
            if (any(map(lambda x: x in self.target_isoform_groups and self.ig_occur_counter[x]==1, igs_in_set))):
                assert (should_have_unique_TIG), "Incorrect ordering of TIGs"
                S[c] = 1
                num_offtarget += self.ppr_classes_num_offtarget[c]
                counter_TIG_isoforms_in_mandatory.update(chain.from_iterable(self.target_isoform_groups & igs_in_set))
                TIGs_in_mandatory.update(self.target_isoform_groups & igs_in_set)
                last_mandatory_col = c
            else:
                should_have_unique_TIG = False
            
        remain_col_soln_space = []
        for col in range(last_mandatory_col+1, len(self.ppr_classes)):
            # Logic: Does the ppr class have a TIG that isn't already accounted for in the mandatory set?
            if (len(self.ppr_classes[col] - TIGs_in_mandatory & self.target_isoform_groups) > 0): 
                remain_col_soln_space.append(col)
        num_redundant = sum(map(lambda x: x-1, counter_TIG_isoforms_in_mandatory.values()))
        remaining_tigs = self.target_isoform_groups - TIGs_in_mandatory

        assert (len(remain_col_soln_space) != 0), "No remaining solution space"
        next_curr_col = remain_col_soln_space[0]
        remain_col_soln_space = [] if (len(remain_col_soln_space)==1) else remain_col_soln_space[1:]

        return (S, next_curr_col, remain_col_soln_space, num_offtarget, num_redundant, remaining_tigs)


    def printPPRClassPatterns(self):
        for i, ppr_class_pattern in enumerate(self.ppr_classes):
            print("\n== %d ==" % i, file=sys.stderr, flush=True)
            for tup in ppr_class_pattern:
                print("\t%s" % ", ".join(list(map(lambda x: x.getCGDBName(), tup))), file=sys.stderr, flush=True)
            print("\nPPRs: %s" % ", ".join(list(map(lambda ppr: ppr.getID(), self.pprs_per_class[i]))), file=sys.stderr, flush=True)


    def findSolution(self):
        assert (self.tempdir != None and os.path.exists(self.tempdir)), "Temporary directory not set or does not exist"
        assert (len(self.ppr_classes) > 0), "No PPR classes with which to find solution"
        self.lowest_num_offtarget = 1e10
        self.lowest_num_redundant = 1e10
        self.lowest_num_unwanted = 1e10
        self.lowest_expected_num_unwanted_from_primers = 1e10
        self.smallest_Npart = 1e10

        try:
            S, start_col, remain_col_soln_space, num_offtarget, num_redundant, remaining_tigs = self.setPPRsMandatoryInSolution()
        except AssertionError as ae:
            if (str(ae) == "No remaining solution space"):
                S = np.zeros((len(self.ppr_classes),), dtype=np.uint16)
                start_col = 0
                remain_col_soln_space = list(range(1,len(self.ppr_classes)))
                num_offtarget = 0
                num_redundant = 0
                remaining_tigs = self.target_isoform_groups
            else:
                raise ae

        self.printPPRClassPatterns()

        #try:
        solutions = self.recurseOnPPRClasses(S, len(self.target_isoform_groups), start_col,
                                             remain_col_soln_space, num_offtarget, num_redundant, remaining_tigs)

        message = "Search timed out. " if (self.time_expired) else ""
        message += "Found a total of %d PPRs combinations solutions" % len(solutions)

        #timed_out = False
        #except TimeoutException as te:
        #    message = str(te)
        #    timed_out = True
        #    if (te.hasSolutions()):
        #        solutions = te.getSolutions()
        #    else:
        #        solutions = []

        solution = None if (len(solutions)==0) else solutions[-1]
        return (solution, self.time_expired, message) # timed_out


    def recurseOnPPRClasses(self, S, max_soln_size, curr_col, remain_col_soln_space, num_offtarget, num_redundant, remaining_tigs):
        Npart_solutions = []

        remain_tigs_covered_by_this_col = set(filter(lambda tig: tig in self.ppr_classes[curr_col], remaining_tigs))
        tigs_already_covered = self.target_isoform_groups - remaining_tigs
        new_num_redundant = num_redundant + len(list(chain.from_iterable(filter(lambda ig: ig in self.ppr_classes[curr_col], tigs_already_covered))))
        new_remaining_tigs = remaining_tigs - remain_tigs_covered_by_this_col
        new_num_offtarget = num_offtarget + self.ppr_classes_num_offtarget[curr_col]
        new_num_unwanted = new_num_offtarget + new_num_redundant

        sum_S = sum(S)

        if (new_num_unwanted <= self.lowest_num_unwanted and (self.lowest_num_unwanted > 0 or self.lowest_expected_num_unwanted_from_primers > 0)):
            if (len(new_remaining_tigs) == 0):
                assert(sum_S < max_soln_size), "Candidate solution has more PPRs than there are TIGs"
                S[curr_col] = 1
                print("\n\t*Searching for primers solution for PPR class combination*", file=sys.stderr, flush=True)
                globally_unverified_Npart_solutions = self.findPDSForPPRClassComboGreedyPath(S)
                if (globally_unverified_Npart_solutions != None and len(globally_unverified_Npart_solutions) > 0):
                    assert (len(Npart_solutions) <= self.smallest_Npart), "Current Npart is larger than the smallest Npart found so far. Should not happen."
                    Npart_solution, expected_num_unwanted_from_primers = verifyNpartSolutionsGlobally(globally_unverified_Npart_solutions, self.transcriptome_vsearch_udb,
                                                                                                      self.transcriptome_ref, self.oligo_thermo, self.complete_TIG_IDs_w_nonlocal_members,
                                                                                                      self.tempdir)

                    # Prevent a solution with new_num_unwanted==self.lowest_num_unwanted and expected_num_unwanted_from_primers==self.lowest_expected_num_unwanted_from_primers
                    # but with more PPRs (high sum_S) from being added as a better solution (eg OS46, chr7)
                    if (Npart_solution != None and expected_num_unwanted_from_primers <= self.lowest_expected_num_unwanted_from_primers):
                        Npart_solutions.append( (new_num_unwanted, new_num_offtarget, new_num_redundant, expected_num_unwanted_from_primers, len(Npart_solution), Npart_solution) )
                        print("\t*NEW ( Num PPRs =", sum_S+1, ", Npart =", len(Npart_solution), ")\t", Npart_solutions[-1][0:-1], file=sys.stderr, flush=True)

                        self.lowest_num_offtarget = min(self.lowest_num_offtarget, new_num_offtarget)
                        self.lowest_num_redundant = min(self.lowest_num_redundant, new_num_redundant)
                        self.lowest_num_unwanted = min(self.lowest_num_unwanted, new_num_unwanted)
                        self.lowest_expected_num_unwanted_from_primers = min(self.lowest_expected_num_unwanted_from_primers, expected_num_unwanted_from_primers)
                        self.smallest_Npart = len(Npart_solution) if (self.smallest_Npart > len(Npart_solution)) else self.smallest_Npart
                    else: 
                        if (expected_num_unwanted_from_primers <= self.lowest_expected_num_unwanted_from_primers):
                            print("\t*Solution found, but failed to verify globally at Npart level", file=sys.stderr, flush=True)
                        else:
                            print("\t*Solution found, but failed to verify globally at Npart level because high expected number unwanted from primers", file=sys.stderr, flush=True)
                        #Npart_solution = verifyNpartSolutionsGlobally(globally_unverified_Npart_solutions, transcriptome_vsearch_udb, transcriptome_ref,
                        #                                              oligo_thermo, self.complete_TIG_IDs_w_nonlocal_members, self.tempdir)

                    if (time.time() >= self.expire_time):
                        self.time_expired = True
                        #raise TimeoutException("Solution search timed out with %d Npart_solutions." % len(Npart_solutions), Npart_solutions)
                else:
                    print("\t*No primers solution of PPR class combination found", file=sys.stderr, flush=True)
                S[curr_col] = 0

            elif (sum_S < max_soln_size): #curr_col+1 < len(self.ppr_classes) and 
                # Get the columns that have remaining TIGs and that won't result in a higher number of unwanted primerings (offtarget and redundant)
                new_remain_col_soln_space = []
                limit_additional_num_unwanted = self.lowest_num_unwanted - new_num_unwanted
                tigs_already_covered_incl_curr_col = tigs_already_covered | remain_tigs_covered_by_this_col
                for remain_col in filter(lambda c: len(self.ppr_classes[c] & new_remaining_tigs) > 0, remain_col_soln_space): # For remaining columns w/ partial soln...
                    assert (S[remain_col] == 0), "S[remain_col] != 0"
                    col_num_unwanted = self.ppr_classes_num_offtarget[remain_col] + \
                      len(list(chain.from_iterable(tigs_already_covered_incl_curr_col & self.ppr_classes[remain_col])))
                    if (col_num_unwanted <= limit_additional_num_unwanted):
                        new_remain_col_soln_space.append(remain_col)

                # Do the columns together have the TIGs needed for a solution?
                soln_space_tigs = set()
                for c in new_remain_col_soln_space:
                    soln_space_tigs.update(self.ppr_classes[c])

                if (new_remaining_tigs.issubset(soln_space_tigs)):
                    S[curr_col] = 1
                    next_col = new_remain_col_soln_space[0]
                    next_soln_space = [] if (len(new_remain_col_soln_space)==1) else new_remain_col_soln_space[1:]
                    assert (S[next_col] == 0), "S[next_col] != 0"
                    other_Npart_solutions = self.recurseOnPPRClasses(S[:], max_soln_size, next_col, next_soln_space, new_num_offtarget, new_num_redundant, new_remaining_tigs)
                    if (len(other_Npart_solutions) > 0):
                        Npart_solutions.extend(other_Npart_solutions)
                    S[curr_col] = 0

        # Also want Npart_solutions that do not include this column's TIG
        # The candidate ppr_classes are ordered by increasing number of offtargets, so don't proceed if no better solution than those already found can be found
        # TODO: Incorporate?: num_offtarget + self.ppr_classes_num_offtarget[curr_col+1] <= self.lowest_num_offtarget and curr_col+1 < len(self.ppr_classes) and 
        if (sum_S < max_soln_size and not self.time_expired): # TODO: should have some logic about continuing searching if non-optimal solution not yet found
            # Get the columns that have remaining TIGs and that won't result in a higher number of unwanted primerings (offtarget and redundant)
            new_remain_col_soln_space = []
            limit_additional_num_unwanted = self.lowest_num_unwanted - (num_offtarget + num_redundant)
            assert (limit_additional_num_unwanted >= 0)

            for remain_col in filter(lambda c: len(self.ppr_classes[c] & remaining_tigs) > 0, remain_col_soln_space): # For remaining columns w/ partial soln...
                assert (S[remain_col] == 0), "S[remain_col] != 0"
                num_unwanted = self.ppr_classes_num_offtarget[remain_col] + \
                  len(list(chain.from_iterable(tigs_already_covered & self.ppr_classes[remain_col])))
                if (num_unwanted <= limit_additional_num_unwanted):
                    new_remain_col_soln_space.append(remain_col)

            # Do the columns' together have the TIGs needed for a solution?
            if (len(new_remain_col_soln_space) > 0):
                soln_space_tigs = set()
                for c in new_remain_col_soln_space:
                    soln_space_tigs.update(self.ppr_classes[c])

                if (remaining_tigs.issubset(soln_space_tigs)):
                    S[curr_col] = 0
                    next_col = new_remain_col_soln_space[0]
                    next_soln_space = [] if (len(new_remain_col_soln_space)==1) else new_remain_col_soln_space[1:]
                    assert (S[next_col] == 0), "S[next_col] != 0"
                    other_Npart_solutions = self.recurseOnPPRClasses(S[:], max_soln_size, next_col,
                                                                     next_soln_space, num_offtarget, num_redundant, remaining_tigs)
                    if (len(other_Npart_solutions) > 0):
                        Npart_solutions.extend(other_Npart_solutions)

        # TODO: still do this?
        Npart_solutions = list(filter(lambda x: x[0] <= self.lowest_num_unwanted or x[4] <= self.smallest_Npart, Npart_solutions))
        
        return Npart_solutions


    def findPDSForPPRClassComboGreedyPath(self, pprclass_pattern):
        '''Tries to return an Npart RT-qSeq experiment that amplifies all TIGs in pprclass_pattern.'''
        # Form M, a conceptual submatrix of all pprclass patterns. 
        # A "column" is a pprclass pattern and rows (the number of which are different for each pprclass) are the pprs that fit that pattern.
        # The "columns" of M, indicated by a column number, are the different ppr_classes, and the elements in M[column_number] are the pprs for that class.
        Npart_soln = []

        init_target_num_equiv_regions = {2:2, 3:4, 4:5, 5:7, 6:9}

        # Construct a OrderedDiGraph, where nodes are the pprs and an edge exists for pprs that do not conflict.
        # OrderedDiGraph is used so that cliques with lowest rank PPRs in columns of M are produced first
        MG = nx.OrderedDiGraph()

        M, M_tigs = {}, {}
        M_global_indices = []
        ppr_my_class = {}
        equiv_region_ID_counts = defaultdict(int)
        all_ppr_indices_w_commonness = []

        # First determine commonness of the equiv regions
        for i in range(len(pprclass_pattern)):
            if (pprclass_pattern[i] == 1):
                for ppr in self.pprs_per_class[i]:
                    ppr_my_class[ppr.getID()] = i
                    fwd_region_ID, rev_region_ID = ppr.getFwdRevIDs()
                    equiv_region_ID_counts[fwd_region_ID] += 1
                    equiv_region_ID_counts[rev_region_ID] += 1

        # Build graph structure to bias search towards using PPRs with the more common equiv regions
        ppr_to_node = {}
        for i in range(len(pprclass_pattern)):
            if (pprclass_pattern[i] == 1):
                # Sort the PPRs by commonness of fwd/rev_region_IDs, to aid efficiency searching below
                sorted_pprs = []
                for ppr in self.pprs_per_class[i]:
                    fwd_region_ID, rev_region_ID = ppr.getFwdRevIDs()
                    sorted_pprs.append( (ppr, equiv_region_ID_counts[fwd_region_ID] + equiv_region_ID_counts[rev_region_ID]) )
                sorted_pprs.sort(key=itemgetter(1), reverse=True)

                for m, (ppr,id_count) in enumerate(sorted_pprs): 
                    all_ppr_indices_w_commonness.append( ((i,m), id_count) )
                    MG.add_node( (i,m), FwdRev_IDs=set(ppr.getFwdRevIDs()) )
                    ppr_to_node[ppr] = (i,m)

                M[i] = list(map(itemgetter(0), sorted_pprs))
                M_tigs[i] = self.ppr_classes[i]
                M_global_indices.append(i)

                #print >> sys.stderr, "len(M[%d]) = %d" % (i, len(M[i]))
                #C = Counter(map(lambda ppr: ppr.getNumOfftargetIsoforms(self.target_isoform_groups), self.pprs_per_class[i]))
                #offtargets_str = ", ".join(map(lambda t: "%d: %d" % t, C.most_common()))
                #print >> sys.stderr, "\tpprclass %d : %d pprs;\tOfftarget: %s" % (i, len(self.pprs_per_class[i]), offtargets_str)

        # Start searches biased to the PPRs with the most commonly used equiv_regions
        all_ppr_indices_w_commonness.sort(key=itemgetter(1), reverse=True)
        for M_index in map(itemgetter(0), all_ppr_indices_w_commonness):
            MG.add_edge( 'B', M_index )

        ppr_lookup = self.olap_set.getEquivPrimerRegionPairsByName()

        # Add the edges to the graph
        for col1, col2 in combinations(M_global_indices, 2):
            for ppr1_index, ppr1 in enumerate(M[col1]):
                for ppr2_index, ppr2 in enumerate(M[col2]):
                    do_conflict, isoform_w_conflict = ppr1.doesConflictOnAnIsoformWith(ppr2, self.target_isoforms)
                    if (not do_conflict):
                        MG.add_edge( (col1,ppr1_index), (col2,ppr2_index) )

        # No solution possible if there are no edges (ie all pprs overlap)
        if (MG.size() == 0):
            return Npart_soln

        # Ideally, want to search for large cliques, but in practice the search time can be too long. So, use a manageable value here as a compromise.
        print("\t\tBeginning Npart search", file=sys.stderr, flush=True)
        start_clique_size = min(nx.algorithms.dag.dag_longest_path_length(MG), 6)
        assert (start_clique_size > 0), "No path in graph MG"
        remaining_Mcols = set(M.keys())
        kept_found_clique = True
        target_num_equiv_regions = None
        #blacklist_pprs = set() 
        while (len(remaining_Mcols) > 0 and kept_found_clique and len(Npart_soln) < self.smallest_Npart):
            print("\t\tMcols remaining = %d" % len(remaining_Mcols), file=sys.stderr, flush=True)
            for clique_size in range(start_clique_size, 0, -1):
                kept_found_clique = False

                assert (clique_size <= len(remaining_Mcols)), "clique size larger than the number of available columns"
                if (clique_size == 1):
                    ppr_index = 0
                    a_col = list(remaining_Mcols)[0]
                    ppr = M[a_col][0]

                    all_pps = []
                    chrom, OS_ID = self.olap_set.getChromosomeAndID()
                    for sop in ppr.getDefaultPrimerPairs():
                        all_pps.append(sop)

                    Npart_soln.append( (set([a_col]), (ppr,), all_pps) )
                    remaining_Mcols.remove(a_col)
                    start_clique_size = 0 if (len(remaining_Mcols)==0) else 1
                    kept_found_clique = True
                    break
                else:
                    # target_num_equiv_regions should only be reset if the clique size changes. Otherwise, the previous (productive) value should be used.
                    #target_num_equiv_regions =  2*clique_size - 2 if (target_num_equiv_regions == None) else target_num_equiv_regions
                    target_num_equiv_regions = init_target_num_equiv_regions[clique_size] if (target_num_equiv_regions == None) else target_num_equiv_regions
                    while (target_num_equiv_regions < 2*clique_size):
                        target_num_equiv_regions += 1
                        print("\t\t\tEvaluating PPR cliques of size %s at efficiency level %d..." % (clique_size, target_num_equiv_regions), file=sys.stderr, flush=True)
                        ppr_cliques_generator = self.findPPRCliquesAtEfficiencyLevel(MG, 'B', clique_size, target_num_equiv_regions)

                        try:
                            #ppr_combinations_w_no_pp = []
                            # Eval PPR cliques in lowest rank priority order
                            begin_time = datetime.datetime.now()
                            for ppr_clique in ppr_cliques_generator:
                                ppr_combo = tuple(map(lambda x: M[x[0]][x[1]], ppr_clique))
                                ppr_IDs_clique = list(map(lambda x:x.getID(), ppr_combo))
                                print("\t\t\t%s" % " ".join(ppr_IDs_clique), file=sys.stderr, flush=True)
                                
                                #if (any(map(lambda x: x <= set(ppr_combo), ppr_combinations_w_no_pp))): # any(map(lambda x: x in blacklist_pprs, ppr_combo)) or 
                                #    print >> sys.stderr, "\t\t\tskipped"
                                #    continue

                                all_ppr_FwdRev_IDs_tuples = set(map(lambda ppr: ppr.getFwdRevIDs(), ppr_combo))
                                possible_implicit_pprs_as_tuples = set(product(map(itemgetter(0), all_ppr_FwdRev_IDs_tuples), map(itemgetter(1), all_ppr_FwdRev_IDs_tuples))) - \
                                                                   all_ppr_FwdRev_IDs_tuples
                                implicit_pprs = list(filter(lambda x: x in ppr_lookup, map(lambda s: "%s%s" % s, possible_implicit_pprs_as_tuples)))
            
                                implicit_ppr_conflict = any(map(lambda x: ppr_lookup[x[0]].doesConflictOnAnIsoformWith(ppr_lookup[x[1]], self.target_isoforms)[0],
                                                                combinations(implicit_pprs,2)))

                                explicit_implicit_ppr_conflict = any(map(lambda x: ppr_lookup[x[0]].doesConflictOnAnIsoformWith(ppr_lookup[x[1]], self.target_isoforms)[0],
                                                                         product(map(lambda ppr: ppr.getID(), ppr_combo), implicit_pprs)))
                                        
                                if (not (implicit_ppr_conflict or explicit_implicit_ppr_conflict or
                                         any(map(lambda x: x in ppr_my_class, implicit_pprs)))): # Last clause means "if didn't implicitly create any PPRs"
                                    ppr_combo_primer_sets, this_ppr_combination_w_no_pp = self.findCandidatePrimerSetsForPPRCombo(ppr_combo)
                                    #ppr_combinations_w_no_pp.extend( this_ppr_combinations_w_no_pp )
                                    #for bad_combo in this_ppr_combinations_w_no_pp:
                                    #    print >> sys.stderr, "\t\t\tDisallowing %s" % (" ".join(map(lambda x:x.getID(), bad_combo)),)

                                    if (len(ppr_combo_primer_sets)>0):
                                        kept_found_clique = True  # TODO: Rename this variable
                                        col_combo = set(map(itemgetter(0), ppr_clique))
                                        Npart_soln.append( (col_combo, ppr_combo, ppr_combo_primer_sets) )   # TODO: is the solution for a clique being cached   TODO
                                        remaining_Mcols -= col_combo
                                        start_clique_size, target_num_equiv_regions = (clique_size,target_num_equiv_regions-1) if (clique_size <= len(remaining_Mcols)) \
                                                                                      else (len(remaining_Mcols),None)

                                        # Remove from MG those nodes that are from any of columns that provided the solution PPRs
                                        MG.remove_nodes_from(list(filter(lambda x:x[0] in col_combo, MG.nodes())))

                                        break
                                    else:
                                        nodes = list(map(lambda ppr: ppr_to_node[ppr], this_ppr_combination_w_no_pp))
                                        nodes.sort(key=itemgetter(0))
                                        removed_an_edge = False
                                        for n1, n2 in zip(nodes[0:-1], nodes[1:]):
                                            if (MG.has_edge(n1,n2)):
                                                MG.remove_edge(n1,n2)
                                                removed_an_edge = True
                                        if (not removed_an_edge):
                                            nodes.sort(key=itemgetter(1), reverse=True)
                                            MG.remove_node( nodes[0] )

                                        ppr_cliques_generator = self.findPPRCliquesAtEfficiencyLevel(MG, 'B', clique_size, target_num_equiv_regions)
                                else:
                                    #print >> sys.stderr, "\t\t\timplicit ppr conflict = %s" % str(implicit_ppr_conflict)
                                    #print >> sys.stderr, "\t\t\texplicit-implicit ppr conflict = %s" % str(explicit_implicit_ppr_conflict)
                                    #implicitly_created_PPR = any(map(lambda x: ppr_my_class.has_key(x), implicit_pprs))
                                    #print >> sys.stderr, "\t\t\timplicitly created PPR = %s" % str(implicitly_created_PPR)
                                    equiv_region_IDs_for_implicit = set(chain.from_iterable(map(lambda frID: ppr_lookup[frID].getFwdRevIDs(), implicit_pprs)))
                                    nodes = []
                                    for ppr in map(lambda x: M[x[0]][x[1]], ppr_clique):
                                        fwd_ID, rev_ID = ppr.getFwdRevIDs()
                                        if (fwd_ID in equiv_region_IDs_for_implicit or rev_ID in equiv_region_IDs_for_implicit):
                                            nodes.append( ppr_to_node[ppr] )

                                    nodes.sort(key=itemgetter(0))
                                    removed_an_edge = False
                                    for n1, n2 in zip(nodes[0:-1], nodes[1:]):
                                        if (MG.has_edge(n1,n2)):
                                            MG.remove_edge(n1,n2)
                                            removed_an_edge = True
                                    if (not removed_an_edge):
                                        nodes.sort(key=itemgetter(1), reverse=True)
                                        MG.remove_node( nodes[0] )

                                    ppr_cliques_generator = self.findPPRCliquesAtEfficiencyLevel(MG, 'B', clique_size, target_num_equiv_regions)
                                    #ppr_combinations_w_no_pp.extend( bad_combo )
                                    #print >> sys.stderr, "\t\t\tDisallowing %s" % (" ".join(map(lambda x:x.getID(), bad_combo)),)

                            if (kept_found_clique):
                                break

                        except StopIteration as si:
                            end_time = datetime.datetime.now()
                            time_delta = end_time - begin_time
                            elapsed_seconds = time_delta.total_seconds() 
                            if (elapsed_seconds > 600): # 600 seconds = 10 minutes
                                print("Dropping to smaller clique size", file=sys.stderr, flush=True)
                                target_num_equiv_regions = None
                                break

                    if (kept_found_clique):
                        break
                    else:
                        target_num_equiv_regions = None

        if (len(remaining_Mcols) > 0):
            assert (len(Npart_soln) >= self.smallest_Npart), "Could not find combinations of PPRs that provided a smaller Npart complete solution"
            Npart_soln = None
        else:
            print("\t\tNpart =", len(Npart_soln), "  Parts' sizes = %s" % ", ".join(list(map(lambda x: str(len(x[1])), Npart_soln))), file=sys.stderr, flush=True)
            Npart_soln = list(map(itemgetter(1,2), Npart_soln))

        return Npart_soln


    def findPPRCliquesAtEfficiencyLevel(self, DG, node, path_len, target_num_equiv_regions, path_equiv_regions=set()):
        '''Efficiency is the number of different F#/R# equiv regions in a PPR clique'''
        if (path_len > 0):
            for neighbor in DG.neighbors(node):
                equiv_regions_incl_neighbor = path_equiv_regions | DG.node[neighbor]['FwdRev_IDs']
                if (len(equiv_regions_incl_neighbor) <= target_num_equiv_regions): #  and DG.out_degree(neighbor)>0
                    for path in self.findPPRCliquesAtEfficiencyLevel(DG, neighbor, path_len-1, target_num_equiv_regions, equiv_regions_incl_neighbor):
                        if (all(map(lambda p: DG.has_edge(node,p), path))):
                            if (node == 'B'):
                                yield path
                            else:
                                yield [node]+path

        elif (len(path_equiv_regions | DG.node[node]['FwdRev_IDs']) == target_num_equiv_regions):
            yield [node]
                

    def findCandidatePrimerSetsForPPRCombo(self, ppr_combo):
        '''Select all candidate primer pairs that can potentially be used for each Primer Pair Region (PPR).
        Paired Primer Regions that share a Forward and/or Reverse primer region need to have the same primer in the common region.
        Paired Primer Regions whose Foward and Reverse primer regions overlap must have primers that do not create overlapping amplicons.
        Paired Primer Regions with no such constraints use their default primer pair.
        '''
        ppr_combo_primer_sets = []
        ppr_combination_w_no_pp = None

        # Get the connected components of PPRs in the ppr_combo, where PPRs in the combo share some Fwd/Rev equiv primer regions and so
        # the exact primer position is constrained by the requirements of more than one PPR.

        # Build graph connecting all forward or reverse equiv primer regions that are in the same PPR or that
        # overlap reverse and forward primer regions in two different PPRs
        equiv_regions_lookup = {}
        eprG = nx.Graph() # "equiv primer regions Graph"
        for ppr in ppr_combo:
            fwd_equiv_region = ppr.getFwdEquivRegion()
            rev_equiv_region = ppr.getRevEquivRegion()
            fwd_region_ID = fwd_equiv_region.getID()
            rev_region_ID = rev_equiv_region.getID()
            equiv_regions_lookup[fwd_region_ID] = fwd_equiv_region
            equiv_regions_lookup[rev_region_ID] = rev_equiv_region
            eprG.add_edge(fwd_region_ID, rev_region_ID, ppr=ppr)

        # Find the connected fwd/rev equiv regions that directly or indirectly constrain the primers that may be used for each.
        CC = list(map(lambda c: (c, len(c)), nx.connected_components(eprG)))
        CC.sort(key=itemgetter(1), reverse=True)
        CC = list(map(itemgetter(0), CC))

        # Find the singleton primers that can potentially be used for each of the connected fwd/rev equiv regions
        ppr_per_cc = []
        expected_isoforms_per_primer_seq = defaultdict(set)
        all_equiv_region_primers = {}
        for cc in CC:
            if (len(cc)==2):  # Case '== 2' corresponds to a single PPR, for which default primers are used below
                erID1, erID2 = tuple(cc)
                ppr = eprG[erID1][erID2]['ppr']
                assert (ppr != None)
                ppr_per_cc.append( set([ppr]) )
            else:
                # Determine the primers that can be considered for each equiv primer region, given constraints.
                equiv_region_primers, this_cc_ppr_combo, all_edges_possible = self.constrainPrimersPerEquivRegion(eprG, cc, equiv_regions_lookup)
                if (all_edges_possible):
                    all_equiv_region_primers.update(equiv_region_primers)
                    ppr_per_cc.append( this_cc_ppr_combo )

                    # Associate expected isoforms to each primer
                    for region_ID, singleton_primers in equiv_region_primers.items():
                        equiv_region_isoform_IDs = equiv_regions_lookup[region_ID].getIsoformsIDs()
                        for primer in singleton_primers:
                            expected_isoforms_per_primer_seq[primer.seq].update(equiv_region_isoform_IDs)
                else:
                    ppr_combination_w_no_pp = this_cc_ppr_combo #.append( this_cc_ppr_combo )
                    break

        # If each connected component of fwd/rev equiv regions has possible solutions,
        # then for each PPR get the primer pairs that can distinguishably amplify intended products
        filtered_equiv_region_primers = None
        primer_num_nontarget_hybs = None
        if (len(ppr_per_cc) == len(CC)): 
            any_cc_is_multi_ppr = any(map(lambda x: len(x)>1, ppr_per_cc))
            if (any_cc_is_multi_ppr):
                ppr_combo_unique_target_isoform_IDs = self.target_isoforms_IDs & set.union(*list(map(methodcaller("getIsoformGroupsDualPrimedIsoformsID"), ppr_combo)))
                filtered_equiv_region_primers, primer_num_nontarget_hybs = \
                    scanIndivPrimersAgainstLocalTranscriptome(self.oligo_thermo, expected_isoforms_per_primer_seq, all_equiv_region_primers, self.antitargets_fasta_name, 
                                                              self.local_transcriptome_fasta, ppr_combo_unique_target_isoform_IDs, self.tempdir)

            if (filtered_equiv_region_primers != None or not any_cc_is_multi_ppr):
                ppr_combo_primer_sets, most_freq_incompat_pprs = self.createPPRComboPrimerSets(ppr_per_cc, filtered_equiv_region_primers, primer_num_nontarget_hybs)
                if (most_freq_incompat_pprs != None):
                    ppr_combination_w_no_pp = most_freq_incompat_pprs
                #    ppr_combinations_w_no_pp.append(most_freq_incompat_pprs)
            else:
                ppr_combination_w_no_pp = set()
                for equiv_region_ID in equiv_regions_w_no_primers:
                    ppr_combination_w_no_pp.update(list(filter(lambda ppr: equiv_region_ID in ppr.getFwdRevIDs(), ppr_combo)))
                    #ppr_combinations_w_no_pp.append( set(filter(lambda ppr: equiv_region_ID in ppr.getID(), ppr_combo)) )
                
        return (ppr_combo_primer_sets, ppr_combination_w_no_pp)


    def createPPRComboPrimerSets(self, ppr_per_cc, filtered_equiv_region_primers, primer_num_nontarget_hybs, rxn_temp_C=57, primer_concentration=1.0*10**-8):
        '''Returns a list that contains lists of primer pairs that are all pairwise thermodynamically compatible and that amplify the TIGs of ppr_combo'''

        thermo_calc = primer3.thermoanalysis.ThermoAnalysis(mv_conc=50, dv_conc=1.5, dntp_conc=0.25, temp_c=rxn_temp_C, dna_conc=50)
        rxn_temp_K = 273.15 + rxn_temp_C
        RT = 1.98717 * rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)
        calc_frac_duplexed = lambda dG: (primer_concentration * exp(-dG/RT))/(1 + primer_concentration * exp(-dG/RT))

        ppr_combo_primer_sets = []

        # Get the allowed and/or reduced set of primer pair start positions (fwd primer 5', rev primer 5') to be considered for each PPR
        pp_for_ppr_fwd_rev_5p = defaultdict(dict)
        #pp_starts_for_ppr = {}
        #amplicon_deviation_from_ideal_5p = {}
        candidates_per_ppr = {}
        for ppr in chain.from_iterable(filter(lambda x: len(x)>1, ppr_per_cc)):
            fwd_region_ID, rev_region_ID = ppr.getFwdRevIDs()
            #amplicon_deviation_from_ideal_5p[ppr] = {}
            candidate_pairs = []
            if (ppr.isDisjoint()):
                #fwd_genomic_starts = set(chain.from_iterable(map(itemgetter(0), filtered_equiv_region_primers[fwd_region_ID])))
                #rev_genomic_starts = set(chain.from_iterable(map(itemgetter(0), filtered_equiv_region_primers[rev_region_ID])))
                primer_starts_data = ppr.getPrimerStartsForIsoformGroupsTuple()
                #pp_starts_for_ppr[ppr] = set(filter(lambda fwd_rev_5p_pos: fwd_rev_5p_pos[0] in fwd_genomic_starts and fwd_rev_5p_pos[1] in rev_genomic_starts, primer_starts_data.keys()))
                #for (fwd_5p_pos, rev_5p_pos) in pp_starts_for_ppr[ppr]:
                #    amplicon_deviation_from_ideal_5p[ppr][(fwd_5p_pos, rev_5p_pos)] = primer_starts_data[(fwd_5p_pos, rev_5p_pos)][2]
                sys.exit(1) # Precompute num_lc for primers, then below add to candidate_pairs and then to sorting criteria.
                for (fwd_primer, rev_primer) in product(filtered_equiv_region_primers[fwd_region_ID], filtered_equiv_region_primers[rev_region_ID]):
                    fwd_5p_pos = fwd_primer.genome_5p_pos[0]
                    rev_5p_pos = rev_primer.genome_5p_pos[0]
                    if ((fwd_5p_pos, rev_5p_pos) in primer_starts_data):
                        #amplicon_deviation_from_ideal_5p[ppr][(fwd_5p_pos, rev_5p_pos)] = primer_starts_data[(fwd_5p_pos, rev_5p_pos)][2]
                        target_amplicon_deviation_len = primer_starts_data[(fwd_5p_pos, rev_5p_pos)][2]
                        sum_indiv_penalties = fwd_primer.thermo_penalty + rev_primer.thermo_penalty
                        sum_aux_penalties = fwd_primer.aux_3p_penalty + rev_primer.aux_3p_penalty
                        pair_penalty_term = self.oligo_thermo.calcPrimerPairThermoDetails(fwd_primer.seq, rev_primer.seq)
                        total_primer_pair_penalty = sum_indiv_penalties + pair_penalty_term
                        # TODO: deviation value not the same as _Dont/_Abutt ?
                        total_num_nontarget_hybs = primer_num_nontarget_hybs[fwd_primer] + primer_num_nontarget_hybs[rev_primer]
                        candidate_pairs.append((fwd_primer, rev_primer, total_primer_pair_penalty, round(total_primer_pair_penalty,5),
                                                sum_aux_penalties, fwd_primer.len+rev_primer.len, target_amplicon_deviation_len, total_num_nontarget_hybs))

                candidate_pairs = list(filter(lambda p: p[3] < 0.01, candidate_pairs))
                shuffle(candidate_pairs)
                candidate_pairs.sort(key=itemgetter(3,7,4,6,5))
                candidates_per_ppr[ppr] = candidate_pairs
            else:
                dist_to_ideal_start_fwd, dist_to_ideal_start_rev = ppr.calcOptFwdRevPrimerPositions()
                filtered_equiv_region_primers[fwd_region_ID].sort(key=lambda x: (x[6],x[7],dist_to_ideal_start_fwd[x[0][0]]))
                filtered_equiv_region_primers[rev_region_ID].sort(key=lambda x: (x[6],x[7],dist_to_ideal_start_rev[x[0][0]]))

                sys.exit(1) # Precompute num_lc for primers, then below add to candidate_pairs and then to sorting criteria.
                for fwd_primer, rev_primer in product(filtered_equiv_region_primers[fwd_region_ID], filtered_equiv_region_primers[rev_region_ID]):
                    if (ppr.areAllowedPrimersForIsoformGroupsTuple(fwd_primer[0], fwd_primer[2], rev_primer[0], rev_primer[2])):
                        fwd_5p_pos = fwd_primer.genome_5p_pos[0]
                        rev_5p_pos = rev_primer.genome_5p_pos[0]
                        #amplicon_deviation_from_ideal_5p[ppr][(fwd_5p_pos, rev_5p_pos)] = dist_to_ideal_start_fwd[fwd_5p_pos] + dist_to_ideal_start_rev[rev_5p_pos]
                        mean_target_amplicon_deviation_len = int((dist_to_ideal_start_fwd[fwd_5p_pos] + dist_to_ideal_start_rev[rev_5p_pos])/2.0)
                        sum_indiv_penalties = fwd_primer.thermo_penalty + rev_primer.thermo_penalty
                        sum_aux_penalties = fwd_primer.aux_3p_penalty + rev_primer.aux_3p_penalty
                        pair_penalty_term = self.oligo_thermo.calcPrimerPairThermoDetails(fwd_primer.seq, rev_primer.seq)
                        total_primer_pair_penalty = sum_indiv_penalties + pair_penalty_term
                        total_num_nontarget_hybs = primer_num_nontarget_hybs[fwd_primer] + primer_num_nontarget_hybs[rev_primer]
                        candidate_pairs.append((fwd_primer, rev_primer, total_primer_pair_penalty, round(total_primer_pair_penalty,5),
                                                sum_aux_penalties, fwd_primer.len+rev_primer.len, mean_target_amplicon_deviation_len, total_num_nontarget_hybs)) 
                        
                candidate_pairs = list(filter(lambda p: p[3] < 0.01, candidate_pairs))
                shuffle(candidate_pairs)
                candidate_pairs.sort(key=itemgetter(3,7,4,6,5))
                candidates_per_ppr[ppr] = candidate_pairs[0:250]

        if (any(map(lambda v: len(v)==0, candidates_per_ppr.values()))): # pp_starts_for_ppr
            return []

        if (len(candidates_per_ppr)>0):
            # Create query_primers, except use query sequence as ID. Record all pprID.pp_num_fwd/rev primers that have the same query sequence.
            query_primers = {}
            query_seq_to_ID_map = defaultdict(list)
            for ppr, primer_pairs in candidates_per_ppr.items():
                ppr_ID = ppr.getID()
                for pp_num, (fwd_primer, rev_primer) in enumerate(map(itemgetter(0,1), primer_pairs)):
                    fwd_id = "%s.%d_fwd" % (ppr_ID, pp_num)
                    rev_id = "%s.%d_rev" % (ppr_ID, pp_num)
                    query_primers[fwd_primer.seq] = fwd_primer.seq     
                    query_primers[rev_primer.seq] = rev_primer.seq     # To prevent including duplicate primers, 
                    query_seq_to_ID_map[fwd_primer.seq].append(fwd_id) # since primer pairs will have a singleton primers in common
                    query_seq_to_ID_map[rev_primer.seq].append(rev_id)

            print("\t\t\tINFO: %d query primers" % len(query_primers.keys()), file=sys.stderr, flush=True)
            sys.exit(1) # Replace with blastn
            all_vsearch_results, primers_w_no_perfect_match = runThoroughVsearch(query_primers, self.transcriptome_vsearch_udb, self.tempdir)

            # First-pass removal of primers that are too promiscuous. Second pass filtering done below.
            unique_isoform_IDs = set()
            for igs_tuples_by_ID in map(lambda ppr: ppr.getIsoformGroupsDualPrimedByID(), candidates_per_ppr.keys()):
                unique_isoform_IDs.update(chain.from_iterable(igs_tuples_by_ID))
            promiscuity_cutoff = 2 * designParams.max_primer_promiscuity + len(unique_isoform_IDs)

            primer_match_counts = Counter([tup[0] for tup in all_vsearch_results]).most_common()
            # TODO: TEMPORARY
            print("\t\t\tpromiscuity_cutoff = %d" % promiscuity_cutoff, file=sys.stderr)
            for primer_count in filter(lambda x: x[1] > promiscuity_cutoff, primer_match_counts):
                print("\t\t\t%s\t%d" % primer_count, file=sys.stderr)
            promiscuous_primers = set(map(itemgetter(0), filter(lambda x: x[1] > promiscuity_cutoff, primer_match_counts)))

            # Expand all_vsearch_results by substituting each query sequence ID with each associated ID.
            print("\t\t\tINFO: expanding vsearch results, ", file=sys.stderr, flush=True, end=''),
            expanded_vsearch_results = []
            num_vsearch_lines_used = 0
            for tup in filter(lambda t: t[0] not in promiscuous_primers, all_vsearch_results):
                num_vsearch_lines_used += 1
                for proper_ID in query_seq_to_ID_map[tup[0]]:
                    new_tup = (proper_ID,) + tup[1:]
                    expanded_vsearch_results.append(new_tup)
            print(" %d -> %d" % (num_vsearch_lines_used, len(expanded_vsearch_results)), file=sys.stderr, flush=True)
            all_vsearch_results = expanded_vsearch_results

            print("\t\t\tINFO: grouping vsearch results by primer pair", file=sys.stderr, flush=True)
            sys.exit(1) # Replace with blastn equivalent
            unverified_pp_grouped_by_ppr, hyb_events_per_primer = groupVsearchResultsByPrimerPair(all_vsearch_results, self.oligo_thermo)

            # Set the one primer to use for a primer 5' start position
            verified_primer_pairs_per_ppr = {}
            fwd_primer_occurrences_for_5p_pos = defaultdict(list)
            rev_primer_occurrences_for_5p_pos = defaultdict(list)
            fwd_primer_for_5p_pos = {}
            rev_primer_for_5p_pos = {}
            for ppr, primer_pairs in candidates_per_ppr.items():
                ppr_ID = ppr.getID()
                print("\t\t\t%s" % ppr_ID, file=sys.stderr, flush=True)
                igs_tuples_by_ID = ppr.getIsoformGroupsDualPrimedByID()
                verified_primer_pairs = filterPrimerPairsByGlobalSpecificityAndPromiscuity(ppr_ID, igs_tuples_by_ID, primer_pairs, unverified_pp_grouped_by_ppr, hyb_events_per_primer,
                                                                                           self.complete_TIG_IDs_w_nonlocal_members, self.oligo_thermo, designParams.max_primer_promiscuity,
                                                                                           self.transcriptome_ref)

                # TODO: handle more cleanly. May require breaking this function up in to multiple functions
                if (len(verified_primer_pairs)==0):
                    return ([],[ppr])

                verified_primer_pairs.sort(key=itemgetter(3,7,4,6,5))
                verified_primer_pairs_per_ppr[ppr] = verified_primer_pairs
                for tup in verified_primer_pairs:
                    fwd_primer, rev_primer = tup[0:2]
                    fwd_5p_pos = fwd_primer.genome_5p_pos[0]
                    rev_5p_pos = rev_primer.genome_5p_pos[0]

                    fwd_primer_occurrences_for_5p_pos[fwd_5p_pos].append(fwd_primer)
                    rev_primer_occurrences_for_5p_pos[rev_5p_pos].append(rev_primer)

            for fwd_5p_pos, primers in fwd_primer_occurrences_for_5p_pos.items():
                fwd_primer_for_5p_pos[fwd_5p_pos] = Counter(primers).most_common(1)[0][0]
            for rev_5p_pos, primers in rev_primer_occurrences_for_5p_pos.items():
                rev_primer_for_5p_pos[rev_5p_pos] = Counter(primers).most_common(1)[0][0]

            for ppr, verified_primer_pairs in verified_primer_pairs_per_ppr.items():
                all_5p_pos_pairs = set()
                for tup in verified_primer_pairs:
                    fwd_primer, rev_primer = tup[0:2]
                    fwd_5p_pos = fwd_primer.genome_5p_pos[0]
                    rev_5p_pos = rev_primer.genome_5p_pos[0]
                    all_5p_pos_pairs.add( (fwd_5p_pos,rev_5p_pos) )

                    if (fwd_primer_for_5p_pos[fwd_5p_pos] == fwd_primer and rev_primer_for_5p_pos[rev_5p_pos] == rev_primer):
                        for fwd_5p_pos, rev_5p_pos in product(fwd_primer.genome_5p_pos, rev_primer.genome_5p_pos):
                            pp_for_ppr_fwd_rev_5p[ppr][(fwd_5p_pos,rev_5p_pos)] = tup[0:-1]

                # Report statistics on paired 5' positions for which primers couldn't be set
                num_missing = len(all_5p_pos_pairs - set(pp_for_ppr_fwd_rev_5p[ppr]))
                print("\t\t\t%s: Cannot set primers for %d of %d 5' start position pairs" % (ppr.getID(), num_missing, len(all_5p_pos_pairs)), file=sys.stderr)

        # Based on the allowed 5' primer start positions for each cc PPR combination, enumerate fully-defined primer combinations and associated penalties
        all_cc_pps = []
        all_ppr_descriptor = "+".join(list(map(methodcaller("getID"), chain.from_iterable(ppr_per_cc))))
        #expected_igs_tuples_by_ID = tuple(chain.from_iterable(map(methodcaller("getIsoformGroupsDualPrimedByID"), chain.from_iterable(ppr_per_cc))))
        debug_ppr_IDs = set(['F61R58', 'F62R58', 'F61R55'])
        for pprs in ppr_per_cc:
            pprs_descriptor = "+".join(list(map(methodcaller("getID"), pprs)))
            print("\t\t\tPPRs: %s" % pprs_descriptor, file=sys.stderr)

            cc_pps = []
            if (len(pprs)==1): 
                ppr = list(pprs)[0]
                sops = ppr.getDefaultPrimerPairs()
                cc_pps = list(map(lambda sop: sop.getAllAnnotPrimers(), sops))
            else:
                for_STG = {}
                pprs_w_equiv_region = defaultdict(set)
                equiv_region_instances = list(chain.from_iterable(map(methodcaller("getFwdRevIDs"), pprs)))
                for ppr in pprs:
                    erID_fwd, erID_rev = ppr.getFwdRevIDs()

                    if (equiv_region_instances.count(erID_fwd)>1): # If equiv region connects 2 or more PPRs...
                        pprs_w_equiv_region[erID_fwd].add(ppr)
                        if (erID_fwd not in for_STG):
                            for_STG[erID_fwd] = defaultdict(list)

                        for (fwd_5p_pos, rev_5p_pos) in pp_for_ppr_fwd_rev_5p[ppr].keys():
                            for_STG[erID_fwd][fwd_5p_pos].append( ((fwd_5p_pos, rev_5p_pos), ppr) )

                    if (equiv_region_instances.count(erID_rev)>1): # If equiv region connects 2 or more PPRs...
                        pprs_w_equiv_region[erID_rev].add(ppr)
                        if (erID_rev not in for_STG):
                            for_STG[erID_rev] = defaultdict(list)
                        
                        for (fwd_5p_pos, rev_5p_pos) in pp_for_ppr_fwd_rev_5p[ppr].keys():
                            for_STG[erID_rev][rev_5p_pos].append( ((fwd_5p_pos, rev_5p_pos), ppr) )

                # Build graph from which spanning tree solutions can be derived
                STG = nx.Graph()
                for erID in for_STG.keys():
                    assert (len(pprs_w_equiv_region[erID]) > 1)
                    for pp_pos_and_ppr in for_STG[erID].values():
                        if (set(map(itemgetter(1), pp_pos_and_ppr)) == pprs_w_equiv_region[erID]):
                            STG.add_nodes_from(pp_pos_and_ppr)
                            STG.add_edges_from(list(filter(lambda tupAB: tupAB[0][1] != tupAB[1][1], combinations(pp_pos_and_ppr,2)))) # Add edge if ppr1 != ppr2

                # Select the PPR with the fewest primer pairs to use for root nodes (for generating all spanning trees)
                root_ppr = Counter(list(map(itemgetter(1), STG.nodes()))).most_common()[-1][0]

                # shuffle the root nodes then sort to put better ones first. Intended effect is to separate primer start positions
                sorted_root_nodes = []
                K = pp_for_ppr_fwd_rev_5p[root_ppr]
                for (fwd_5p_pos,rev_5p_pos), ppr in map(itemgetter(0,1), filter(lambda x: x[1]==root_ppr, STG.nodes())):
                    sorted_root_nodes.append( (fwd_5p_pos,rev_5p_pos) + K[(fwd_5p_pos,rev_5p_pos)][2:5] )
                shuffle(sorted_root_nodes)
                sorted_root_nodes.sort(key=itemgetter(2,3)) # ,4

                all_valid_spanning_trees = []
                backup_spanning_trees = []
                counter = 0

                # DEBUGGING
                ppr_IDs = set(map(lambda x:x.ID, pprs))
                #if (ppr_IDs == debug_ppr_IDs):
                #    pdb.set_trace()

                for (fwd_5p_pos, rev_5p_pos) in map(itemgetter(0,1), sorted_root_nodes):
                    root_node = ((fwd_5p_pos, rev_5p_pos), root_ppr)
                    #root_node_spanning_trees = find_all_spanning_trees(STG, root_node, len(pprs))
                    #valid_spanning_trees = findAllSpanningTrees(STG, root_node, len(pprs))
                    valid_spanning_trees = findAllSpanningTrees2(STG, [root_node], set([root_node[1]])) # set(), 

                    valid_spanning_trees_w_penalties = []
                    for vst in valid_spanning_trees:
                        vst = (root_node,) + vst
                        # n = (primers_descr,ppr)
                        try:
                            qc_subG = STG.subgraph(vst)
                            assert (nx.is_connected(qc_subG) and len(qc_subG)==len(pprs)), "spanning tree not connected and/or correct size"
                        except AssertionError as ae:
                            print(ae, file=sys.stderr)
                            print("%s,%d\t%d" % (nx.is_connected(qc_subG),len(qc_subG),len(pprs)), file=sys.stderr)
                            print(ppr_IDs, file=sys.stderr)
                            sys.exit(1)
                            
                        annot_pps_w_penalties = list(map(lambda n: n[1].getFwdRevIDs() + pp_for_ppr_fwd_rev_5p[n[1]][n[0]], vst)) # vst.nodes()

                        # Calc sum of the primer interaction penalties
                        already_computed = set(list(map(itemgetter(2,3), annot_pps_w_penalties)) + list(map(itemgetter(3,2), annot_pps_w_penalties)))
                        all_primers_annot = set(chain.from_iterable(map(lambda x:[(x[0],x[2]),(x[1],x[3])], annot_pps_w_penalties)))
                        interact_thermo_penalties = []
                        for annot_pp in filter(lambda pApB: (pApB[0][1],pApB[1][1]) not in already_computed, combinations(all_primers_annot, 2)):
                            p = self.oligo_thermo.calcPrimerPairThermoDetails(annot_pp[0][1].seq, annot_pp[1][1].seq)
                            interact_thermo_penalties.append(p)
                                
                        if (max(interact_thermo_penalties) < 0.01):
                            sys.exit(1) # Add num_lc logic here, and to sort criteria below
                            # Calc sum of the primer interaction penalties
                            total_thermo_penalties = sum(map(itemgetter(4), annot_pps_w_penalties)) + sum(interact_thermo_penalties)
                            total_aux_penalties = sum(map(itemgetter(5), annot_pps_w_penalties))
                            total_deviation = sum(map(itemgetter(6), annot_pps_w_penalties))
                            total_promiscuity = sum(map(itemgetter(7), annot_pps_w_penalties))

                            valid_spanning_trees_w_penalties.append( (all_primers_annot, round(total_thermo_penalties,4), total_aux_penalties, total_deviation, total_promiscuity) )

                    # Select the best solution for this root node
                    if (len(valid_spanning_trees_w_penalties) > 0):
                        valid_spanning_trees_w_penalties.sort(key=itemgetter(1,2,4,3))
                        any_adds_to_diversity = False
                        for candidate_spanning_tree in valid_spanning_trees_w_penalties:
                            any_adds_to_diversity = addsToDiversity(list(map(itemgetter(0), all_valid_spanning_trees)), candidate_spanning_tree[0])
                            if (any_adds_to_diversity):
                                all_valid_spanning_trees.append(candidate_spanning_tree)    # valid_spanning_trees_w_penalties[0])
                                break

                        if (not any_adds_to_diversity):
                            backup_spanning_trees.append(valid_spanning_trees_w_penalties[0])

                    STG.remove_node(root_node)
                        
                    #print("%d %d %d" % (counter, len(all_valid_spanning_trees), len(backup_spanning_trees)), file=sys.stderr, flush=True)
                    counter += 1
                
                    if (len(all_valid_spanning_trees) >= 15):
                        all_valid_spanning_trees.sort(key=itemgetter(1,2,4,3))
                        break
                    elif (len(backup_spanning_trees) >= 50):
                        all_valid_spanning_trees.sort(key=itemgetter(1,2,4,3))
                        backup_spanning_trees.sort(key=itemgetter(1,2,4,3))
                        all_valid_spanning_trees += backup_spanning_trees

                if (len(all_valid_spanning_trees) > 0):
                    diversity_of_spanning_trees = selectDiversityOfPrimers(list(map(itemgetter(0), all_valid_spanning_trees)), 15, 2)
                    cc_pps.extend( diversity_of_spanning_trees )
                else:
                    pdb.set_trace()
                    
            all_cc_pps.append(cc_pps)
            
        #multi_index_generator = self.generateMultiIndexesInOrderSmallestSumRanks( map(lambda x: len(x), all_cc_pps) )
        #multi_index_generator = self.generateMultiIndexesRandomly( map(lambda x: len(x)-1, all_cc_pps) )

        compat_primers_cache = set()
        num_rotations = np.zeros((len(all_cc_pps),), 'f')
        duplicate_soln_found = False
        max_solns = max(map(len, all_cc_pps))
        rotation_ratio = np.zeros((len(all_cc_pps),), 'f')
        incompat_index_pair_occurs = []
        # TODO: need better criterion than '5'
        while (len(ppr_per_cc) == len(all_cc_pps) and len(ppr_combo_primer_sets) < 5 and not duplicate_soln_found and
               not np.all(rotation_ratio > 2.0) and np.max(num_rotations) <= 2*max_solns): 
            #try:
            #    multi_index = multi_index_generator.next()
            #except StopIteration as si:
            #    break
            #print >> sys.stderr, multi_index
            
            # The primers in each individual cc are compatible. Determine here the compatibility of the primers in the different cc's
            incompat_primers_found = False
            #for ((cc1_index,ps1_index), (cc2_index,ps2_index)) in combinations(multi_index,2):
            for (cc1_index, cc2_index) in combinations(range(len(all_cc_pps)),2):
                primers_set1 = all_cc_pps[cc1_index][0] #[ps1_index]
                primers_set2 = all_cc_pps[cc2_index][0] #[ps2_index]
                
                # If, on a target isoform, a fwd primer is upstream of a rev primer or a rev primer is downstream of a fwd
                # primer less than the distance of a max confusion amplicon length, then the primer sets are positionally incompatible.
                for (erID1,primer1), (erID2,primer2) in product(primers_set1, primers_set2):                    
                    if (erID1[0] != erID2[0]):
                        fwd_primer, rev_primer = (primer1,primer2) if (erID1[0]=='F') else (primer2,primer1)
                        for fwd_genome_5p_pos, rev_genome_5p_pos in product(fwd_primer.genome_5p_pos, rev_primer.genome_5p_pos):
                            if ((self.context_strand == '+' and fwd_genome_5p_pos < rev_genome_5p_pos) or (self.context_strand == '-' and fwd_genome_5p_pos > rev_genome_5p_pos)):
                                for target_isoform in self.target_isoforms:
                                    seqlen = target_isoform.getSequenceLength(fwd_genome_5p_pos, rev_genome_5p_pos)
                                    if (seqlen != None and seqlen < designParams.confusion_amplicon_max_len):
                                        incompat_index_pair_occurs.append( (cc1_index, cc2_index) )
                                        incompat_primers_found = True
                                        break
                            if (incompat_primers_found):
                                break
                    if (incompat_primers_found):
                        break

                if (not incompat_primers_found): # No primer interactions were found to created unwanted amplicons on target isoforms. Now check for thermodynamic compatibility.
                    for (erID1,primer1), (erID2,primer2) in product(primers_set1, primers_set2):
                        if ((primer1,primer2) in compat_primers_cache):
                            continue
                        else:
                            interact_thermo_penalty = self.oligo_thermo.calcPrimerPairThermoDetails(primer1.seq,primer2.seq)
                            if (interact_thermo_penalty < 0.01): # or (regions_different_type and interact_thermo_penalty < 0.05)):
                                compat_primers_cache.add( (primer1,primer2) )
                            else:
                                incompat_index_pair_occurs.append( (cc1_index, cc2_index) )
                                incompat_primers_found = True
                                break

                if (incompat_primers_found):
                    rotation_ratio_1 = rotation_ratio[cc1_index] # num_rotations[cc1_index]/(len(all_cc_pps[cc1_index]) + len(ppr_combo_primer_sets))
                    rotation_ratio_2 = rotation_ratio[cc2_index] # num_rotations[cc2_index]/(len(all_cc_pps[cc2_index]) + len(ppr_combo_primer_sets))

                    if (rotation_ratio_1 < rotation_ratio_2 or (rotation_ratio_1 == rotation_ratio_2 and len(all_cc_pps[cc1_index]) >= len(all_cc_pps[cc2_index]))):
                        if (len(all_cc_pps[cc1_index]) > 1):
                            ps_to_move = all_cc_pps[cc1_index].pop(0) # ps1_index) # Move to the end of the list
                            all_cc_pps[cc1_index].append(ps_to_move)
                        num_rotations[cc1_index] += 1
                        rotation_ratio[cc1_index] = num_rotations[cc1_index]/len(all_cc_pps[cc1_index])
                    else:
                        if (len(all_cc_pps[cc2_index]) > 1):
                            ps_to_move = all_cc_pps[cc2_index].pop(0) # ps2_index)
                            all_cc_pps[cc2_index].append(ps_to_move)
                        num_rotations[cc2_index] += 1
                        rotation_ratio[cc2_index] = num_rotations[cc2_index]/len(all_cc_pps[cc2_index])

                        #multi_index_generator = self.generateMultiIndexesInOrderSmallestSumRanks( map(lambda x: len(x), all_cc_pps) )
                    break

            if (not incompat_primers_found):
                chrom, OS_ID = self.olap_set.getChromosomeAndID()
                # mi = (cc_index,primer_set_index)
                all_cc_primers = list(chain.from_iterable( map(lambda index: all_cc_pps[index][0], range(len(all_cc_pps))) )) # mi[1], multiindex
                new_sop = SetOfPrimers(chrom, OS_ID, all_ppr_descriptor, len(ppr_combo_primer_sets), all_cc_primers) # expected_igs_tuples_by_ID, 
                duplicate_soln_found = any(map(lambda sop: new_sop.isSameAs(sop), ppr_combo_primer_sets))

                if (not duplicate_soln_found):
                    print("\t\t\tINFO: Found solution #%d" % (len(ppr_combo_primer_sets)+1,), file=sys.stderr, flush=True)
                    ppr_combo_primer_sets.append( new_sop )

                    # To enable primer diversity in solutions, remove any new multi_index that includes an index used in generating this solution...
                    if (len(ppr_combo_primer_sets) < 5):
                        for cc_index in range(len(all_cc_pps)): #, ps_index in multi_index: 
                            if (len(all_cc_pps[cc_index]) > 1):
                                ps_to_move = all_cc_pps[cc_index].pop(0) # ps_index # Move to the end of the list
                                all_cc_pps[cc_index].append(ps_to_move)
                            num_rotations[cc_index] += 1
                            rotation_ratio[cc_index] = num_rotations[cc_index]/len(all_cc_pps[cc_index])

                        #multi_index_generator = self.generateMultiIndexesInOrderSmallestSumRanks( map(lambda x: len(x), all_cc_pps) )
                        #multi_index_generator = self.generateMultiIndexesRandomly( map(lambda x: len(x)-1, all_cc_pps) )

        if (len(ppr_combo_primer_sets)==0):
            C = Counter(incompat_index_pair_occurs)
            most_common = C.most_common(1)
            index1, index2 = most_common[0][0]
            most_freq_incompat_pprs = ppr_per_cc[index1] | ppr_per_cc[index2]
            print("\t\t\tINFO: No solutions found", file=sys.stderr, flush=True)
        else:
            most_freq_incompat_pprs = None
            
        return ppr_combo_primer_sets, most_freq_incompat_pprs


    def generatePrimersInOrderSmallestSumRanks(self, equiv_region_primers, ordered_equiv_region_IDs, G):
        '''equiv_region_primers is a dictionary of list of primers for each equiv region, with primers listed in rank order of goodness.
        This method returns list of tuples (equiv region ID, primer) that are composed of one primer
        from each list and ordered such that the tuples with lower sum of primer ranks are returned first.'''

        L = list(map(lambda erID: equiv_region_primers[erID], ordered_equiv_region_IDs))

        DG = nx.DiGraph()
        for i,j in zip(range(len(L)-1), range(1,len(L))):
            for (i_ind,primer_i), (j_ind,primer_j) in product(enumerate(L[i]),enumerate(L[j])):
                if (G.has_edge(primer_i,primer_j)):
                    DG.add_edge( primer_i, primer_j, weight=j_ind )
            
        for j_ind, primer_j in enumerate(L[0]):
            DG.add_edge('B', primer_j, weight=j_ind)

        for primer_i in L[-1]:
            DG.add_edge(primer_i, 'E', weight=0)

        if (nx.has_path(DG, 'B', 'E')):
            compat_primers = set()
            for node in G.nodes():
                compat_primers.add( (node,node) )
            for (n1,n2) in G.edges():
                compat_primers.add( (n1,n2) )
                compat_primers.add( (n2,n1) )

            for path in nx.shortest_simple_paths(DG, 'B', 'E', weight='weight'):
                path_primers = path[1:-1]
                reqd_edges = set(product(path_primers,path_primers))
                if (reqd_edges <= compat_primers):
                    assert (len(ordered_equiv_region_IDs) == len(path_primers))
                    yield zip(ordered_equiv_region_IDs, path_primers)
                else:
                    num_missing = len(reqd_edges - compat_primers)
                    #print >> sys.stderr, "Need %d" % num_missing
                    

    def generateMultiIndexesRandomly(self, ranges):
        i_values = range(len(ranges))
        while (True):
            yield list(map(lambda i: (i,randint(0,ranges[i])), i_values))


    def generateMultiIndexesInOrderSmallestSumRanks(self, ranges):
        if (0 in ranges):
            raise StopIteration
        
        L = list(map(lambda n: list(range(n)), ranges))

        DG = nx.DiGraph()
        for i,j in zip(range(len(L)-1), range(1,len(L))):
            for i_ind, j_ind in product(L[i],L[j]):
                DG.add_edge( (i,i_ind), (j,j_ind), weight=j_ind )
            
        for j_ind in L[0]:
            DG.add_edge('B', (0,j_ind), weight=j_ind)

        for i_ind in L[-1]:
            DG.add_edge( (len(L)-1,i_ind), 'E', weight=0)

        if (nx.has_path(DG, 'B', 'E')):
            for path in nx.shortest_simple_paths(DG, 'B', 'E', weight='weight'):
                yield path[1:-1]


    def scanIndivPrimersAgainstLocalTranscriptome_DEPRACATED(self, expected_isoforms_per_primer_seq, ppr_combo, all_equiv_region_primers):
        '''Evaluate each primer from each set of equiv region candidates individually, without regard to a priming partner. Filter out any primers that have 
        unexpected/unwanted hybridization potential. Only return results for which frac_duplexed of primer to target > 0.01. primer and target sequences are 
        returned in 5'->3' orientation.'''
        bad_primer_seqs = set()

        local_transcriptome_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=self.tempdir, delete=True)
        isoforms_seqs = self.olap_set.getAllOverlappingIsoformsSeqs()
        for isoform_ID, isoform_seq in isoforms_seqs.items():
            local_transcriptome_fasta.write(">%s\n%s\n" % (isoform_ID, isoform_seq))
        local_transcriptome_fasta.flush()

        ppr_unique_target_isoform_IDs = set()
        for ppr in ppr_combo:
            igs_tuples_isoforms_ID = ppr.getIsoformGroupsDualPrimedIsoformsID()
            ppr_unique_target_isoform_IDs.update(igs_tuples_isoforms_ID)
        ppr_unique_target_isoform_IDs &= self.target_isoforms_IDs

        # First, filter out primers that align to antitarget mRNAs
        primers = {}
        for equiv_region_ID, equiv_region_primers in all_equiv_region_primers.items():
            for primer_num, singleton_primer in enumerate(equiv_region_primers):
                primer_id = "%s.%d" % (equiv_region_ID, primer_num)
                primers[primer_id] = singleton_primer.seq
        ids_of_antitarget_matching_primers = scanPrimersAgainstAntitargets(primers, self.antitargets_fasta_name, self.tempdir, 'B')

        # Second, align to local transcriptome those primers not matching to antitargets
        primers_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=self.tempdir, delete=True)
        for primer_id, primer_seq in primers.items():
            if (primer_id not in ids_of_antitarget_matching_primers):
                primers_fasta.write(">%s\n%s\n" % (primer_id,primer_seq))
            else:
                bad_primer_seqs.add(primer_seq)
        primers_fasta.flush()

        ip_vsearch = tempfile.NamedTemporaryFile(mode='w+t', suffix=".vsearch", dir=self.tempdir, delete=True)
        vsearch_cmd = ['vsearch', '--usearch_global', primers_fasta.name, '-db', local_transcriptome_fasta.name,
                       '--id', '0.75', '--minwordmatches', '2', '--wordlength', '5', '--dbmask', 'none', '--qmask', 'none', '--strand', 'both',
                       '--userout', ip_vsearch.name, '--userfields', 'query+qstrand+qlo+qhi+qrow+target+tstrand+tilo+tihi+trow+mism+opens',
                       '--maxaccepts', '1000000', '--maxgaps', '1', '--maxdiffs', '4', '--maxrejects', '4096', '--threads', '3', '--maxseqlength', '1000000', '--query_cov', '0.93']

        try:
            vsearch_output = check_output(vsearch_cmd, stderr=subprocess.STDOUT)
        except CalledProcessError as cpe:
            pdb.set_trace()

        ip_vsearch.flush()
        ip_vsearch.seek(0)

        # Compile local transcriptome matches for each primer
        local_primer_matches = defaultdict(set)
        primed_isoforms_counts = defaultdict(Counter)
        for line in ip_vsearch:
            query, qstrand, qlo, qhi, qrow, target, tstrand, tilo, tihi, trow, mism, opens = line.strip().split('\t')
            equiv_region_ID, primer_num = query.rsplit('.',1)
            primer_fwdrev = "fwd" if (equiv_region_ID[0] == 'F') else "rev"
            match_fwdrev = "fwd" if (qstrand=='+') else "rev"
            query_len = int(qlo) if (qstrand=='-') else int(qhi)

            if (opens == '0'):
                query_seq, template_seq = qrow, trow
            else:
                query_seq = qrow.replace('-','')
                template_seq = trow.replace('-','')

            if (len(query_seq) == query_len): # Skip partial matches at the 5'/3' ends of mRNAs
                frac_duplexed = 0
                # indels in the 3' end or more than 1bp mismatch -> effectively no primer extension
                # Point here is to catch all even low-likelihood priming hybridizations
                if (qstrand == '+'): # and sum(map(lambda x:x[0]!=x[1], zip(qrow[-4:],trow[-4:]))) <= 1): #  and '-' not in qrow[-6:] and '-' not in trow[-6:] and 
                    template_seq = template_seq[::-1].translate(DNA_complement_table)
                    any_hyb = self.oligo_thermo.calcHeterodimer(query_seq, template_seq)
                    end_hyb = self.oligo_thermo.calcEndStability(query_seq, template_seq)
                    min_dG = min(any_hyb.dg, end_hyb.dg)
                    frac_duplexed = self.oligo_thermo.calcFracDuplexed(min_dG)
                    
                else: #if (qstrand == '-' and sum(map(lambda x:x[0]!=x[1], zip(qrow[0:4],trow[0:4]))) <= 1): # and '-' not in qrow[0:6] and '-' not in trow[0:6] 
                    query_seq = query_seq[::-1].translate(DNA_complement_table)
                    any_hyb = self.oligo_thermo.calcHeterodimer(query_seq, template_seq)
                    end_hyb = self.oligo_thermo.calcEndStability(query_seq, template_seq)
                    min_dG = min(any_hyb.dg, end_hyb.dg)
                    frac_duplexed = self.oligo_thermo.calcFracDuplexed(min_dG)

                if (frac_duplexed >= 0.01):
                    data_tup = (primer_fwdrev, match_fwdrev, query_seq, target, int(tilo), int(tihi), template_seq, opens != "0", frac_duplexed)
                    local_primer_matches[query_seq].add(data_tup)
                    if (target in ppr_unique_target_isoform_IDs):
                        primed_isoforms_counts[query_seq][target] += 1

        ip_vsearch.close()
        primers_fasta.close()
        local_transcriptome_fasta.close()

        # Identify "bad" primers that prime more than once on *any* dual-primed PPR target isoform OR that prime on an unintended *dual-primed PPR target isoform*
        # This logic should catch primers that also act oppositely (eg forward primer acting as a reverse primer)
        for primer_seq, primer_target_counts in primed_isoforms_counts.items():
            actual_targets_primed = set(primer_target_counts.keys())
            if (primer_target_counts.most_common(1)[0][1] > 1 or len(actual_targets_primed - expected_isoforms_per_primer_seq[primer_seq])>0):
                bad_primer_seqs.add(primer_seq)
                    
        # Filter out "bad" singleton primers
        filtered_equiv_region_primers = {}
        equiv_regions_w_no_primers = set()
        for equiv_region_ID in all_equiv_region_primers.keys():
            filtered_equiv_region_primers[equiv_region_ID] = list(filter(lambda sp: sp.seq not in bad_primer_seqs, all_equiv_region_primers[equiv_region_ID]))
            if (len(filtered_equiv_region_primers[equiv_region_ID]) == 0):
                equiv_regions_w_no_primers.add(equiv_region_ID)

        if (any(map(lambda v: len(v)==0, filtered_equiv_region_primers.values()))):
            filtered_equiv_region_primers = None
            
        return (filtered_equiv_region_primers, local_primer_matches, equiv_regions_w_no_primers)


    def buildPrimerConnectionsGraph2(self, equiv_region_primers, ppr_combo):
        ordered_equiv_region_IDs = equiv_region_primers.keys()
        ordered_equiv_region_IDs.sort() # puts the Fwd regions (F#) before the Rev regions (R#)

        eprpDG = OrderedDiGraph()
        eprpDG.graph['ppr_combo'] = ppr_combo
        eprpDG.graph['disallowed_primer_combos'] = set()
        eprpDG.graph['pair_penalty'] = {}
        eprpDG.graph['ordered_equiv_region_IDs'] = ordered_equiv_region_IDs

        for equiv_region_ID, primers in equiv_region_primers.items():
            eprpDG.add_nodes_from(primers, equiv_region_ID=equiv_region_ID)
        eprpDG.add_nodes_from(['B', 'E'])

        # Critical: Add edges in rank priority order. Reason for use of the OrderedDiGraph data structure.
        for i in range(len(ordered_equiv_region_IDs)-1):
            id1 = ordered_equiv_region_IDs[i]
            id2 = ordered_equiv_region_IDs[i+1]
            eprpDG.add_edges_from(product(equiv_region_primers[id1], equiv_region_primers[id2]))

        eprpDG.add_edges_from(product('B', equiv_region_primers[ordered_equiv_region_IDs[0]]))
        eprpDG.add_edges_from(product(equiv_region_primers[ordered_equiv_region_IDs[-1]], 'E'))

        assert (nx.has_path(eprpDG,'B','E')), "eprpDG does not have a solution"

        return eprpDG


    def evalPrimersAgainstGlobalTranscriptome_DEPRACATED(self, ppr_combo, all_equiv_region_primers, transcriptome_ref, oligo_thermo, tempdir):
        '''Return a list of verified primer pairs for each PPR '''

        # Perform vsearch of primers against complete transcriptome
        primers_fasta = tempfile.NamedTemporaryFile(suffix=".fa", dir=tempdir, delete=True)
        for equiv_region_ID, equiv_region_primers in all_equiv_region_primers.items():
            for primer_num, singleton_primer in enumerate(equiv_region_primers):
                primers_fasta.write(">%s.%d\n%s\n" % (equiv_region_ID, primer_num, singleton_primer.seq))
        primers_fasta.flush()

            
        # Gather match data for each primer, for each chunk of the transcriptome indexed for usearch
        grouped_primer_matches = {}
        for udb_name in glob.glob("/raid1/projects/CGDB/models/indexes/usearch/CGDB_*.udb"):
            ip_usearch = tempfile.NamedTemporaryFile(suffix=".usearch", dir=tempdir, delete=True)
            usearch_cmd = ['usearch', '-usearch_global', primers_fasta.name, '-threads', '2', '-db', udb_name, '-id', '0.8', '-maxaccepts', '0', '-maxrejects', '256',
                           '-userout', ip_usearch.name, '-userfields', 'query+qstrand+qrow+target+tstrand+trow+id+mism+opens', '-strand', 'both', '-maxgaps', '1', '-maxdiffs', '3']
            try:
                usearch_output = check_output(usearch_cmd, stderr=subprocess.STDOUT)
            except CalledProcessError as cpe:
                pdb.set_trace()

            ip_usearch.flush()
            ip_usearch.seek(0)

            for line in ip_usearch:
                query, qstrand, qrow, target, tstrand, trow, frac_id, mism, opens = line.strip().split('\t')
                equiv_region_ID, primer_num = query.rsplit('.',1)
                fwdrev = "fwd" if (equiv_region_ID[0] == 'F') else "rev"
                is_fwd = fwdrev == "fwd"

                # usearch global doesn't give template sequence alignment start/stop
                query_seq = qrow.replace('-','') if ('-' in qrow) else qrow
                template_seq = trow.replace('-','') if ('-' in trow) else trow

                frac_duplexed = 0
                if (is_fwd and qstrand == '+'):
                    # indels in the 3' end or more than 1bp mismatch -> effectively no primer extension
                    if ('-' not in qrow[-6:] and '-' not in trow[-6:] and sum(map(lambda x:x[0]!=x[1], zip(qrow[-6:],trow[-6:]))) <= 1):
                        template_seq_rc = template_seq[::-1].translate(DNA_complement_table)
                        any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq_rc)
                        end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq_rc)
                        min_dG = min(any_hyb.dg, end_hyb.dg)
                        frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)

                elif (not is_fwd and qstrand == '-'):
                    # indels in the 3' end or more than 1bp mismatch -> effectively no primer extension
                    if ('-' not in qrow[0:6] and '-' not in trow[0:6] and sum(map(lambda x:x[0]!=x[1], zip(qrow[0:6],trow[0:6]))) <= 1):
                        query_seq = query_seq[::-1].translate(DNA_complement_table)
                        any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                        end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                        min_dG = min(any_hyb.dg, end_hyb.dg)
                        frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)

                if (frac_duplexed >= 0.01):
                    tilo = str(transcriptome_ref[target]).index(template_seq)
                    tihi = tilo + len(template_seq)
                    data_tup = (fwdrev, qstrand, qrow, target, tilo, tihi, trow, opens!="0")
                    if (not grouped_primer_matches.has_key(equiv_region_ID)):
                        grouped_primer_matches[equiv_region_ID] = defaultdict(list)
                    grouped_primer_matches[equiv_region_ID][int(primer_num)].append(data_tup)

            ip_usearch.close()
        primers_fasta.close()

            
        # Evaluate candidate primer pairs for each PPR
        ppr_pp = {}
        for ppr in ppr_combo:
            ppr_pp[ppr] = []
            expected_igs_tuple_by_ID = ppr.getIsoformGroupsDualPrimedByID()
            fwd_region_ID, rev_region_ID = ppr.getFwdRevIDs()
            usearch_all_fwd_primer = grouped_primer_matches[fwd_region_ID]
            usearch_all_rev_primer = grouped_primer_matches[rev_region_ID]
            for fi, usearch_fwd_primer_results in usearch_all_fwd_primer.items():
                fwd_singleton_primer = all_equiv_region_primers[fwd_region_ID][fi]
                for ri, usearch_rev_primer_results in usearch_all_rev_primer.items():
                    rev_singleton_primer = all_equiv_region_primers[rev_region_ID][ri]
                    usearch_results = usearch_fwd_primer_results + usearch_rev_primer_results
                    verify_success, categorized_isoform_groups = \
                        ppr.verifyIsoformGroupsAreDistinctBySequencing2((fwd_singleton_primer, rev_singleton_primer), usearch_results, oligo_thermo,
                                                                        expected_igs_tuple_by_ID, self.target_isoforms_IDs,
                                                                        transcriptome_ref, self.complete_TIG_IDs_w_nonlocal_members)
                    if (verify_success):
                        pair_penalty_term = oligo_thermo.calcPrimerPairThermoDetails(fwd_singleton_primer.seq, rev_singleton_primer.seq)
                        pair_penalty = fwd_singleton_primer.thermo_penalty + rev_singleton_primer.thermo_penalty + pair_penalty_term
                        ppr_pp[ppr].append( ((fwd_singleton_primer, rev_singleton_primer, pair_penalty), categorized_isoform_groups) )

        return ppr_pp


    def constrainPrimersPerEquivRegion(self, eprG, cc, equiv_regions_lookup):
        all_edges_possible = True
        equiv_region_primers = None
        this_cc_ppr_combo = set()

        # First, compile all of the allowed forward and reverse primer specifications 
        # for each equiv primer region, for each individual ppr in which it occurs
        subG = eprG.subgraph(cc)

        allowed_primer_specs_per_equiv_region = defaultdict(list)

        for region_ID_A, region_ID_B, edge_data in subG.edges(data=True):
            fwd_region_ID, rev_region_ID = (region_ID_A, region_ID_B) if (region_ID_A[0]=='F' and region_ID_B[0]=='R') else (region_ID_B, region_ID_A)
            assert ((fwd_region_ID[0]=='F' and rev_region_ID[0]=='R')), "Should get a Fwd and Rev region. Got %s and %s" % (fwd_region_ID, rev_region_ID)
            ppr = edge_data['ppr']
            assert (ppr != None)

            this_cc_ppr_combo.add(ppr)

            fwd_primer_specs_this_ppr = set()
            for fwd_genomic_5ps, fwd_primer_lens in ppr.getDecoupledPrimerSpecsForIsoformGroupsTuple("Fwd"):
                fwd_primer_specs_this_ppr.update( list(product(fwd_genomic_5ps, fwd_primer_lens)) )
            allowed_primer_specs_per_equiv_region[fwd_region_ID].append(fwd_primer_specs_this_ppr)

            rev_primer_specs_this_ppr = set()
            for rev_genomic_5ps, rev_primer_lens in ppr.getDecoupledPrimerSpecsForIsoformGroupsTuple("Rev"): 
               rev_primer_specs_this_ppr.update( list(product(rev_genomic_5ps, rev_primer_lens)) )
            allowed_primer_specs_per_equiv_region[rev_region_ID].append(rev_primer_specs_this_ppr)

        # Intersect the allowed primer specs for each equiv primer region.
        common_allowed_primer_specs_per_equiv_region = {}
        for region_ID, list_of_spec_sets in allowed_primer_specs_per_equiv_region.items():
            common_allowed_primer_specs_per_equiv_region[region_ID] = set.intersection(*list_of_spec_sets) if (len(list_of_spec_sets) > 1) else list_of_spec_sets[0]
            if (len(common_allowed_primer_specs_per_equiv_region[region_ID]) == 0):
                all_edges_possible = False
                break

        if (all_edges_possible):
            # Get primer sequences for the allowed primer start position+lengths in each fwd/rev region
            equiv_region_primers = defaultdict(list)
            for equiv_region_ID in cc: 
                equiv_region = equiv_regions_lookup[equiv_region_ID]
                fwdrev = equiv_region.getType()
                genomic_starts_and_lens = defaultdict(list)
                for genomic_5p_start, primer_len in common_allowed_primer_specs_per_equiv_region[equiv_region_ID]:
                    primer_seq = equiv_region.getSubseq(genomic_5p_start, primer_len)
                    too_many_repeats = hasTooManyRepeats(primer_seq)
                    
                    Tm, frac_duplexed, has_hairpin, penalty = self.oligo_thermo.calcPrimerThermoDetails(fwdrev, primer_seq)
                    # Not currently used:
                    # (any(map(lambda s: s in primer_seq, designParams.disallowed_primer_subseqs)) or
                    # any(map(lambda s: primer_seq.endswith(s), designParams.disallowed_primer_3p_subseqs))))
                    primer_fail_QC = (too_many_repeats or has_hairpin or penalty > 0.01 or Tm < designParams.min_primer_Tm or Tm > designParams.max_primer_Tm)
                    if (not primer_fail_QC):
                        alt_genomic_5p_pos, alt_genomic_positions = equiv_region.getAltGenomicPositions(genomic_5p_start, primer_len)
                        nuc_imbalance = abs(3 - primer_seq[-6:].count('C') - primer_seq[-6:].count('G'))
                        AT_3p = 1 if (primer_seq[-1] == 'A' or primer_seq[-1] == 'T') else 0
                        aux_3p_penalty = nuc_imbalance + AT_3p
                        singleton_primer = PrimerSingleton(alt_genomic_5p_pos, primer_seq, primer_len, Tm, frac_duplexed, alt_genomic_positions, round(penalty,4), aux_3p_penalty)
                        equiv_region_primers[equiv_region_ID].append( singleton_primer )

                if (len(equiv_region_primers[equiv_region_ID]) == 0):
                    all_edges_possible = False
                    equiv_region_primers = None
                    break
                
        return equiv_region_primers, this_cc_ppr_combo, all_edges_possible


    def generateRegionsPrimerCombinationWithGraph_ORIG(self, eprpDG, node='B', li=0):
        if (node != 'E'):
            next_region_ID = eprpDG.graph['ordered_equiv_region_IDs'][li]
            for neighbor in filter(lambda n: eprpDG.node[n]['region_ID']==next_region_ID, eprpDG.neighbors(node)):
                for path in self.generateRegionsPrimerCombinationWithGraph(eprpDG, neighbor, li+1):
                    if (all(map(lambda p: eprpDG.has_edge(node,p) or node=='B' or p=='E', path))):
                        new_path = [node]+path
                        if (new_path[0] == 'B' and new_path[-1] == 'E'):
                            new_path = new_path[1:-1]
                            try:
                                assert (node=='B'), "First node is not B"
                                assert (len(new_path) == len(eprpDG.graph['ordered_equiv_region_IDs'])-1), "Path does not include node for each equiv region"
                            except AssertionError as ae:
                                #print >> sys.stderr, ae.message
                                pdb.set_trace()

                            yield dict(zip(eprpDG.graph['ordered_equiv_region_IDs'], new_path))
                        else:
                            try:
                                assert (node!='B'), "First node is B"
                                assert (new_path[-1] == 'E'), "Last node is not E"
                            except AssertionError as ae:
                                #print >> sys.stderr, ae.message
                                pdb.set_trace()

                            yield new_path
        else:
            yield ['E']

