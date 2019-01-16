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


class EquivPrimerTargetsRegionPair(object):
    def __init__(self, olap_set_ID, equiv_region_fwd, equiv_region_rev):
        self.olap_set_ID = olap_set_ID
        self.ID = "%s%s" % (equiv_region_fwd.getID(), equiv_region_rev.getID())
        self.starts = None
        self.stops = None

        self.isoform_groups_dual_primed = None # tuple of isoform group(s)
        self.local_target_isoform_groups_dual_primed_by_ID = None
        self.local_nontarget_isoform_groups_dual_primed_by_ID = None
        
        self.equiv_region_fwd = equiv_region_fwd
        self.equiv_region_rev = equiv_region_rev

        self.default_primer_pairs = None

        # key is a tuple of isoform groups, and the value is a list of (fwd primer start, rev primer start) tuples in
        # genomic coordinates that will distinguish the igs with paired-end sequencing
        self.starts_for_igs = None # IMPORTANT: has a different format depending on the (boolean) value of self.regions_merge_or_abutt

        self.ranked_igs_tuples = None

        assert (equiv_region_fwd.getChromosome() == equiv_region_rev.getChromosome())
        self.chromosome = equiv_region_fwd.getChromosome()

        assert (equiv_region_fwd.getStrand() == equiv_region_rev.getStrand())
        self.strand = equiv_region_fwd.getStrand() # The strand of the target isoforms
        
        # 5p and 3p are relative to the strand of the target isoforms, not the fwd and rev primer regions in their strandwise context
        # HERE: need to confirm that this logic doesn't leave dangling starts/stops that are not recorded in self.starts/self.stops.
        if (self.strand=='+'):
            #most_5p_fwd_genomic_positions, most_3p_fwd_genomic_positions = self.equiv_region_fwd.getStartsStops() # HERE: Don't need these two
            #most_5p_rev_genomic_positions, most_3p_rev_genomic_positions = self.equiv_region_rev.getStartsStops()    
            self.starts = self.equiv_region_fwd.getStarts()
            self.stops = self.equiv_region_rev.getStops()
        else:
            #most_3p_fwd_genomic_positions, most_5p_fwd_genomic_positions = self.equiv_region_fwd.getStartsStops() # HERE: Don't need these two
            #most_3p_rev_genomic_positions, most_5p_rev_genomic_positions = self.equiv_region_rev.getStartsStops()    
            self.starts = self.equiv_region_rev.getStarts()
            self.stops = self.equiv_region_fwd.getStops()

        # The combination of fwd/rev primer region & +/- chromosome strand must now be considered together.
        # The 1st/2nd strand nomenclature is used here in the context of RT-qSeq, where the 1st strand is
        # the cDNA created by reverse transcription from an RNA template and the 2nd strand is the cDNA 
        # created from priming of the 1st strand cDNA.
        opposite_strand = '-' if (self.strand == '+') else '+'
        fwd_primer_2nd_strand_isoforms = equiv_region_fwd.getIsoforms(from_strand=self.strand)
        rev_primer_1st_strand_isoforms = equiv_region_rev.getIsoforms(from_strand=self.strand)

        # All isoforms with the fwd and rev primer, regardless of whether the resulting amplicon is a legal length or not.
        self.isoforms_dual_primed = fwd_primer_2nd_strand_isoforms.intersection(rev_primer_1st_strand_isoforms)
        assert (len(self.isoforms_dual_primed) > 0), "PPR does not have any dual-primed isoforms"
            
        # Values: False, "merge", "abutt" # TODO: abutt no longer used? DELETE.
        self.regions_merge_or_abutt = self.determineIfPrimerRegionsMergeOrAbutt(self.isoforms_dual_primed)

        # Only set if fwd/rev equiv regions are merged
        self.genomic_positions = None # Is a merge of the fwd/rev_equiv_region.genomic_positions if regions_merge_or_abutt isn't 'False'
        self.nuc_seq = None           
        self.local_index_to_genomic_positions = None
        self.isoform_genomic_paths = None

        if (self.regions_merge_or_abutt):
            self.unifyGenomicPositions()
            self.setNucSeq()


    def isDisjoint(self):
        return self.regions_merge_or_abutt==False


    def getRankedIGsTuples(self):
        return self.ranked_igs_tuples


    def determineIfPrimerRegionsMergeOrAbutt(self, isoforms_for_eval):
        regions_merge_or_abutt = False

        fwd_genomic_positions = self.equiv_region_fwd.getGenomicPositions(as_set=True)
        rev_genomic_positions = self.equiv_region_rev.getGenomicPositions(as_set=True)

        fwd_primer_5p_positions = self.equiv_region_fwd.getPrimer5pPositions()
        rev_primer_5p_positions = self.equiv_region_rev.getPrimer5pPositions()

        # TODO: Is it possible for two regions to both abutt and merge (some isoforms overlap and others abutt)?
        #if (not fwd_genomic_positions.isdisjoint(rev_genomic_positions)):         # TODO: delete?
        if (not fwd_primer_5p_positions.isdisjoint(rev_primer_5p_positions)):
            # Regions have the potential to be merged. Need to confirm that their overlap regions on each isoform are identical.
            overlaps_identical_on_isoform = []
            dont_overlap_on_isoform = []
            for isoform in isoforms_for_eval:
                fwd_start, fwd_stop, rev_start, rev_stop = self.getStartsStopsForIsoform(isoform)
                l,r = max(fwd_start, rev_start), min(fwd_stop, rev_stop)
                
                dont_overlap_on_isoform.append( fwd_start > rev_stop or fwd_stop < rev_start )

                isoform_genomic_positions = isoform.getAllGenomicCoords()

                fwd_shared_pos = fwd_genomic_positions & isoform_genomic_positions
                rev_shared_pos = rev_genomic_positions & isoform_genomic_positions
            
                fwd_shared_pos_lr = set(filter(lambda x: x>=l and x<=r, fwd_shared_pos))
                rev_shared_pos_lr = set(filter(lambda x: x>=l and x<=r, rev_shared_pos))

                overlaps_identical_on_isoform.append( len(fwd_shared_pos_lr) > 0 and fwd_shared_pos_lr == rev_shared_pos_lr )
                
            if (all(overlaps_identical_on_isoform) and not any(dont_overlap_on_isoform)):
                regions_merge_or_abutt = "merge"

        return regions_merge_or_abutt


    def unifyGenomicPositions(self):
        self.isoform_genomic_paths = {}
        self.local_index_to_genomic_positions = defaultdict(set)
        self.genomic_positions = nx.Graph()

        merged_lengths = set()
        trunc_merged_lengths = set()
        unique_trunc_merged_path_genomic_positions = set()

        for isoform in self.isoforms_dual_primed:
            fwd_genomic_positions = self.equiv_region_fwd.getIsoformSortedGenomicPositions(isoform)
            rev_genomic_positions = self.equiv_region_rev.getIsoformSortedGenomicPositions(isoform)

            l,r = (fwd_genomic_positions[0],rev_genomic_positions[-1]) if (self.strand=="+") else (rev_genomic_positions[0],fwd_genomic_positions[-1])

            merged_path = tuple(sorted(set(fwd_genomic_positions + rev_genomic_positions)))
            trunc_merged_path = tuple(filter(lambda x: x>=l and x<=r, merged_path))
            unique_trunc_merged_path_genomic_positions.add( trunc_merged_path )
            self.isoform_genomic_paths[isoform] = trunc_merged_path

            merged_lengths.add(len(merged_path))
            trunc_merged_lengths.add(len(trunc_merged_path))

            self.genomic_positions.add_edges_from(zip(merged_path[0:-1],merged_path[1:]))
            for pos in fwd_genomic_positions:
                self.genomic_positions.node[pos]["in_fwd_region"] = True
            for pos in rev_genomic_positions:
                self.genomic_positions.node[pos]["in_rev_region"] = True

        try:
            assert (len(merged_lengths)==1), "All merged paths not the same length"
            assert (len(trunc_merged_lengths)==1), "All truncated merged paths not the same length"
        except AssertionError as ae:
            print(ae.message, file=sys.stderr, flush=True)
            pdb.set_trace()
            
        trunc_merged_length = list(trunc_merged_lengths)[0]
        
        # Associate with each genomic position its local sequence index
        local_indices = range(trunc_merged_length) if (self.strand=='+') else range(trunc_merged_length-1,-1,-1)
        for path in unique_trunc_merged_path_genomic_positions:
            assert (len(local_indices) == len(path))
            for local_index, genomic_position in zip(local_indices, path):
                self.genomic_positions.node[genomic_position]["local_index"] = local_index
                self.local_index_to_genomic_positions[local_index].add(genomic_position)

        for n in self.genomic_positions.nodes():
            assert ("local_index" in self.genomic_positions.node[genomic_position] or
                    "in_fwd_region" in self.genomic_positions.node[genomic_position] or
                    "in_rev_region" in self.genomic_positions.node[genomic_position])


    def setNucSeq(self):
        isoform_subseqs_in_region = set()
        for isoform, sorted_positions in self.isoform_genomic_paths.items():
            extrema_start, extrema_stop = sorted_positions[0], sorted_positions[-1]
            if (isoform.getStrand() == "+"):
                extrema_seq = isoform.getSequence(None, extrema_start, extrema_stop)
                if (self.strand == "-"):
                    extrema_seq = extrema_seq[::-1].translate(DNA_complement_table)
            else:
                extrema_seq = isoform.getSequence(None, extrema_stop, extrema_start)
                if (self.strand == "+"):
                    extrema_seq = extrema_seq[::-1].translate(DNA_complement_table)
            isoform_subseqs_in_region.add(extrema_seq)

        assert (len(isoform_subseqs_in_region)==1), "Equiv region should only have one corresponding nucleotide sequence"
        self.nuc_seq = list(isoform_subseqs_in_region)[0] 
        

    def calcOptFwdRevPrimerPositions(self):
        assert(not self.isDisjoint()), "Method not intended for disjoint PPRs"

        half_amplicon_opt_len = int(designParams.amplicon_opt_len/2)

        if (self.equiv_region_fwd.getLength() < self.equiv_region_rev.getLength()):
            shorter_equiv_region, longer_equiv_region = (self.equiv_region_fwd, self.equiv_region_rev)
        else:
            shorter_equiv_region, longer_equiv_region = (self.equiv_region_rev, self.equiv_region_fwd)

        overlap_positions = list(map(itemgetter(0), filter(lambda nd: nd[1].get("in_fwd_region") and nd[1].get("in_rev_region"), self.genomic_positions.nodes(data=True))))
        overlap_positions.sort()
        olap_mid_genomic = overlap_positions[int(len(overlap_positions)/2)]

        shorter_local_pos_of_olap_mid = shorter_equiv_region.getLocalPosition(olap_mid_genomic)
        longer_local_pos_of_olap_mid = longer_equiv_region.getLocalPosition(olap_mid_genomic)

        shorter_is_fwd = shorter_equiv_region.getID()[0] == 'F'
        if (shorter_is_fwd):
            if (shorter_local_pos_of_olap_mid - half_amplicon_opt_len >= 0):
                shorter_local_ideal_pos = shorter_local_pos_of_olap_mid - half_amplicon_opt_len
                shorter_deficiency = 0
            else:
                shorter_local_ideal_pos = 0
                shorter_deficiency = abs(shorter_local_pos_of_olap_mid - half_amplicon_opt_len)

            if (longer_local_pos_of_olap_mid + half_amplicon_opt_len + shorter_deficiency < longer_equiv_region.getLength()):
                longer_local_ideal_pos = longer_local_pos_of_olap_mid + half_amplicon_opt_len + shorter_deficiency
            else:
                longer_local_ideal_pos = longer_equiv_region.getLength() - 1

        else:
            if (shorter_local_pos_of_olap_mid + half_amplicon_opt_len < shorter_equiv_region.getLength()):
                shorter_local_ideal_pos = shorter_local_pos_of_olap_mid + half_amplicon_opt_len
                shorter_deficiency = 0
            else:
                shorter_local_ideal_pos = shorter_equiv_region.getLength() - 1
                shorter_deficiency = shorter_local_pos_of_olap_mid + half_amplicon_opt_len - shorter_equiv_region.getLength()

            if (longer_local_pos_of_olap_mid - half_amplicon_opt_len - shorter_deficiency >= 0):
                longer_local_ideal_pos = longer_local_pos_of_olap_mid - half_amplicon_opt_len - shorter_deficiency
            else:
                longer_local_ideal_pos = 0
            
        shorter_pos_mapping = shorter_equiv_region.getMappingGenomicToLocalPositions()
        longer_pos_mapping = longer_equiv_region.getMappingGenomicToLocalPositions()

        shorter_dist_to_ideal = {}
        for genomic_pos, local_pos in shorter_pos_mapping.items():
            shorter_dist_to_ideal[genomic_pos] = abs(local_pos - shorter_local_ideal_pos)

        longer_dist_to_ideal = {}
        for genomic_pos, local_pos in longer_pos_mapping.items():
            longer_dist_to_ideal[genomic_pos] = abs(local_pos - longer_local_ideal_pos)

        if (shorter_is_fwd):
            return (shorter_dist_to_ideal, longer_dist_to_ideal)
        else:
            return (longer_dist_to_ideal, shorter_dist_to_ideal)


    def areAllowedPrimersForIsoformGroupsTuple(self, fwd_genomic_5p_pos_tup, fwd_primer_len, rev_genomic_5p_pos_tup, rev_primer_len):
        '''fwd_genomic_5p_pos_tup and rev_genomic_5p_pos_tup are tuples of potentially multiple start positions. Such cases correspond
        to the situation in which the same primer sequence corresponds to branching genomic positions within the same equiv region.'''
        all_is_allowed = []

        if (self.regions_merge_or_abutt):
            amplicon_lengths = set()
            for isoform in chain.from_iterable(self.isoform_groups_dual_primed):
                isoform_genomic_coords = isoform.getAllGenomicCoords()
                fwd_genomic_5p_pos = list(set(fwd_genomic_5p_pos_tup) & isoform_genomic_coords)
                rev_genomic_5p_pos = list(set(rev_genomic_5p_pos_tup) & isoform_genomic_coords)
                assert (len(fwd_genomic_5p_pos)==1 and len(rev_genomic_5p_pos)==1), "Should have only 1 start/stop"

                fwd_genomic_5p_pos = fwd_genomic_5p_pos[0]
                rev_genomic_5p_pos = rev_genomic_5p_pos[0]

                if ((isoform.getStrand() == "+" and fwd_genomic_5p_pos < rev_genomic_5p_pos) or
                    (isoform.getStrand() == "-" and fwd_genomic_5p_pos > rev_genomic_5p_pos)):
                    amplicon_len = isoform.getSequenceLength(fwd_genomic_5p_pos, rev_genomic_5p_pos)
                    amplicon_lengths.add(amplicon_len)
                    is_allowed = amplicon_len >= designParams.amplicon_min_len and amplicon_len <= designParams.amplicon_max_len
                    #if (not is_allowed): # DEBUG
                    #    if (amplicon_len >= designParams.amplicon_max_len):
                    #        print >> sys.stderr, "amplicon too long"
                    #    else:
                    #        print >> sys.stderr, "amplicon too short"
                else:
                    #print >> sys.stderr, "fwd primer upstream of rev primer" # DEBUG
                    is_allowed = False
                all_is_allowed.append( is_allowed )

            assert(len(amplicon_lengths)<=1), "Should only have at most one amplicon length for all isoforms for merge/abutt PPR"

        else:
            for fwd_genomic_5p_pos, rev_genomic_5p_pos in product(fwd_genomic_5p_pos_tup, rev_genomic_5p_pos_tup):
                fwd_rev_lens = self.starts_for_igs[self.isoform_groups_dual_primed].get((fwd_genomic_5p_pos, rev_genomic_5p_pos), None)
                pos_and_len_allowed_for_ig = fwd_rev_lens != None and fwd_primer_len in fwd_rev_lens[0] and rev_primer_len in fwd_rev_lens[1]
                if (pos_and_len_allowed_for_ig):
                    for isoform in chain.from_iterable(self.isoform_groups_dual_primed):
                        isoform_genomic_coords = isoform.getAllGenomicCoords()
                        fwd_genomic_5p_pos = list(set(fwd_genomic_5p_pos_tup) & isoform_genomic_coords)
                        rev_genomic_5p_pos = list(set(rev_genomic_5p_pos_tup) & isoform_genomic_coords)
                        assert (len(fwd_genomic_5p_pos)==1 and len(rev_genomic_5p_pos)==1), "Should have only 1 start/stop"

                        fwd_genomic_5p_pos = fwd_genomic_5p_pos[0]
                        rev_genomic_5p_pos = rev_genomic_5p_pos[0]

                        if ((isoform.getStrand() == "+" and fwd_genomic_5p_pos < rev_genomic_5p_pos) or
                            (isoform.getStrand() == "-" and fwd_genomic_5p_pos > rev_genomic_5p_pos)):
                            amplicon_len = isoform.getSequenceLength(fwd_genomic_5p_pos, rev_genomic_5p_pos)
                            is_allowed = amplicon_len >= designParams.amplicon_min_len and amplicon_len <= designParams.amplicon_max_len
                            if (not is_allowed): # DEBUG
                                if (amplicon_len >= designParams.amplicon_max_len):
                                    print("amplicon too long", file=sys.stderr, flush=True)
                                else:
                                    print("amplicon too short", file=sys.stderr, flush=True)
                        else:
                            #print >> sys.stderr, "fwd primer upstream of rev primer" # DEBUG
                            is_allowed = False
                        all_is_allowed.append( is_allowed )
                else:
                    #print >> sys.stderr, "primer 5p position combination not allowed" # DEBUG
                    all_is_allowed.append( False )
                    break
                
        #is_allowed = all(all_is_allowed)
        #assert (is_allowed or not any(all_is_allowed)), "The multiple fwd and rev primer genomic start positions should all yield the same product"

        return all(all_is_allowed)


    def getPrimerStartsForIsoformGroupsTuple(self):
        return self.starts_for_igs[self.isoform_groups_dual_primed]


    def getDecoupledPrimerSpecsForIsoformGroupsTuple(self, fwdrev):
        assert (fwdrev == "Fwd" or fwdrev == "Rev")
        primer_specs = []
        primer_starts_data = self.getPrimerStartsForIsoformGroupsTuple()

        if (self.regions_merge_or_abutt):
            primer_5p_pos_and_len = primer_starts_data[0] if (fwdrev == "Fwd") else primer_starts_data[1]
            for tup in primer_5p_pos_and_len:
                primer_specs.append( (tuple(map(itemgetter(0), tup[0])), tup[1]) )
        else:
            index = 0 if (fwdrev == "Fwd") else 1
            accum_primer_specs = defaultdict(set)
            for start_pos, primer_lens in primer_starts_data.items():
                accum_primer_specs[primer_lens[index]].add(start_pos[index])
            for primer_lens_tuple, primer_5p_positions in accum_primer_specs.items():
                primer_specs.append( (primer_5p_positions, primer_lens_tuple) )

        return primer_specs


    def getDecoupledPrimerStartsForIsoformGroupsTuple(self, fwdrev):
        assert (fwdrev == "Fwd" or fwdrev == "Rev")
        primer_starts_data = self.getPrimerStartsForIsoformGroupsTuple()

        if (self.regions_merge_or_abutt):
            if (fwdrev == "Fwd"):
                primer_starts = set(chain.from_iterable(map(itemgetter(0), map(itemgetter(0), primer_starts_data[0]))))
            else:
                primer_starts = set(chain.from_iterable(map(itemgetter(0), map(itemgetter(0), primer_starts_data[1]))))
        else:
            index = 0 if (fwdrev == "Fwd") else 1
            primer_starts = set(map(itemgetter(index), primer_starts_data.keys()))

        return primer_starts


    def setIsoformGroupingsPrimerPositions(self, genome_ref, target_isoform_groups):
        #possible_primer_lengths = range(designParams.min_primer_len, designParams.max_primer_len+1) # TODO: delete
        read_len = designParams.read_len

        if (self.regions_merge_or_abutt and len(self.nuc_seq) >= designParams.amplicon_min_len):
            fwd_primer_region_genomic_positions = self.equiv_region_fwd.getGenomicPositions(as_set=True)
            rev_primer_region_genomic_positions = self.equiv_region_rev.getGenomicPositions(as_set=True)

            # Trim the Fwd and Rev regions by most_pos_needed
            if (self.regions_merge_or_abutt == "abutt"):
                most_pos_needed = designParams.amplicon_max_len - designParams.min_primer_len
            else:
                most_pos_needed = designParams.amplicon_max_len + \
                                  len(fwd_primer_region_genomic_positions & rev_primer_region_genomic_positions) - \
                                  designParams.min_primer_len

            # Get the primer specs
            ig = tuple(sorted(self.isoforms_dual_primed, key=methodcaller("getCGDBName")))
            if (ig in target_isoform_groups):
                ig_by_ID = tuple(map(methodcaller("getCGDBName"), ig))
                self.starts_for_igs = {}
                self.starts_for_igs[(ig,)] = self.getDescriptorsForLegalPrimers(most_pos_needed) 
                self.ranked_igs_tuples = [((ig,), (ig_by_ID,), ())]

        elif (not self.regions_merge_or_abutt):
            self.starts_for_igs = defaultdict(dict)
            target_isoforms = set(chain.from_iterable(target_isoform_groups))
            
            # Constrain the primer genomic positions to only what is maximally potentially needed, based on read 
            # length, so that the whole fwd & rev regions aren't searched for primer start positions
            intervene_seq_lens = []
            common_isoforms_intervene_seq = []
            for isoform in self.isoforms_dual_primed:
                assert (self.strand == isoform.getStrand())
                fwd_start, fwd_stop, rev_start, rev_stop = self.getStartsStopsForIsoform(isoform)

                #isoform_genomic_coords = isoform.getAllGenomicCoords()
                if (self.strand == '+'):
                    #fwd_primer_region_termini = self.equiv_region_fwd.getStops(as_set=True) & isoform_genomic_coords
                    #rev_primer_region_termini = self.equiv_region_rev.getStarts(as_set=True) & isoform_genomic_coords
                    #most_3p_fwd_genomic_position, most_5p_rev_genomic_position = min(fwd_primer_region_termini), max(rev_primer_region_termini)
                    most_3p_fwd_genomic_position, most_5p_rev_genomic_position = fwd_stop, rev_start
                else:
                    #fwd_primer_region_termini = self.equiv_region_fwd.getStarts(as_set=True) & isoform_genomic_coords
                    #rev_primer_region_termini = self.equiv_region_rev.getStops(as_set=True) & isoform_genomic_coords
                    #most_3p_fwd_genomic_position, most_5p_rev_genomic_position = max(fwd_primer_region_termini), min(rev_primer_region_termini)
                    most_3p_fwd_genomic_position, most_5p_rev_genomic_position = fwd_start, rev_stop

                # Cache the dual primed isoforms and each of their associated intervening sequences.
                # There are special cases in which the fwd/rev equiv regions overlap on some isoforms.
                # The overlap amount is the negative of intervene_seq_len, and intervene_seq will be an empty string
                try:
                    intervene_seq_plus = isoform.getSequence(genome_ref, most_3p_fwd_genomic_position, most_5p_rev_genomic_position)
                    # Trim the 1bp of equiv_regions on each side
                    intervene_seq_len = len(intervene_seq_plus) - 2
                    intervene_seq = intervene_seq_plus[1:-1] if (len(intervene_seq_plus) > 2) else ""
                except AssertionError as ae:
                    if (str(ae) == "mRNA 5' position > 3' position"):
                        intervene_seq_len = isoform.getNegativeSequenceLength(most_3p_fwd_genomic_position, most_5p_rev_genomic_position)
                        intervene_seq = ''
                    else:
                        pdb.set_trace()
                        raise ae

                common_isoforms_intervene_seq.append( (isoform, intervene_seq_len, intervene_seq) )
                intervene_seq_lens.append(intervene_seq_len)

            shortest_intervene_seq_len = min(intervene_seq_lens) if (min(intervene_seq_lens) >= 0) else 0

            # IMPORTANT: The intervening sequence length constraint used here must be the same as is used in couldPotentiallyProductivelyPairWith()
            most_pos_needed = designParams.amplicon_max_len - shortest_intervene_seq_len - designParams.min_primer_len

            if (most_pos_needed > 0):
                fwd_nuc_seq = self.equiv_region_fwd.getNucSeq()
                rev_nuc_seq = self.equiv_region_rev.getNucSeq()

                fwd_primer_5p_positions_and_lens, rev_primer_5p_positions_and_lens = self.getDescriptorsForLegalPrimers(most_pos_needed)
                assert (fwd_primer_5p_positions_and_lens != None and rev_primer_5p_positions_and_lens != None)

                # Group the fwd primer/rev primer start positions that yield the same isoform groups (as revealed by paired-end sequencing)
                for (fwd_genomic_5ps, fwd_lens), (rev_genomic_5ps, rev_lens) in product(fwd_primer_5p_positions_and_lens, rev_primer_5p_positions_and_lens):
                    for ((fwd_genomic_5p, fwd_region_amplicon_len), (rev_genomic_5p, rev_region_amplicon_len)) in product(fwd_genomic_5ps, rev_genomic_5ps):
                        isoforms_grouped_by_paired_reads_short = defaultdict(list)
                        isoforms_grouped_by_paired_reads_long = defaultdict(list)
                        amplicon_lens = {}
                        for isoform, intervene_seq_len, intervene_seq in common_isoforms_intervene_seq:
                            if (intervene_seq_len < 0 and abs(intervene_seq_len) > fwd_region_amplicon_len):
                                continue

                            amplicon_len = fwd_region_amplicon_len + intervene_seq_len + rev_region_amplicon_len
                            #confirm_amplicon_len = isoform.getSequenceLength(fwd_genomic_5p, rev_genomic_5p)
                            #try:
                            #    assert (amplicon_len == confirm_amplicon_len), "Inconsistent amplicon lengths, computed vs actual"
                            #except AssertionError:
                            #    pdb.set_trace()
                                
                            if (amplicon_len < designParams.amplicon_min_len and isoform in target_isoforms):
                                # Disallow candidate primer positions to avoid possible confusion by actual sequencing of amplicons near min length limit.
                                break 

                            #if (amplicon_len > designParams.amplicon_max_len and isoform in target_isoforms):
                            #    pdb.set_trace() # Indicates a problem # Now don't think so. Such cases will be handled below with igs_short/igs_long
                                
                            if (amplicon_len <= designParams.confusion_amplicon_max_len): 
                                if (intervene_seq_len < 0):
                                    assert (intervene_seq == '')
                                    amplicon = fwd_nuc_seq[-fwd_region_amplicon_len:intervene_seq_len] + rev_nuc_seq[0:rev_region_amplicon_len]
                                else:
                                    amplicon = fwd_nuc_seq[-fwd_region_amplicon_len:] + intervene_seq + rev_nuc_seq[0:rev_region_amplicon_len]

                                paired_reads = (amplicon[0:read_len], amplicon[-read_len:])
                                isoforms_grouped_by_paired_reads_long[paired_reads].append( isoform )
                                if (amplicon_len <= designParams.amplicon_max_len):
                                    isoforms_grouped_by_paired_reads_short[paired_reads].append( isoform )
                                    amplicon_lens[isoform] = amplicon_len

                        igs_short = set()
                        for isoforms_w_same_read_pair in isoforms_grouped_by_paired_reads_short.values():
                            ig = sorted(isoforms_w_same_read_pair, key=methodcaller("getCGDBName"))
                            igs_short.add(tuple(ig))

                        igs_long = set()
                        for isoforms_w_same_read_pair in isoforms_grouped_by_paired_reads_long.values():
                            ig = sorted(isoforms_w_same_read_pair, key=methodcaller("getCGDBName"))
                            igs_long.add(tuple(ig))

                        # Checks to see if there is a target isoform group composed of (short) sequenceably-distinguishable amplicons that is
                        # robust to the presence of longer amplicons that could avoid size selection and get clustered and sequenced.
                        distinct_tigs = list(filter(lambda tig: tig in igs_short and tig in igs_long, target_isoform_groups))
                        if (len(distinct_tigs)>0):
                            mean_target_amplicon_len_deviation = np.mean(list(map(lambda i: abs(amplicon_lens[i]-designParams.amplicon_opt_len), chain.from_iterable(distinct_tigs))))
                            igs_long = sorted(igs_long)
                            self.starts_for_igs[tuple(igs_long)][(fwd_genomic_5p, rev_genomic_5p)] = (fwd_lens, rev_lens, int(mean_target_amplicon_len_deviation))

                numtigs_to_tigs = defaultdict(set)
                igs_tuples_w_rank_criteria = []
                for igs_tuple, primers_search_space in self.starts_for_igs.items():
                    ntigs_by_ID, tigs_by_ID = [], []
                    for ig in igs_tuple:
                        ig_by_ID = tuple(map(methodcaller("getCGDBName"), ig))
                        if (ig in target_isoform_groups):
                            tigs_by_ID.append(ig_by_ID)
                        else:
                            ntigs_by_ID.append(ig_by_ID)

                    if (len(tigs_by_ID) > 0):
                        igs_tuples_w_rank_criteria.append( (igs_tuple, tuple(tigs_by_ID), tuple(ntigs_by_ID), len(tigs_by_ID), len(ntigs_by_ID), len(primers_search_space)) )

                # The order in which to consider the igs_tuple to use for this PPR. Want to prioritize those that have the most TIGs.
                # igs_tuples with fewer TIGs will contain a subset of the TIGs of those of larger igs_tuples.
                #
                # REMOVED: Second ranking criterion: For igs_tuples with the same number of TIGs, want to use those with more non-TIGs because that implies
                # primers closer to edges of the equiv primer regions, and so increase the chances of picking up novel splicing. 
                #
                # Third ranking criterion: Larger search space implies close to edges of equiv primer regions.
                igs_tuples_w_rank_criteria.sort(key=itemgetter(3,5), reverse=True) # 2,
                self.ranked_igs_tuples = list(map(itemgetter(0,1,2), igs_tuples_w_rank_criteria))


    def getDescriptorsForLegalPrimers(self, most_pos_needed):
        '''Returns a list of tuples (fwd primer 5p genomic position, rev primer 5p genomic position, fwd primer lengths, rev primer lengths)
        specifying primer pairs that are contained within the allowed position ranges.'''

        fwd_primer_5p_positions_and_lens = self.equiv_region_fwd.getDescriptorsForLegalPrimers(most_pos_needed)
        rev_primer_5p_positions_and_lens = self.equiv_region_rev.getDescriptorsForLegalPrimers(most_pos_needed)

        return (fwd_primer_5p_positions_and_lens, rev_primer_5p_positions_and_lens)


    def getID(self):
        return self.ID


    def getFwdID(self):
        return self.equiv_region_fwd.getID()


    def getRevID(self):
        return self.equiv_region_rev.getID()


    def getFwdRevIDs(self):
        return (self.equiv_region_fwd.getID(), self.equiv_region_rev.getID())


    def getFwdEquivRegion(self):
        return self.equiv_region_fwd


    def getRevEquivRegion(self):
        return self.equiv_region_rev


    def getStartsStops(self):
        return (self.starts, self.stops)


    def getStrand(self):
        return self.strand


    def getPrimedIsoformGroupsTuple(self):
        return self.isoform_groups_dual_primed


    def myFwdOtherRevSharePos(self, other_ppr):
        my_fwd_genomic_positions = self.equiv_region_fwd.getGenomicPositions(as_set=True)
        other_rev_region = other_ppr.getRevEquivRegion()
        other_rev_genomic_positions = other_rev_region.getGenomicPositions(as_set=True)
        return not my_fwd_genomic_positions.isdisjoint(other_rev_genomic_positions)


    def myFwdOtherRevHaveCompatiblePrimers(self, other_ppr):
        # Incorporate code idea from constrainPrimersPerEquivRegion(). Can make evaluation on just the fwd/rev genomic 5' positions.
        sys.exit(1)


    def getStartsStopsForIsoform(self, isoform):
        fwd_start, fwd_stop = self.equiv_region_fwd.getIsoformStartStop(isoform)
        rev_start, rev_stop = self.equiv_region_rev.getIsoformStartStop(isoform)
        return fwd_start, fwd_stop, rev_start, rev_stop


    def getStartsStopsForIsoform_DEPRACATED(self, isoform):
        fwd_start, fwd_stop, rev_start, rev_stop = None, None, None, None

        isoform_genomic_positions = isoform.getAllGenomicCoords()

        fwd_genomic_positions_graph = self.equiv_region_fwd.getGenomicPositions()
        rev_genomic_positions_graph = self.equiv_region_rev.getGenomicPositions()

        # Get the start/stop that defines the contiguous genomic positions of the fwd/rev regions in this isoform
        fwd_starts, fwd_stops = self.equiv_region_fwd.getStartsStops()
        rev_starts, rev_stops = self.equiv_region_rev.getStartsStops()

        for start,stop in product(fwd_starts,fwd_stops):
            if (start in isoform_genomic_positions and stop in isoform_genomic_positions):
                path = nx.shortest_path(fwd_genomic_positions_graph,start,stop)
                if (isoform.areAllSequential(path)):
                    assert (fwd_start == None and fwd_stop == None), "There should only be one contiguous path"
                    fwd_start, fwd_stop = start, stop

        for start,stop in product(rev_starts,rev_stops):
            if (start in isoform_genomic_positions and stop in isoform_genomic_positions):
                path = nx.shortest_path(rev_genomic_positions_graph,start,stop)
                if (isoform.areAllSequential(path)):
                    assert (rev_start == None and rev_stop == None), "There should only be one contiguous path"
                    rev_start, rev_stop = start, stop

        termini_tup = (fwd_start, fwd_stop, rev_start, rev_stop)

        assert (all(map(lambda x: x != None, termini_tup))), "Start/stops for fwd/rev equiv primer regions were not all set"
        
        return termini_tup


    def doesConflictOnAnIsoformWith(self, other_ppr, isoforms_of_interest):
        '''Unless the Fwd or Rev primer regions are the same between self and other_ppr, there is a conflict if
        amplicons primed on a common isoform overlap.'''
        do_conflict = False
        isoform_w_conflict = None
        
        # If the Fwd this & Rev other OR Fwd other & Rev this have a common isoforms that is an isoform of interest,
        # then there is the possibility of unwanted crosspriming between the PPRs. The crosspriming does exist if
        # Fwd-Rev of this overlaps Fwd-Rev of other
        this_fwd_genomic_positions = self.equiv_region_fwd.getGenomicPositions(as_set=True)
        this_rev_genomic_positions = self.equiv_region_rev.getGenomicPositions(as_set=True)
        other_fwd_genomic_positions = other_ppr.equiv_region_fwd.getGenomicPositions(as_set=True)
        other_rev_genomic_positions = other_ppr.equiv_region_rev.getGenomicPositions(as_set=True)

        this_fwd_isoforms = self.equiv_region_fwd.getIsoforms(from_strand=self.strand)
        this_rev_isoforms = self.equiv_region_rev.getIsoforms(from_strand=self.strand)
        other_fwd_isoforms = other_ppr.equiv_region_fwd.getIsoforms(from_strand=self.strand)
        other_rev_isoforms = other_ppr.equiv_region_rev.getIsoforms(from_strand=self.strand)

        # Fwd this vs Rev other
        fwd_starts_this, fwd_stops_this = self.equiv_region_fwd.getStartsStops()
        rev_starts_other, rev_stops_other = other_ppr.equiv_region_rev.getStartsStops()
        for isoform in isoforms_of_interest & this_fwd_isoforms & other_rev_isoforms:
            isoform_genomic_coords = isoform.getAllGenomicCoords()
            if (self.strand=='+'):
                fwd_stop = min(fwd_stops_this & isoform_genomic_coords)
                rev_start = max(rev_starts_other & isoform_genomic_coords)

                seqlen = 0 if (fwd_stop >= rev_start) else isoform.getSequenceLength(fwd_stop, rev_start)
                if (seqlen > 0 and seqlen < designParams.confusion_amplicon_max_len):
                    do_conflict = True

            else:
                fwd_start = max(fwd_starts_this & isoform_genomic_coords)
                rev_stop = min(rev_stops_other & isoform_genomic_coords)

                seqlen = 0 if (fwd_start <= rev_stop) else isoform.getSequenceLength(fwd_start, rev_stop)
                if (seqlen > 0 and seqlen < designParams.confusion_amplicon_max_len):
                    do_conflict = True

            if (do_conflict):
                isoform_w_conflict = isoform
                break
                    
        if (not do_conflict):
            # Fwd other vs Rev this
            fwd_starts_other, fwd_stops_other = other_ppr.equiv_region_fwd.getStartsStops()
            rev_starts_this, rev_stops_this = self.equiv_region_rev.getStartsStops()
            for isoform in isoforms_of_interest & other_fwd_isoforms & this_rev_isoforms:
                isoform_genomic_coords = isoform.getAllGenomicCoords()
                if (self.strand=='+'):
                    fwd_stop = min(fwd_stops_other & isoform_genomic_coords)
                    rev_start = max(rev_starts_this & isoform_genomic_coords)

                    seqlen = 0 if (fwd_stop >= rev_start) else isoform.getSequenceLength(fwd_stop, rev_start)
                    if (seqlen > 0 and seqlen < designParams.confusion_amplicon_max_len):
                        do_conflict = True

                else:
                    fwd_start = max(fwd_starts_other & isoform_genomic_coords)
                    rev_stop = min(rev_stops_this & isoform_genomic_coords)

                    seqlen = 0 if (fwd_start <= rev_stop) else isoform.getSequenceLength(fwd_start, rev_stop)
                    if (seqlen > 0 and seqlen < designParams.confusion_amplicon_max_len):
                        do_conflict = True

                if (do_conflict):
                    isoform_w_conflict = isoform
                    break

        return (do_conflict, isoform_w_conflict)


    def couldAmplifyAnyTargetGroups(self, target_isoform_groups):
        '''Returns True if any isoform group in target_isoforms exactly matches a group of isoforms distinguishingly primed'''
        return self.starts_for_igs != None and len(self.starts_for_igs.keys()) > 0 # Because an igs tuple is only added to self.starts_for_igs if it contains a target isoform group


    def getNumOfftargetIsoforms(self, target_isoform_groups):
        return len(set(chain.from_iterable(self.isoform_groups_dual_primed)) - set(chain.from_iterable(target_isoform_groups)))


    def getSetOfPrimedIsoformGroups(self):
        return set(self.isoform_groups_dual_primed)


    def getIsoformGroupsDualPrimedByID(self):
        assert (self.isoform_groups_dual_primed != None)
        return set(map(lambda ig: tuple(map(lambda i: i.getCGDBName(), ig)), self.isoform_groups_dual_primed))


    def getIsoformGroupsDualPrimedIsoformsID(self):
        return set(map(methodcaller("getCGDBName"), chain.from_iterable(self.isoform_groups_dual_primed)))


    def getIsoformGroupsDualPrimed(self):
        return self.isoform_groups_dual_primed


    def getDualPrimedIsoformsIDs(self):
        return set(map(lambda i: i.getCGDBName(), self.isoforms_dual_primed))


    def designPairedPrimers(self, oligo_thermo, igs_tuple, ppr_unique_target_isoform_IDs, fwd_primer_cache, rev_primer_cache,
                            antitargets_fasta_name, local_transcriptome_fasta, tempdir):
        assert (len(self.starts_for_igs) != None), "Isoform groups and primer starts not set for PPR %s" % self.ID
        assert (len(self.starts_for_igs) > 0), "No isoform groups and primer starts for PPR %s" % self.ID
        
        if (self.regions_merge_or_abutt):
            primer3_global_args, primer3_seq_args = self.parameterizePrimer3ForIGsTuple_AbuttOrOverlap()
            primer_pairs = self.getPrimerPairForIGsTuple_AbuttOrOverlap(oligo_thermo, primer3_global_args, primer3_seq_args, fwd_primer_cache, rev_primer_cache,
                                                                        antitargets_fasta_name, local_transcriptome_fasta, ppr_unique_target_isoform_IDs, tempdir)
        else:
            primer_pairs = self.getPrimerPairForIGsTuple_DontAbuttOrOverlap(oligo_thermo, igs_tuple, fwd_primer_cache, rev_primer_cache, 
                                                                            antitargets_fasta_name, local_transcriptome_fasta, ppr_unique_target_isoform_IDs, tempdir)

        return primer_pairs


    def parameterizePrimer3ForIGsTuple_AbuttOrOverlap(self):
        # Set the primer3-py Global arguments
        primer3_global_args = {"PRIMER_TASK" : "generic",
                               "PRIMER_PICK_LEFT_PRIMER" : 1,
                               "PRIMER_PICK_INTERNAL_OLIGO" : 0,
                               "PRIMER_PICK_RIGHT_PRIMER" : 1,
                               "PRIMER_MIN_THREE_PRIME_DISTANCE" : 2,
                               "PRIMER_OPT_SIZE" : designParams.opt_primer_len,
                               "PRIMER_MIN_SIZE" : designParams.min_primer_len,
                               "PRIMER_MAX_SIZE" : designParams.max_primer_len,
                               "PRIMER_EXPLAIN_FLAG" : 1,
                               "PRIMER_NUM_RETURN" : 500,
                               "PRIMER_SALT_MONOVALENT" : designParams.monovalent_salt_conc,
                               "PRIMER_SALT_DIVALENT" : designParams.divalent_salt_conc,
                               "PRIMER_DNTP_CONC" : designParams.dntp_conc,
                               "PRIMER_DNA_CONC" : designParams.input_primer_conc_nM,
                               "PRIMER_MIN_TM" : designParams.min_primer_Tm,
                               "PRIMER_MAX_TM" : designParams.max_primer_Tm,
                               "PRIMER_OPT_TM" : designParams.opt_primer_Tm,
                               "PRIMER_TM_FORMULA" : 1,
                               "PRIMER_SALT_CORRECTIONS" : 1,
                               "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT" : 1,
                               "PRIMER_PRODUCT_SIZE_RANGE" : (designParams.amplicon_min_len, designParams.amplicon_max_len)}
                        
        if (len(self.nuc_seq) >= 10000):
            primer3_global_args["PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT"] = 0
            primer3_global_args["PRIMER_MAX_TEMPLATE_MISPRIMING"] = 12
        else:
            primer3_global_args["PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT"] = 1
            primer3_global_args["PRIMER_MAX_TEMPLATE_MISPRIMING_TH"] = designParams.min_primer_Tm - 8


        # Set the primer3-py Sequence arguments
        primerable_fwd_positions = list(n for n in self.genomic_positions if 'in_fwd_region' in self.genomic_positions.node[n] \
                                        and 'local_index' in self.genomic_positions.node[n])

        primerable_rev_positions = list(n for n in self.genomic_positions if 'in_rev_region' in self.genomic_positions.node[n] \
                                        and 'local_index' in self.genomic_positions.node[n])

        fwd_min_genomic_position, fwd_max_genomic_position = min(primerable_fwd_positions), max(primerable_fwd_positions)
        rev_min_genomic_position, rev_max_genomic_position = min(primerable_rev_positions), max(primerable_rev_positions)

        # The local start/stop are in the 5'->3' context of the nucleotide sequence from which primers are to be found
        if (self.strand == "+"):
            fwd_primer_local_start, fwd_primer_local_stop = self.genomic_positions.node[fwd_min_genomic_position]['local_index'], \
                                                            self.genomic_positions.node[fwd_max_genomic_position]['local_index']
            rev_primer_local_start, rev_primer_local_stop = self.genomic_positions.node[rev_min_genomic_position]['local_index'], \
                                                            self.genomic_positions.node[rev_max_genomic_position]['local_index']
        else:
            fwd_primer_local_start, fwd_primer_local_stop = self.genomic_positions.node[fwd_max_genomic_position]['local_index'], \
                                                            self.genomic_positions.node[fwd_min_genomic_position]['local_index']
            rev_primer_local_start, rev_primer_local_stop = self.genomic_positions.node[rev_max_genomic_position]['local_index'], \
                                                            self.genomic_positions.node[rev_min_genomic_position]['local_index']

        assert (fwd_primer_local_start <= rev_primer_local_start), "Fwd primer start position is upstream of Rev primer start position"

        fwd_primer_region_spec = (fwd_primer_local_start, fwd_primer_local_stop-fwd_primer_local_start+1)
        rev_primer_region_spec = (rev_primer_local_start, rev_primer_local_stop-rev_primer_local_start+1)

        primer3_seq_args = {"SEQUENCE_ID" : self.ID,
                            "SEQUENCE_TEMPLATE" : self.nuc_seq,
                            "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST" : fwd_primer_region_spec + rev_primer_region_spec}

        # If the forward and reverse primer regions are large enough, want to pick primers that are on both sides of the midpoint of their overlap.
        # Purpose is to get centered amplicons
        if (self.equiv_region_fwd.getLength() > 2*designParams.max_primer_len and self.equiv_region_rev.getLength() > 2*designParams.max_primer_len):
            genomic_overlap_positions = map(itemgetter(0), filter(lambda nd: nd[1].get("in_fwd_region") and nd[1].get("in_rev_region"), self.genomic_positions.nodes(data=True)))
            local_overlap_positions = list(map(lambda i: self.genomic_positions.node[i]["local_index"], genomic_overlap_positions))
            local_overlap_positions.sort()
            olap_mid_local = local_overlap_positions[int(len(local_overlap_positions)/2)]
            primer3_seq_args["SEQUENCE_TARGET"] = [(olap_mid_local,1)]

        return (primer3_global_args, primer3_seq_args)


    def getPrimerPairForIGsTuple_AbuttOrOverlap(self, oligo_thermo, primer3_global_args, primer3_seq_args, fwd_primer_cache, rev_primer_cache,
                                                antitargets_fasta_name, local_transcriptome_fasta, ppr_unique_target_isoform_IDs, tempdir):
        max_num_pairs_to_return = 25
        precandidate_pairs = []
        candidate_pairs = []
        primer_pairs = []
        fwd_erID, rev_erID = self.getFwdRevIDs()
        all_equiv_region_primers = {fwd_erID:set(), rev_erID:set()}
        AT = set(['a', 'A', 't', 'T'])
        GC = set(['c','C','g','G'])
        acgt = set(['a','c','g','t'])
        
        design_results = designPrimers(primer3_seq_args, primer3_global_args)
        if ("PRIMER_ERROR" in design_results):
            print("Primer3 Error: %s" % design_results["PRIMER_ERROR"], file=sys.stderr, flush=True)
            pdb.set_trace()
        else:
            dist_to_ideal_start_fwd, dist_to_ideal_start_rev = self.calcOptFwdRevPrimerPositions()

            num_pps_returned = design_results["PRIMER_PAIR_NUM_RETURNED"]
            for i in range(num_pps_returned):
                fwd_primer_5p_pos, fwd_primer_len = design_results["PRIMER_LEFT_%d" % i]
                rev_primer_5p_pos, rev_primer_len = design_results["PRIMER_RIGHT_%d" % i]
                alt_fwd_genomic_5p_pos = tuple(sorted(self.local_index_to_genomic_positions[fwd_primer_5p_pos]))
                alt_rev_genomic_5p_pos = tuple(sorted(self.local_index_to_genomic_positions[rev_primer_5p_pos]))

                if ((alt_fwd_genomic_5p_pos, fwd_primer_len) not in fwd_primer_cache):
                    fwd_seq = design_results["PRIMER_LEFT_%d_SEQUENCE" % i]
                    too_many_repeats = hasTooManyRepeats(fwd_seq)

                    fwd_Tm, fwd_frac_duplexed, fwd_has_hairpin, fwd_penalty = oligo_thermo.calcPrimerThermoDetails("fwd", fwd_seq)
                    num_GC = sum([x in GC for x in fwd_seq[-6:]])

                    if (num_GC == 6 or too_many_repeats or fwd_has_hairpin or fwd_penalty > 0.01 or fwd_Tm < designParams.min_primer_Tm or fwd_Tm > designParams.max_primer_Tm):
                        fwd_primer_cache[(alt_fwd_genomic_5p_pos, fwd_primer_len)] = None
                    else:
                        alt_fwd_genomic_positions = set.union(*map(lambda x: self.local_index_to_genomic_positions[x],
                                                                   range(fwd_primer_5p_pos, fwd_primer_5p_pos+fwd_primer_len)))
                        alt_fwd_genomic_positions = sorted(alt_fwd_genomic_positions)
                        nuc_imbalance = abs(3 - num_GC)
                        AT_3p = 1 if (fwd_seq[-1] in AT) else 0
                        aux_3p_penalty = nuc_imbalance + AT_3p
                        fwd_primer = PrimerSingleton(alt_fwd_genomic_5p_pos, fwd_seq, int(fwd_primer_len), fwd_Tm, fwd_frac_duplexed, tuple(alt_fwd_genomic_positions),
                                                     round(fwd_penalty,4), aux_3p_penalty)
                        fwd_primer_cache[(alt_fwd_genomic_5p_pos, fwd_primer_len)] = fwd_primer


                if ((alt_rev_genomic_5p_pos, rev_primer_len) not in rev_primer_cache):
                    rev_seq = design_results["PRIMER_RIGHT_%d_SEQUENCE" % i] 
                    too_many_repeats = hasTooManyRepeats(rev_seq)

                    rev_Tm, rev_frac_duplexed, rev_has_hairpin, rev_penalty = oligo_thermo.calcPrimerThermoDetails("rev", rev_seq)
                    num_GC = sum([x in GC for x in rev_seq[-6:]])
                    # Not currently used:
                    #(any(map(lambda s: s in rev_seq, designParams.disallowed_primer_subseqs)) or
                    # any(map(lambda s: rev_seq.endswith(s), designParams.disallowed_primer_3p_subseqs)))):
                    if (num_GC == 6 or too_many_repeats or rev_has_hairpin or rev_penalty > 0.01 or rev_Tm < designParams.min_primer_Tm or rev_Tm > designParams.max_primer_Tm):
                        rev_primer_cache[(alt_rev_genomic_5p_pos, rev_primer_len)] = None
                    else:
                        alt_rev_genomic_positions = set.union(*map(lambda x: self.local_index_to_genomic_positions[x],
                                                                   range(rev_primer_5p_pos-rev_primer_len+1, rev_primer_5p_pos+1)))
                        alt_rev_genomic_positions = sorted(alt_rev_genomic_positions)
                        nuc_imbalance = abs(3 - num_GC)
                        AT_3p = 1 if (rev_seq[-1] in AT) else 0
                        aux_3p_penalty = nuc_imbalance + AT_3p
                        rev_primer = PrimerSingleton(alt_rev_genomic_5p_pos, rev_seq, int(rev_primer_len), rev_Tm, rev_frac_duplexed, tuple(alt_rev_genomic_positions),
                                                     round(rev_penalty,4), aux_3p_penalty)
                        rev_primer_cache[(alt_rev_genomic_5p_pos, rev_primer_len)] = rev_primer

                fwd_primer = fwd_primer_cache[(alt_fwd_genomic_5p_pos, fwd_primer_len)]
                rev_primer = rev_primer_cache[(alt_rev_genomic_5p_pos, rev_primer_len)]

                if (fwd_primer != None and rev_primer != None):
                    precandidate_pairs.append( (fwd_primer,rev_primer) )
                    all_equiv_region_primers[rev_erID].add(rev_primer)
                    all_equiv_region_primers[fwd_erID].add(fwd_primer)

        # Associate expected isoforms to each primer
        expected_isoforms_per_primer_seq = defaultdict(set)
        fwd_equiv_region_isoform_IDs = self.equiv_region_fwd.getIsoformsIDs()
        for fwd_primer in all_equiv_region_primers[fwd_erID]:
            expected_isoforms_per_primer_seq[fwd_primer.seq].update(fwd_equiv_region_isoform_IDs)

        rev_equiv_region_isoform_IDs = self.equiv_region_rev.getIsoformsIDs()
        for rev_primer in all_equiv_region_primers[rev_erID]:
            expected_isoforms_per_primer_seq[rev_primer.seq].update(rev_equiv_region_isoform_IDs)

        print("Initial: %d fwd primers, %d rev primers" % (len(all_equiv_region_primers[fwd_erID]), len(all_equiv_region_primers[rev_erID])), file=sys.stderr)
        if (all(map(lambda v: len(v)>0, all_equiv_region_primers.values()))):
            filtered_equiv_region_primers, primer_num_nontarget_hybs = scanIndivPrimersAgainstLocalTranscriptome(oligo_thermo, expected_isoforms_per_primer_seq, all_equiv_region_primers,
                                                                                                                 antitargets_fasta_name, local_transcriptome_fasta,
                                                                                                                 ppr_unique_target_isoform_IDs, tempdir)
            if (filtered_equiv_region_primers != None):
                print("After antitarget and local transcriptome filtering: %d fwd primers, %d rev primers" %
                      (len(filtered_equiv_region_primers[fwd_erID]), len(filtered_equiv_region_primers[rev_erID])), file=sys.stderr)
                filtered_equiv_region_primers = scanIndivPrimersAgainstAdapterTags(filtered_equiv_region_primers, tempdir)
                print("After adapter tag filtering: %d fwd primers, %d rev primers" %
                      (len(filtered_equiv_region_primers[fwd_erID]), len(filtered_equiv_region_primers[rev_erID])), file=sys.stderr)
            else:
                print("After antitarget and local transcriptome filtering: 0 fwd primers or 0 rev primers", file=sys.stderr)

            if (filtered_equiv_region_primers != None):
                # For primers that pass local scan, compute aux data here and form complete candidate_pairs
                for fwd_primer,rev_primer in precandidate_pairs:
                    if (fwd_primer in filtered_equiv_region_primers[fwd_erID] and rev_primer in filtered_equiv_region_primers[rev_erID]):
                        num_lc = sum([x in acgt for x in fwd_primer.seq]) + sum([x in acgt for x in rev_primer.seq]) # Number of primer positions that are repeat masked
                        sum_individ_penalties = fwd_primer.thermo_penalty + rev_primer.thermo_penalty
                        sum_aux_penalties = fwd_primer.aux_3p_penalty + rev_primer.aux_3p_penalty
                        pair_penalty_term = oligo_thermo.calcPrimerPairThermoDetails(fwd_primer.seq, rev_primer.seq)
                        total_primer_pair_penalty = sum_individ_penalties + pair_penalty_term
                        cumul_dist_to_ideal = dist_to_ideal_start_fwd[fwd_primer[0][0]] + dist_to_ideal_start_rev[rev_primer[0][0]]
                        total_num_nontarget_hybs = primer_num_nontarget_hybs[fwd_primer] + primer_num_nontarget_hybs[rev_primer]
                        candidate_pairs.append( (fwd_primer, rev_primer, total_primer_pair_penalty, round(total_primer_pair_penalty,5), sum_aux_penalties,
                                                 fwd_primer.len+rev_primer.len, cumul_dist_to_ideal, total_num_nontarget_hybs, num_lc) )

        candidate_pairs = list(filter(lambda p: p[3] < 0.01, candidate_pairs))
        print("Compiled %d candidate pairs" % len(candidate_pairs), file=sys.stderr, flush=True)
        if (len(candidate_pairs) > 0):
            candidate_pairs.sort(key=itemgetter(8,3,7,4,6,5))
            annot_pps = list(map(lambda x: ((fwd_erID, x[0]), (rev_erID, x[1])), candidate_pairs))
            diverse_annot_pps = selectDiversityOfPrimers(annot_pps, max_num_pairs_to_return, 1)
        else:
            diverse_annot_pps = []

        print("Selected %d candidate pairs" % len(diverse_annot_pps), file=sys.stderr, flush=True)
        return list(map(lambda x: (x[0][1], x[1][1]), diverse_annot_pps))


    def getPrimerPairForIGsTuple_DontAbuttOrOverlap(self, oligo_thermo, igs_tuple, fwd_primer_cache, rev_primer_cache,
                                                    antitargets_fasta_name, local_transcriptome_fasta, ppr_unique_target_isoform_IDs, tempdir):
        max_num_pairs_to_return = 25
        candidate_pairs = []
        primer_pairs = []
        fwd_erID, rev_erID = self.getFwdRevIDs()
        all_equiv_region_primers = {fwd_erID:set(), rev_erID:set()}
        AT = set(['a', 'A', 't', 'T'])
        GC = set(['c','C','g','G'])
        acgt = set(['a','c','g','t'])

        # Compile the thermodynamics for each unique primer sequence
        unique_fwd_primer_pos = set()
        unique_rev_primer_pos = set()
        for ((fwd_5p_genomic_start, rev_5p_genomic_start), (fwd_primer_lens, rev_primer_lens, mean_target_amplicon_len_deviation)) in self.starts_for_igs[igs_tuple].items():
            unique_fwd_primer_pos.update( product([fwd_5p_genomic_start], fwd_primer_lens) )
            unique_rev_primer_pos.update( product([rev_5p_genomic_start], rev_primer_lens) )
        
        for fwd_5p_genomic_start, fwd_primer_len in unique_fwd_primer_pos:
            if ((fwd_5p_genomic_start, fwd_primer_len) not in fwd_primer_cache):
                fwd_seq = self.equiv_region_fwd.getSubseq(fwd_5p_genomic_start, fwd_primer_len)
                too_many_repeats = hasTooManyRepeats(fwd_seq)

                fwd_Tm, fwd_frac_duplexed, fwd_has_hairpin, fwd_penalty = oligo_thermo.calcPrimerThermoDetails("fwd", fwd_seq)
                num_GC = sum([x in GC for x in fwd_seq[-6:]])
                # Not currently used:
                #(any(map(lambda s: s in fwd_seq, designParams.disallowed_primer_subseqs)) or
                # any(map(lambda s: fwd_seq.endswith(s), designParams.disallowed_primer_3p_subseqs)))):
                if (num_GC == 6 or too_many_repeats or fwd_has_hairpin or fwd_penalty > 0.01 or fwd_Tm < designParams.min_primer_Tm or fwd_Tm > designParams.max_primer_Tm):
                    fwd_primer_cache[(fwd_5p_genomic_start, fwd_primer_len)] = None
                else:
                    alt_fwd_genomic_5p_pos, alt_fwd_genomic_positions = self.equiv_region_fwd.getAltGenomicPositions(fwd_5p_genomic_start, fwd_primer_len)
                    nuc_imbalance = abs(3 - num_GC)
                    AT_3p = 1 if (fwd_seq[-1] in AT) else 0
                    aux_3p_penalty = nuc_imbalance + AT_3p
                    fwd_primer = PrimerSingleton(alt_fwd_genomic_5p_pos, fwd_seq, fwd_primer_len, fwd_Tm, fwd_frac_duplexed, alt_fwd_genomic_positions,
                                                 round(fwd_penalty,4), aux_3p_penalty)
                    fwd_primer_cache[(fwd_5p_genomic_start, fwd_primer_len)] = fwd_primer
                    all_equiv_region_primers[fwd_erID].add( fwd_primer )

            elif (fwd_primer_cache[(fwd_5p_genomic_start, fwd_primer_len)] != None):
                all_equiv_region_primers[fwd_erID].add( fwd_primer_cache[(fwd_5p_genomic_start, fwd_primer_len)] )

        for rev_5p_genomic_start, rev_primer_len in unique_rev_primer_pos:
            if ((rev_5p_genomic_start, rev_primer_len) not in rev_primer_cache):
                rev_seq = self.equiv_region_rev.getSubseq(rev_5p_genomic_start, rev_primer_len)
                too_many_repeats = hasTooManyRepeats(rev_seq)

                rev_Tm, rev_frac_duplexed, rev_has_hairpin, rev_penalty = oligo_thermo.calcPrimerThermoDetails("rev", rev_seq)
                num_GC = sum([x in GC for x in rev_seq[-6:]])
                # Not currently used:
                #(any(map(lambda s: s in rev_seq, designParams.disallowed_primer_subseqs)) or
                #     any(map(lambda s: rev_seq.endswith(s), designParams.disallowed_primer_3p_subseqs)))):
                if (num_GC == 6 or too_many_repeats or rev_has_hairpin or rev_penalty > 0.01 or rev_Tm < designParams.min_primer_Tm or rev_Tm > designParams.max_primer_Tm):
                    rev_primer_cache[(rev_5p_genomic_start, rev_primer_len)] = None
                else:
                    alt_rev_genomic_5p_pos, alt_rev_genomic_positions = self.equiv_region_rev.getAltGenomicPositions(rev_5p_genomic_start, rev_primer_len)
                    nuc_imbalance = abs(3 - num_GC)
                    AT_3p = 1 if (rev_seq[-1] in AT) else 0
                    aux_3p_penalty = nuc_imbalance + AT_3p
                    rev_primer = PrimerSingleton(alt_rev_genomic_5p_pos, rev_seq, rev_primer_len, rev_Tm, rev_frac_duplexed, alt_rev_genomic_positions,
                                                 round(rev_penalty,4), aux_3p_penalty)
                    rev_primer_cache[(rev_5p_genomic_start, rev_primer_len)] = rev_primer
                    all_equiv_region_primers[rev_erID].add( rev_primer )

            elif (rev_primer_cache[(rev_5p_genomic_start, rev_primer_len)] != None):
                all_equiv_region_primers[rev_erID].add( rev_primer_cache[(rev_5p_genomic_start, rev_primer_len)] )

        # Associate expected isoforms to each primer
        expected_isoforms_per_primer_seq = defaultdict(set)
        fwd_equiv_region_isoform_IDs = self.equiv_region_fwd.getIsoformsIDs()
        for fwd_primer in all_equiv_region_primers[fwd_erID]:
            expected_isoforms_per_primer_seq[fwd_primer.seq].update(fwd_equiv_region_isoform_IDs)

        rev_equiv_region_isoform_IDs = self.equiv_region_rev.getIsoformsIDs()
        for rev_primer in all_equiv_region_primers[rev_erID]:
            expected_isoforms_per_primer_seq[rev_primer.seq].update(rev_equiv_region_isoform_IDs)

        print("Initial: %d fwd primers, %d rev primers" % (len(all_equiv_region_primers[fwd_erID]), len(all_equiv_region_primers[rev_erID])), file=sys.stderr)
        if (all(map(lambda v: len(v)>0, all_equiv_region_primers.values()))):
            filtered_equiv_region_primers, primer_num_nontarget_hybs = scanIndivPrimersAgainstLocalTranscriptome(oligo_thermo, expected_isoforms_per_primer_seq,
                                                                                                                 all_equiv_region_primers, antitargets_fasta_name,
                                                                                                                 local_transcriptome_fasta, ppr_unique_target_isoform_IDs, tempdir)
            if (filtered_equiv_region_primers != None):
                print("After antitarget and local transcriptome filtering: %d fwd primers, %d rev primers" %
                      (len(filtered_equiv_region_primers[fwd_erID]), len(filtered_equiv_region_primers[rev_erID])), file=sys.stderr)
                filtered_equiv_region_primers = scanIndivPrimersAgainstAdapterTags(filtered_equiv_region_primers, tempdir)
                if (filtered_equiv_region_primers != None):
                    print("After adapter tag filtering: %d fwd primers, %d rev primers" %
                          (len(filtered_equiv_region_primers[fwd_erID]), len(filtered_equiv_region_primers[rev_erID])), file=sys.stderr)
                else:
                    print("After antitarget and local transcriptome filtering: 0 fwd primers or 0 rev primers", file=sys.stderr)

            else:
                print("After antitarget and local transcriptome filtering: 0 fwd primers or 0 rev primers", file=sys.stderr)
                
            if (filtered_equiv_region_primers != None):
                print("Begin compiling candidate pairs", file=sys.stderr, flush=True)
                for ((fwd_5p_genomic_start, rev_5p_genomic_start), (fwd_primer_lens, rev_primer_lens, mean_target_amplicon_len_deviation)) in self.starts_for_igs[igs_tuple].items():
                    for fwd_primer_len in fwd_primer_lens:
                        fwd_primer = fwd_primer_cache[(fwd_5p_genomic_start, fwd_primer_len)]
                        if (fwd_primer != None and fwd_primer in filtered_equiv_region_primers[fwd_erID]):
                            num_fwd_lc = sum([x in acgt for x in fwd_primer.seq])
                            for rev_primer_len in rev_primer_lens:
                                rev_primer = rev_primer_cache[(rev_5p_genomic_start, rev_primer_len)]
                                if (rev_primer != None and rev_primer in filtered_equiv_region_primers[rev_erID]):
                                    num_lc = num_fwd_lc + sum([x in acgt for x in rev_primer.seq]) # Number of primer positions that are repeat masked
                                    sum_indiv_penalties = fwd_primer.thermo_penalty + rev_primer.thermo_penalty
                                    sum_aux_penalties = fwd_primer.aux_3p_penalty + rev_primer.aux_3p_penalty
                                    pair_penalty_term = oligo_thermo.calcPrimerPairThermoDetails(fwd_primer.seq, rev_primer.seq)
                                    total_primer_pair_penalty = sum_indiv_penalties + pair_penalty_term
                                    total_num_nontarget_hybs = primer_num_nontarget_hybs[fwd_primer] + primer_num_nontarget_hybs[rev_primer]
                                    candidate_pairs.append((fwd_primer, rev_primer, total_primer_pair_penalty, round(total_primer_pair_penalty,5), sum_aux_penalties,
                                                            fwd_primer.len+rev_primer.len, mean_target_amplicon_len_deviation, total_num_nontarget_hybs, num_lc))

        candidate_pairs = list(filter(lambda p: p[3] < 0.01, candidate_pairs))
        print("Compiled %d candidate pairs" % len(candidate_pairs), file=sys.stderr, flush=True)

        if (len(candidate_pairs) > 0):
            candidate_pairs.sort(key=itemgetter(8,3,7,4,6,5))
            annot_pps = list(map(lambda x: ((fwd_erID, x[0]), (rev_erID, x[1])), candidate_pairs))
            diverse_annot_pps = selectDiversityOfPrimers(annot_pps, max_num_pairs_to_return, 1)
        else:
            diverse_annot_pps = []

        # TODO: apply
        # Restrict search to those primers approx 50 fwd x 50 rev primer start pos combinations yielding the lowest deviations
        #deviations = []
        #for ((fwd_5p_genomic_start, rev_5p_genomic_start), (fwd_primer_lens, rev_primer_lens, mean_target_amplicon_len_deviation)) in self.starts_for_igs[igs_tuple].items():
        #    deviations.append(mean_target_amplicon_len_deviation)
        #deviations.sort()
        #deviation_cutoff = deviations[2499] if (len(deviations)>2500) else deviations[-1]
        #if (mean_target_amplicon_len_deviation <= deviation_cutoff):

        print("Selected %d candidate pairs" % len(diverse_annot_pps), file=sys.stderr, flush=True)
        return list(map(lambda x: (x[0][1], x[1][1]), diverse_annot_pps))

    
    def hasDefaultPrimerPairs(self):
        return self.default_primer_pairs != None


    def clearDefaultPrimerPairs(self):
        assert (self.isoform_groups_dual_primed != None and self.default_primer_pairs != None), "Clearing already-cleared default primer pairs"
        self.isoform_groups_dual_primed = None
        self.default_primer_pairs = None


    def setDefaultPrimerPairs(self, igs_tuple, tigs_by_ID, ntigs_by_ID, primer_pairs):
        '''primer_pairs is a list of tuples (primer_pair, mean_target_product_len_deviation, new_igs_tuples_by_ID)'''
        assert (isinstance(igs_tuple, tuple)), "igs_tuple is not a tuple"
        assert (all(map(lambda t: isinstance(t,tuple), igs_tuple))),  "igs_tuple contains an element that is not a tuple"

        self.isoform_groups_dual_primed = igs_tuple
        self.local_target_isoform_groups_dual_primed_by_ID = tigs_by_ID
        self.local_nontarget_isoform_groups_dual_primed_by_ID = ntigs_by_ID
        
        if (self.default_primer_pairs == None):
            self.default_primer_pairs = {igs_tuple:[]}
        else:
            assert (igs_tuple not in self.default_primer_pairs), "Default primer pair already set for isoform groups tuple"

        fwd_rev_IDs = self.getFwdRevIDs()
        for counter, pp in enumerate(primer_pairs):
            new_sop = SetOfPrimers(self.chromosome, self.olap_set_ID, self.ID, counter, zip(fwd_rev_IDs, pp[0:2]))
            new_sop.setIsoformGroupsFromCategorized(tigs_by_ID, ntigs_by_ID, pp[-1])  # pp[-1] is new_igs_tuples_by_ID
            self.default_primer_pairs[igs_tuple].append(new_sop)


    def getDefaultPrimerPairs(self):
        return self.default_primer_pairs[self.isoform_groups_dual_primed]


    def getIsoformsIDsFormatted(self):
        ret_line = ""

        if (self.isoform_groups_dual_primed != None):
            line_elems = []
            for ig in self.isoform_groups_dual_primed:
                line_elems.append( "(%s)" % ",".join(map(lambda x:x.getCGDBName(), ig)) )
            isoform_ID_groups = ", ".join(line_elems)
            ret_line += "%s%s Fwd+Rev: %s\n" % (self.olap_set_ID, self.ID, isoform_ID_groups)

        isoform_names = sorted(map(lambda x:x.getCGDBName(), self.isoforms_primed_only_2nd_strand))
        if (isinstance(isoform_names, list)):
            isoform_names = ", ".join( isoform_names )
        ret_line += "%s%s Only 2nd strand: %s\n" % (self.olap_set_ID, self.ID, isoform_names)

        isoform_names = sorted(map(lambda x:x.getCGDBName(), self.isoforms_primed_only_1st_strand))
        if (isinstance(isoform_names, list)):
            isoform_names = ", ".join( isoform_names )
        ret_line += "%s%s Only 1st strand: %s\n" % (self.olap_set_ID, self.ID, isoform_names)

        return ret_line


    def getAsBED12Line(self, label_prefix, score=0, itemRgb="0,0,0"):
        all_genomic_positions = list(self.equiv_region_fwd.getGenomicPositions(as_set=True) | self.equiv_region_rev.getGenomicPositions(as_set=True))
        all_genomic_positions.sort()
        start = all_genomic_positions[0]
        stop = all_genomic_positions[-1]
        
        block_sizes = []
        block_starts = []
        for k, g in groupby(enumerate(all_genomic_positions), lambda ix:ix[0]-ix[1]):
            group = list(map(itemgetter(1), g))
            block_starts.append( str(group[0] - start) )
            block_sizes.append( str(group[-1] - group[0] + 1) )
        
        label = "%s%s" % (label_prefix, self.ID)
        bed12_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s" % \
          (self.chromosome, start-1, stop, label, score, ".", start, start, itemRgb, len(block_starts), ','.join(block_sizes), ','.join(block_starts))

        return bed12_line
