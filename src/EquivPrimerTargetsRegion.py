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


# 1) Contiguous set(s) of positions that
# 2) have one nucleotide sequence
# 3) contains the primer 5' positions + primer lengths that are consistent with 1)-2)
class EquivPrimerTargetsRegion(object):
    def __init__(self, olap_set_ID, ID, chromosome, strand, sorted_primer_positions_per_isoform,
                 nuc_seq, isoforms_signature, rev_fwd, opposite_overlapping_isoforms): # TODO DELETE equiv_region_position_path_graph, 
        assert (rev_fwd in ["rev", "fwd"])
        self.olap_set_ID = olap_set_ID
        self.ID = ID
        #self.genomic_positions = equiv_region_position_path_graph # TODO: don't really need. Replace with self.genomic_position_to_local_index = {} DELETE
        self.genomic_position_to_local_index = {}
        self.sorted_primer_positions_per_isoform = sorted_primer_positions_per_isoform
        self.length = None
        self.starts = None
        self.stops = None
        self.rev_fwd = rev_fwd
        self.isoforms_signature = set(isoforms_signature)
        self.opposite_overlapping_isoforms = opposite_overlapping_isoforms
        self.chromosome = chromosome
        self.strand = strand
        self.nuc_seq = nuc_seq  # In the 5'->3' orientation according to self.strand (i.e. reverse complemented if self.strand == '-')
        self.primer_5p_positions_and_lens = None
        self.local_index_to_genomic_positions = defaultdict(set)

        self.length = len(nuc_seq)
        local_indices = list(range(self.length)) if (self.strand=='+') else list(range(self.length-1,-1,-1))

        self.starts = set()
        self.stops = set()
        self.start_stop_per_isoform = {}
        self.unique_sorted_genomic_positions = set()
        for isoform, sorted_genomic_positions in sorted_primer_positions_per_isoform.items():
            self.unique_sorted_genomic_positions.add( sorted_genomic_positions )
            start, stop = sorted_genomic_positions[0], sorted_genomic_positions[-1]
            self.starts.add(start)
            self.stops.add(stop)
            self.start_stop_per_isoform[isoform] = (start, stop)
            
            # Associate genomic positions with local sequence indices.
            # TODO: Implementation has redundancy, could be improved
            assert (len(local_indices) == len(sorted_genomic_positions))
            for local_index, genomic_position in zip(local_indices, sorted_genomic_positions):
                self.genomic_position_to_local_index[genomic_position] = local_index
                #self.genomic_positions.node[genomic_position]["local_index"] = local_index DELETE
                self.local_index_to_genomic_positions[local_index].add(genomic_position)


    def getAltGenomicPositions(self, primer_5p_genomic_position, primer_len):
        local_5p_index = self.genomic_position_to_local_index[primer_5p_genomic_position]
        alt_genomic_5p_pos = sorted(self.local_index_to_genomic_positions[local_5p_index])
        
        if (self.rev_fwd == "fwd"):
            alt_genomic_pos = set.union(*map(lambda x: self.local_index_to_genomic_positions[x], range(local_5p_index, local_5p_index+primer_len)))
        else:
            alt_genomic_pos = set.union(*map(lambda x: self.local_index_to_genomic_positions[x], range(local_5p_index-primer_len+1, local_5p_index+1)))
        alt_genomic_pos = sorted(alt_genomic_pos)

        return (tuple(alt_genomic_5p_pos), tuple(alt_genomic_pos))


    def getID(self):
        return self.ID


    def setID(self, ID):
        self.ID = ID


    def getLength(self):
        return self.length


    def getIsoforms(self, from_strand):
        if (from_strand == self.strand):
            return self.isoforms_signature
        else:
            return self.opposite_overlapping_isoforms


    def getIsoformsIDs(self):
        return set(map(lambda x:x.getCGDBName(), self.isoforms_signature))
        

    def getIsoformsIDsFormatted(self):
        isoform_names = list(map(lambda x:x.getCGDBName(), self.isoforms_signature))
        if (isinstance(isoform_names, list)):
            isoform_names = ", ".join( isoform_names )
        return "%s%s: %s" % (self.olap_set_ID, self.ID, isoform_names)


    def getLocalPosition(self, genomic_position):
        return self.genomic_position_to_local_index[genomic_position]


    def getGenomicPositions(self, as_set=False):
        if (as_set):
            return set(self.genomic_position_to_local_index.keys())
        else:
            return list(self.genomic_position_to_local_index.keys())


    def getMappingGenomicToLocalPositions(self):
        return self.genomic_position_to_local_index

        
    def getChromosome(self):
        return self.chromosome

    
    def getStrand(self):
        return self.strand


    def getType(self):
        return self.rev_fwd
    

    def getStarts(self):
        return self.starts


    def getStops(self):
        return self.stops            


    def getStartsStops(self):
        return (self.starts, self.stops)


    def getExtremalStartStop(self):
        return (min(self.starts), max(self.stops))


    def getIsoformStartStop(self, isoform):
        return self.start_stop_per_isoform[isoform]


    def getIsoformSortedGenomicPositions(self, isoform):
        return self.sorted_primer_positions_per_isoform[isoform]


    def hasPositionDependentPrimerData(self):
        return len(self.primer_5p_positions_and_lens) > 0


    def setPositionDependentPrimerData(self, possible_primer_lengths):
        '''For each genomic position, set the allowed primer lengths and its distance into the equiv region'''
        self.primer_5p_positions_and_lens = []
        
        possible_primer_lengths.sort(reverse=True)
        possible_primer_lengths = tuple(possible_primer_lengths)
        min_primer_len = possible_primer_lengths[-1]

        if ((self.strand == "+" and self.rev_fwd == "fwd") or (self.strand == "-" and self.rev_fwd == "rev")):
            amount_in_region = range(self.length,0,-1)
            S = set()
            for path in self.unique_sorted_genomic_positions:
                S.update(zip(path, amount_in_region))

        elif ((self.strand == "+" and self.rev_fwd == "rev") or (self.strand == "-" and self.rev_fwd == "fwd")):
            amount_in_region = range(1,self.length+1)
            S = set()
            for path in self.unique_sorted_genomic_positions:
                S.update(zip(path, amount_in_region))

        else:
            print("ERROR: encountered illegal case in setPositionDependentPrimerData(). Exiting.", file=sys.stderr, flush=True)
            sys.exit(1)

        assert (len(S) > 0)

        q = self.length + 1
        for m, primer_len in enumerate(possible_primer_lengths):
            self.primer_5p_positions_and_lens.append( (tuple(filter(lambda x: x[1]>=primer_len and x[1]<q, S)), possible_primer_lengths[m:]) )
            q = primer_len


    def getPrimer5pPositions(self):
        primer_5p_positions = set()
        for pos_and_amount in self.primer_5p_positions_and_lens:
            primer_5p_positions.update( set(map(itemgetter(0), pos_and_amount[0])) )
        return primer_5p_positions


    def getDescriptorsForLegalPrimers(self, most_pos_needed):
        if (most_pos_needed >= self.length):
            genomic_pos_and_amount = self.primer_5p_positions_and_lens
        else:
            cutoff = self.length - most_pos_needed
            genomic_pos_and_amount = []
            for pos_and_amount, primer_lengths in self.primer_5p_positions_and_lens:
                genomic_pos_and_amount.append( (tuple(filter(lambda x: x[1] <= most_pos_needed, pos_and_amount)), primer_lengths) )
        return genomic_pos_and_amount


    def getNucSeq(self):
        assert (self.nuc_seq != None and len(self.nuc_seq)>0)
        return self.nuc_seq


    def getSubseq(self, genome_5p_pos, subseq_len):
        local_5p_pos = self.genomic_position_to_local_index[genome_5p_pos]
        if (self.rev_fwd == "fwd"):
            assert (local_5p_pos + subseq_len <= len(self.nuc_seq)), "Error calculation. Primer extends off end of template."
            primer_seq = self.nuc_seq[local_5p_pos:local_5p_pos+subseq_len]
        else:
            assert (local_5p_pos - subseq_len + 1 >= 0), "Error calculation. Primer extends off beginning of template."
            primer_seq = self.nuc_seq[local_5p_pos-subseq_len+1:local_5p_pos+1][::-1].translate(DNA_complement_table)

        return primer_seq


    def hasNucSeq(self):
        return self.nuc_seq != None


    #def setNucSeq(self, genome_ref):
    #    '''Sets the nucleotide sequence 5'->3' for the appropriate strand.'''
    #    assert (self.nuc_seq == None), "Nucleotide sequence already set"

    #    # Any start/stop pair will suffice. All should yield the same nucleotide sequence
    #    paths = map(lambda p: (len(p),p), nx.all_simple_paths(self.genomic_positions, self.starts[0], self.stops[0]))
    #    paths.sort(key=itemgetter(0), reverse=True)
    #    longest_path = paths[0][1]
        
    #    region_tuples = []
    #    for k, g in groupby(enumerate(longest_path), lambda (i,x):i-x):
    #        group = map(itemgetter(1), g)
    #        region_tuples.append( (group[0], group[-1]) )

    #    # Cast to int (from uint) prevents pyfaidx AssertionError
    #    nuc_seq = str("".join(map(lambda t: genome_ref[self.chromosome][int(t[0])-1:int(t[1])], region_tuples)))

    #    if (self.strand == '-'):
    #        self.nuc_seq = nuc_seq[::-1].translate(DNA_complement_table)
    #    else:
    #        self.nuc_seq = str(nuc_seq)


    def couldPotentiallyProductivelyPairWith(self, equiv_primer_region, target_isoforms):
        can_pair = False
        assert (self.rev_fwd != equiv_primer_region.getType())
        assert (self.chromosome == equiv_primer_region.getChromosome())
        assert (self.strand == equiv_primer_region.getStrand())
        
        common_target_isoforms = self.isoforms_signature & equiv_primer_region.getIsoforms(from_strand=self.strand) & target_isoforms
        both_regions_have_primers = self.hasPositionDependentPrimerData() and equiv_primer_region.hasPositionDependentPrimerData()

        if (len(common_target_isoforms) > 0 and both_regions_have_primers):
            fwd_genomic_positions = self.getGenomicPositions(as_set=True) if (self.rev_fwd == "fwd") else equiv_primer_region.getGenomicPositions(as_set=True)
            rev_genomic_positions = equiv_primer_region.getGenomicPositions(as_set=True) if (self.rev_fwd == "fwd") else self.getGenomicPositions(as_set=True)
            fwd_starts, fwd_stops = self.getStartsStops() if (self.rev_fwd == "fwd") else equiv_primer_region.getStartsStops()
            rev_starts, rev_stops = equiv_primer_region.getStartsStops() if (self.rev_fwd == "fwd") else self.getStartsStops()
            
            if (len(fwd_genomic_positions & rev_genomic_positions) > 0):
                # Find the number of genomic positions in the combined fwd and rev genomic regions. Calculation logic needs to account for
                # regions overlapping in such a way that some positions of one or both regions cannot be used for primers (eg Fwd region start
                # after Rev region stop on + strand prevents Rev region positions that are before the Fwd region start from being used for primers).

                if (self.strand == '+'):
                    fwd_starts_in_rev = rev_genomic_positions & fwd_starts
                    allowed_rev_positions = set([x for x in rev_genomic_positions if x>=min(fwd_starts_in_rev)]) if (len(fwd_starts_in_rev)>0) else rev_genomic_positions

                    rev_stops_in_fwd = fwd_genomic_positions & rev_stops
                    allowed_fwd_positions = set([x for x in fwd_genomic_positions if x<=max(rev_stops_in_fwd)]) if (len(rev_stops_in_fwd)>0) else fwd_genomic_positions
                else:
                    fwd_stops_in_rev = rev_genomic_positions & fwd_stops
                    allowed_rev_positions = set([x for x in rev_genomic_positions if x<=max(fwd_stops_in_rev)]) if (len(fwd_stops_in_rev)>0) else rev_genomic_positions

                    rev_starts_in_fwd = fwd_genomic_positions & rev_starts
                    allowed_fwd_positions = set([x for x in fwd_genomic_positions if x>=min(rev_starts_in_fwd)]) if (len(rev_starts_in_fwd)>0) else fwd_genomic_positions

                can_pair = len(allowed_fwd_positions | allowed_rev_positions) >= designParams.amplicon_min_len
            else:
                # Confirm that there is an intervening sequence that fits length constraints
                max_intervene_len = designParams.amplicon_max_len - 2*designParams.min_primer_len
                for isoform in common_target_isoforms:
                    genome_coords = isoform.getAllGenomicCoords()
                    if (self.strand == '+'):
                        in_isoform = filter(lambda fr: fr[0]<fr[1] and fr[0] in genome_coords and fr[1] in genome_coords, product(fwd_stops, rev_starts))
                    else:
                        in_isoform = filter(lambda fr: fr[0]>fr[1] and fr[0] in genome_coords and fr[1] in genome_coords, product(fwd_starts, rev_stops))

                    intervene_lens = map(lambda fr: isoform.getSequenceLength(fr[0],fr[1]), in_isoform)
                    can_pair = any(map(lambda j: j <= max_intervene_len, intervene_lens))
                    if (can_pair):
                        break   # All that we need to know is that at least one target isoform can be primed to produce a sequenceable amplicon

        return can_pair


    def getAsBED12Line(self, label_prefix, score=0, itemRgb="0,0,0"):
        genomic_positions = self.getGenomicPositions(as_set=True)
        genomic_positions = sorted(genomic_positions)
        start = genomic_positions[0]
        stop = genomic_positions[-1]
        
        block_sizes = []
        block_starts = []
        for k, g in groupby(enumerate(genomic_positions), lambda ix:ix[0]-ix[1]):
            group = list(map(itemgetter(1), g))
            block_starts.append( str(group[0] - start) )
            block_sizes.append( str(group[-1] - group[0] + 1) )
        
        # The "strand" field for the BED line depends on the strand of the isoforms and whether this is for fwd or rev primers
        strand = "+" if ((self.strand=="+" and self.rev_fwd=="fwd") or (self.strand=="-" and self.rev_fwd=="rev")) else "-"

        label = "%s%s" % (label_prefix, self.ID)
        bed12_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s" % \
          (self.chromosome, start-1, stop, label, score, strand, start, start, itemRgb, len(block_starts), ','.join(block_sizes), ','.join(block_starts))

        return bed12_line
