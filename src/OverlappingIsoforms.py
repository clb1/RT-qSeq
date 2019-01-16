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

class OverlappingIsoforms(object):
    def __init__(self, ID, all_overlapping_isoforms, target_isoform_groups, complete_TIG_IDs_w_nonlocal_members):
        self.ID = ID
        self.region_start = 1e15
        self.region_stop = -1
        self.isoforms = all_overlapping_isoforms
        self.isoforms_seqs = None
        self.target_isoform_groups = target_isoform_groups
        self.complete_TIG_IDs_w_nonlocal_members = complete_TIG_IDs_w_nonlocal_members
        self.target_isoforms = set(chain.from_iterable(target_isoform_groups))
        self.primerable_target_isoform_groups = set()
        self.chromosome = None
        self.strand = None
        self.genomic_position_data = {}
        self.equivalence_regions_rev = [] # For reverse primers
        self.equivalence_regions_fwd = [] # For forward primers
        self.fwd_rev_equivalence_regions_pairs = []
        self.primers_designs_sets = None
        
        strands = set()
        chromosomes = set()
        for isoform in self.target_isoforms:
            strands.add(isoform.getStrand())
            chromosomes.add(isoform.getChromosome())

        assert (len(strands)==1)
        assert (len(chromosomes)==1)

        self.strand = list(strands)[0]
        self.chromosome = list(chromosomes)[0]

        self.initGenomicPositionData()

        assert (len(self.target_isoform_groups) > 0), "No target isoform group(s) for %s" % self.ID


    def getChromosomeAndID(self):
        return (self.chromosome, self.ID)
    
    def getID(self):
        return self.ID

    
    def setID(self, ID):
        self.ID = ID


    def getStrand(self):
        return self.strand


    def getRegionStartStop(self):
        return (self.region_start, self.region_stop)


    def getSizeDescription(self):
        assert (not designParams.is_directional or all(map(lambda i: self.strand==i.getStrand(), self.isoforms))),\
          "Directionality constraint not applied to isoforms in olap set"
        return (len(self.target_isoform_groups), len(self.target_isoforms), len(self.isoforms))


    def getAllOverlappingIsoforms(self):
        return self.isoforms


    def getAllOverlappingIsoformsSeqs(self):
        assert (self.isoforms_seqs != None), "Isoforms' sequences have not been set"
        return self.isoforms_seqs


    def getAllOverlappingIsoformsIDs(self):
        return list(map(methodcaller("getCGDBName"), self.isoforms))


    def getTargetIsoformGroups(self, selector):
        if (selector=="all"):
            return self.target_isoform_groups
        elif (selector=="primerable"):
            return self.primerable_target_isoform_groups
        else:
            print("ERROR: unrecognized target isoform groups selector -> %s" % selector, file=sys.stderr, flush=True)
            sys.exit(1)
            

    def getNumTargetIsoformGroups(self, selector):
        if (selector=="all"):
            return len(self.target_isoform_groups)
        elif (selector=="primerable"):
            return len(self.primerable_target_isoform_groups)
        else:
            print("ERROR: unrecognized target isoform groups selector -> %s" % selector, file=sys.stderr, flush=True)
            sys.exit(1)


    def getGlobalTIGsOfLocalTIGs(self):
        return self.complete_TIG_IDs_w_nonlocal_members


    def getEquivPrimerRegionPairs(self):
        return self.fwd_rev_equivalence_regions_pairs
        

    def checkDirectionality(self):
        assert (not designParams.is_directional or all(map(lambda i: self.strand==i.getStrand(), self.isoforms))),\
          "Directionality constraint not applied to isoforms in olap set"
        

    def initGenomicPositionData(self):
        connected_genomic_positions = set()

        assert (not designParams.is_directional or all(map(lambda i: self.strand==i.getStrand(), self.isoforms))),\
          "Directionality constraint not applied to isoforms in olap set"

        genome_positions = set()
        for isoform in self.isoforms:
            start, stop = isoform.getStartStop()
            self.region_start = min(start, self.region_start)
            self.region_stop = max(stop, self.region_stop)
                
            for chrom, exon_start, exon_stop in isoform.getExonsAsTuples():
                exon_positions = range(exon_start,exon_stop+1)
                genome_positions.update(exon_positions)
                for pos in exon_positions:
                    if (pos not in self.genomic_position_data):
                        self.genomic_position_data[pos] = GenomicPositionData(pos)
                    self.genomic_position_data[pos].addIsoform(isoform)

        # Group the genomic positions that have the same set of isoforms
        #grouped_positions = defaultdict(list)
        #for pos in genome_positions:
        #    isoforms = self.genomic_position_data[pos].getIsoforms()
        #    isoforms_signature = tuple(sorted(isoforms, key=methodcaller("getCGDBName")))
        #    grouped_positions[isoforms_signature].append(pos)
                
        # For each distinct isoform set, get the group(s) contiguous of genomic positions.
        # These are the contiguous regions with constant isoforms association (ie same "signature" of isoforms).
        #self.contig_same_isoforms_region = []
        #for positions_with_same_isoforms in grouped_positions.values():
        #    assert (len(positions_with_same_isoforms) == len(set(positions_with_same_isoforms))), "Shouldn't have accumulated duplicate positions"
        #    positions_with_same_isoforms = sorted(positions_with_same_isoforms)
        #    for k, g in groupby(enumerate(positions_with_same_isoforms), lambda (i,x):i-x):
        #        self.contig_same_isoforms_region.append( map(itemgetter(1), g) )


    #def setPrimersDesignsSets(self, PDSs):
    #    assert (len(PDSs) > 0)
    #    self.primers_designs_sets = PDSs
    #    for experiment_number, pds in enumerate(PDSs):
    #        pds.setExperimentNumber(experiment_number)
    #        pds.annotateIsoformsWithAmplicons()
            

    def writeExperimentsDesigns(self, op_designs):
        for pds in self.primers_designs_sets:
            pds.writeDesign(op_designs)


    def setIsoformSequences(self, genome_ref):
        for isoform in self.isoforms:
            isoform.setSequence(genome_ref)


    def setIsoformsCoordCorrespondData(self, set_rna_to_genome=False, set_genome_to_rna=True):
        for isoform in self.isoforms:
            isoform.correspondmRNACoordsToGenomicCoords(set_rna_to_genome, set_genome_to_rna)


    def unsetIsoformsCoordCorrespondData(self, clear_rna_to_genome=False, clear_genome_to_rna=True):
        for isoform in self.isoforms:
            isoform.clearCorrespondmRNACoordsToGenomicCoords(clear_rna_to_genome, clear_genome_to_rna)


    def setIsoformsSeqs(self, transcriptome_ref):
        assert (self.isoforms_seqs == None), "Isoforms' sequences already set"
        self.isoforms_seqs = {}
        for isoform in self.isoforms:
            isoform_ID = isoform.getCGDBName()
            self.isoforms_seqs[isoform_ID] = str(transcriptome_ref[isoform_ID])


    def unsetIsoformsSeqs(self):
        self.isoforms_seqs = None


    def setPrimerGroupsPerPosition(self):
        to_remove = []

        for pos, gpd in self.genomic_position_data.items():
            gpd.setPrimerGroups()
            if (gpd.hasNoData()):
                to_remove.append(pos)

        # Drop positions that support neither fwd nor rev primers
        for pos in to_remove:
            del self.genomic_position_data[pos]


    def getCorrespondingNucSeq(self, genome_ref, genomic_positions, reverse_complement):
        genomic_positions = sorted(genomic_positions)

        nuc_seq = ''
        for k, g in groupby(enumerate(genomic_positions), lambda ix:ix[0]-ix[1]):
            group = list(map(itemgetter(1), g))
            nuc_seq += genome_ref[self.chromosome][int(group[0])-1:int(group[-1])]
        nuc_seq = str(nuc_seq)
        
        if (reverse_complement):
            nuc_seq = nuc_seq[::-1].translate(DNA_complement_table)
                
        return nuc_seq


    def getPositionsAsBED12Line(self, positions, label, score=0, itemRgb="0,0,0"):
        genomic_positions = sorted(positions)
        start = genomic_positions[0]
        stop = genomic_positions[-1]
        
        block_sizes = []
        block_starts = []
        for k, g in groupby(enumerate(genomic_positions), lambda ix:ix[0]-ix[1]):
            group = list(map(itemgetter(1), g))
            block_starts.append( str(group[0] - start) )
            block_sizes.append( str(group[-1] - group[0] + 1) )

        strand = "."

        bed12_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s" % \
          (self.chromosome, start-1, stop, label, score, strand, start, start, itemRgb, len(block_starts), ','.join(block_sizes), ','.join(block_starts))

        return bed12_line


    # 1) For each isoform signature, put all 5' primer positions in a graph and add edges between adjacent ones.
    # 2) Each CC of this graph is an EquivRegion, but an incomplete EquivRegion. The CC is just the primer 5p positions.
    #     For each 5p position, get the full primer positions. Use all, then, to form the EquivRegion.
    def buildEquivPrimerTargetsRegions(self, genome_ref):
        opposite_strand_isoforms = list(filter(lambda i: i.getStrand() != self.strand, self.isoforms))

        # Across all contiguous regions having the same isoform signature, 
        # compile genomic positions with same primer group, for all primer groups.
        for rev_fwd, label in [("rev", "R"), ("fwd", "F")]:
            unnamed_equivalence_regions = []
            known_adjacent = defaultdict(set)
            #primer_position_paths_per_5p = defaultdict(list)
            primer_5p_pos_graph_per_ig = defaultdict(nx.Graph)
            primer_path_per_isoform = {}
            for genome_5p_pos, gpd in self.genomic_position_data.items():
                for primer_seq, isoforms_signature, primer_position_paths in gpd.getPrimerGroups(rev_fwd):
                    primer_5p_pos_graph_per_ig[isoforms_signature].add_node(genome_5p_pos)
                    for primer_position_path in map(itemgetter(1), primer_position_paths):
                        known_adjacent[isoforms_signature].add(primer_position_path[0:2]) # Cache known primer 5' position and adjacent position
                    for isoform, primer_path in primer_position_paths:
                        primer_path_per_isoform[(genome_5p_pos,isoform)] = primer_path

            for isoforms_signature, G in primer_5p_pos_graph_per_ig.items():
                for pos_i, pos_j in known_adjacent[isoforms_signature]:
                    if (pos_i in G and pos_j in G and all(map(lambda x: x.areSequential(pos_i, pos_j), isoforms_signature))):
                        G.add_edge(pos_i, pos_j)

            # Form a primer EquivalenceRegion for each disjoint set of contiguous primer positions for each isoforms_signature.
            for isoforms_signature, genomic_positions_graph in primer_5p_pos_graph_per_ig.items():
                # Get connected components in Graph for each signature. Each CC is an equiv region
                for contiguous_primer_start_positions in nx.connected_component_subgraphs(genomic_positions_graph):
                    primer_positions_per_isoform = defaultdict(set)
                    for (genome_5p_pos, isoform) in product(contiguous_primer_start_positions.nodes(), isoforms_signature):
                        p = primer_path_per_isoform[(genome_5p_pos, isoform)]
                        primer_positions_per_isoform[isoform].update( p )

                    sorted_primer_positions_per_isoform = {}
                    isoform_subseqs_in_region = set()
                    for isoform, primer_positions in primer_positions_per_isoform.items():
                        reverse_complement = isoform.getStrand() != self.strand
                        sorted_primer_positions = sorted(primer_positions)

                        transcript_positions = isoform.convertToTranscriptCoords(sorted_primer_positions)
                        assert (all(map(lambda x: abs(x[1]-x[0])==1, zip(transcript_positions[0:-1], transcript_positions[1:])))), "Possible primer positions not all contiguous on isoform"

                        sorted_primer_positions_per_isoform[isoform] = tuple(sorted_primer_positions)
                        primer_src_seq = self.getCorrespondingNucSeq(genome_ref, sorted_primer_positions, reverse_complement)
                        isoform_subseqs_in_region.add(primer_src_seq)

                        # TODO DELETE?
                        #extrema_start, extrema_stop = sorted_primer_positions[0], sorted_primer_positions[-1]
                        #if (isoform.getStrand() == "+"):
                        #    extrema_seq = isoform.getSequence(None, extrema_start, extrema_stop)
                        #    if (self.strand == "-"):
                        #        extrema_seq = extrema_seq[::-1].translate(DNA_complement_table)
                        #else:
                        #    extrema_seq = isoform.getSequence(None, extrema_stop, extrema_start)
                        #    if (self.strand == "+"):
                        #        extrema_seq = extrema_seq[::-1].translate(DNA_complement_table)
                        #isoform_subseqs_in_region.add(extrema_seq)

                    assert (len(isoform_subseqs_in_region)==1), "Equiv region should only have one corresponding nucleotide sequence"
                        
                    region_sequence = list(isoform_subseqs_in_region)[0] 
                    if (self.strand == '-'):
                        region_sequence = region_sequence[::-1].translate(DNA_complement_table)

                    primer_positions_for_region = set(chain.from_iterable(sorted_primer_positions_per_isoform.values()))
                    opposite_overlapping_isoforms = set(filter(lambda i: i.numSharedPositions(primer_positions_for_region) >= designParams.min_primer_len,
                                                               opposite_strand_isoforms))

                    # equiv_region_position_path_graph, TODO DELETE
                    equiv_set = EquivPrimerTargetsRegion(self.ID, "NOT SET", self.chromosome, self.strand, sorted_primer_positions_per_isoform,
                                                         region_sequence, isoforms_signature, rev_fwd, opposite_overlapping_isoforms)

                    unnamed_equivalence_regions.append( (equiv_set, equiv_set.getExtremalStartStop()) )


            # Instantiate a canonical ordering and naming of the equiv_sets
            unnamed_equivalence_regions.sort(key=itemgetter(1))
            for counter, equiv_set in enumerate(map(itemgetter(0), unnamed_equivalence_regions), 1):
                regionID = "%s%d" % (label, counter)
                equiv_set.setID(regionID)
                if (rev_fwd == "rev"):
                    self.equivalence_regions_rev.append( equiv_set )
                else:
                    self.equivalence_regions_fwd.append( equiv_set )


    def pairPrimerRegions(self, genome_ref):
        possible_primer_lengths = list(range(designParams.min_primer_len, designParams.max_primer_len+1))

        for equiv_region_fwd in self.equivalence_regions_fwd:
            equiv_region_fwd.setPositionDependentPrimerData(possible_primer_lengths)

        for equiv_region_rev in self.equivalence_regions_rev:
            equiv_region_rev.setPositionDependentPrimerData(possible_primer_lengths)

        for equiv_region_fwd in self.equivalence_regions_fwd:
            for equiv_region_rev in self.equivalence_regions_rev:
                if (equiv_region_fwd.couldPotentiallyProductivelyPairWith(equiv_region_rev, self.target_isoforms)):
                    candidate_pair = EquivPrimerTargetsRegionPair(self.ID, equiv_region_fwd, equiv_region_rev)
                    candidate_pair.setIsoformGroupingsPrimerPositions(genome_ref, self.target_isoform_groups)
                    if (candidate_pair.couldAmplifyAnyTargetGroups(self.target_isoform_groups)):
                        self.fwd_rev_equivalence_regions_pairs.append(candidate_pair)


    def designDefaultPrimersForAllPairedPrimerRegions(self, oligo_thermo, tempdir, transcriptome_vsearch_udb, transcriptome_ref, antitargets_fasta_name):
        '''Delays the verification procedure until after the primer pairs have been found for all PPRs.'''
        max_promiscuous = designParams.max_primer_promiscuity
        num_success = 0

        olap_set_target_isoforms_IDs = set(map(lambda i: i.getCGDBName(), self.target_isoforms))

        # PrimerSingleton objects that pass on thermodynamic grounds.
        # Indexed like primer_caches["F#"|"R#"][(primer_genomic_5p, primer_len)] = PrimerSingleton or None if primer sequence fails on thermodynamic grounds
        #                                                                          or primes multiple times on target isoform
        # NOTE: For "_Abutt" case, primer_genomic_5p is instead a tuple of alt_primer_genomic_5p
        primer_caches = defaultdict(dict)

        local_transcriptome_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
        isoforms_seqs = self.getAllOverlappingIsoformsSeqs()
        for isoform_ID, isoform_seq in isoforms_seqs.items():
            local_transcriptome_fasta.write(">%s\n%s\n" % (isoform_ID, isoform_seq))
        local_transcriptome_fasta.flush()
        try:
            makeblastdb_output = check_output(['makeblastdb', '-dbtype', 'nucl', '-in', local_transcriptome_fasta.name], stderr=subprocess.STDOUT)
        except CalledProcessError as cpe:
            pdb.set_trace()

        # Collect candidate default primer pairs for all PPRs
        candidates_per_ppr = {}
        for ppr_num, ppr in enumerate(self.fwd_rev_equivalence_regions_pairs):
            fwd_erID, rev_erID = ppr.getFwdRevIDs()
            ppr_ID = ppr.getID() 
            descr = "disjoint" if ppr.isDisjoint() else "merged"
            print("\nPPR %s (%s), %d of %d" % (ppr_ID, descr, ppr_num+1, len(self.fwd_rev_equivalence_regions_pairs)), file=sys.stderr, flush=True)
        
            for igs_tuple, tigs_by_ID, ntigs_by_ID in ppr.getRankedIGsTuples():
               ppr_unique_target_isoform_IDs = olap_set_target_isoforms_IDs & set(map(methodcaller("getCGDBName"), chain.from_iterable(igs_tuple)))
               primer_pairs = ppr.designPairedPrimers(oligo_thermo, igs_tuple, ppr_unique_target_isoform_IDs, primer_caches[fwd_erID], primer_caches[rev_erID],
                                                      antitargets_fasta_name, local_transcriptome_fasta, tempdir)
               if (len(primer_pairs) > 0):
                   candidates_per_ppr[ppr] = (primer_pairs, igs_tuple, tigs_by_ID, ntigs_by_ID)
                   break

        local_transcriptome_fasta.close()
        os.remove("%s.nhr" % local_transcriptome_fasta.name)
        os.remove("%s.nin" % local_transcriptome_fasta.name)
        os.remove("%s.nsq" % local_transcriptome_fasta.name)

        # Evaluate candidate primers en masse by performing one transcriptome-wide search with all of them, and then evaluate on a per-PPR basis
        if (len(candidates_per_ppr) > 0):
            query_primers = {}
            for ppr, (primer_pairs, igs_tuple, tigs_by_ID, ntigs_by_ID) in candidates_per_ppr.items():
                ppr_ID = ppr.getID()
                for pp_num, (fwd_primer, rev_primer) in enumerate(primer_pairs):
                    fwd_id = "%s.%d_fwd" % (ppr_ID, pp_num)
                    rev_id = "%s.%d_rev" % (ppr_ID, pp_num)
                    query_primers[fwd_id] = fwd_primer.seq
                    query_primers[rev_id] = rev_primer.seq
                    
            print("INFO: running blastn and grouping results by primer pair", file=sys.stderr, flush=True)
            unverified_pp_grouped_by_ppr, hyb_events_per_primer = runThoroughBlastAndGroupResultsByPrimerPair(query_primers, "/raid1/projects/CGDB/models/indexes/BLAST/CGDBv2.0.fa", oligo_thermo, tempdir) # TODO: remove hardcoded path

            # Confirm that the candidate primer pairs amplified expected isoform-distinguishing amplicons
            print("INFO: Confirming amplified isoforms groups per primer pair", file=sys.stderr, flush=True)
            for ppr, (primer_pairs, igs_tuple, tigs_by_ID, ntigs_by_ID) in candidates_per_ppr.items():
                ppr_ID = ppr.getID()
                ig_tuples_by_ID = set(map(lambda ig: tuple(map(methodcaller("getCGDBName"), ig)), igs_tuple))
                verified_primer_pairs = filterPrimerPairsByGlobalSpecificityAndPromiscuity(ppr_ID, ig_tuples_by_ID, query_primers, primer_pairs, unverified_pp_grouped_by_ppr,
                                                                                           hyb_events_per_primer, self.complete_TIG_IDs_w_nonlocal_members, oligo_thermo,
                                                                                           max_promiscuous, transcriptome_ref)

                if (len(verified_primer_pairs)>0):
                    ppr.setDefaultPrimerPairs(igs_tuple, tigs_by_ID, ntigs_by_ID, verified_primer_pairs)
                    num_success += 1
                    print("- %s\tDefault primer pairs found (n = %d)" % (ppr_ID, len(verified_primer_pairs)), file=sys.stderr, flush=True)
                else:
                    print("- %s\tNo default primer pairs found" % ppr_ID, file=sys.stderr, flush=True)

        return (num_success, len(self.fwd_rev_equivalence_regions_pairs))


    def setPrimerableTargetIsoformGroups(self):
        results_descr_strings = []
        self.primerable_target_isoform_groups.clear()
        primerable_target_isoform_groups_str = []
        unprimerable_target_isoform_groups_str = []

        for ppr in filter(lambda ppr: ppr.hasDefaultPrimerPairs(), self.fwd_rev_equivalence_regions_pairs):
            igs = ppr.getSetOfPrimedIsoformGroups()
            self.primerable_target_isoform_groups.update( self.target_isoform_groups & igs )

        unprimerable_target_isoform_groups = self.target_isoform_groups - self.primerable_target_isoform_groups
        counter = 1
        if (len(unprimerable_target_isoform_groups) > 0):
            results_descr_strings.append("Cannot target %d of %d isoform groups for isoform overlap set %s:" % \
              (len(unprimerable_target_isoform_groups), len(self.target_isoform_groups), self.ID))
            for ig in unprimerable_target_isoform_groups:
                unprimerable_target_isoform_groups_str.append(" ".join(map(methodcaller("getCGDBName"), ig)))
                results_descr_strings.append("\t%d : %s" % (counter, unprimerable_target_isoform_groups_str[-1]))
                counter += 1
            results_descr_strings.append("\t------------------------------------------------------------------")
        else:
            results_descr_strings.append( "Can target all isoform groups:" )

        for ig in self.primerable_target_isoform_groups:
            primerable_target_isoform_groups_str.append(" ".join(map(methodcaller("getCGDBName"), ig)))
            results_descr_strings.append("\t%d : %s" % (counter, primerable_target_isoform_groups_str[-1]))
            counter += 1    

        return ("\n".join(results_descr_strings), len(unprimerable_target_isoform_groups), len(self.target_isoform_groups),
                " | ".join(primerable_target_isoform_groups_str), " | ".join(unprimerable_target_isoform_groups_str))


    def writeEquivalenceRegions(self, op, write_singletons, write_paired): # op_details, 
        if (write_singletons):
            for region in self.equivalence_regions_rev:
                line = region.getAsBED12Line(self.ID, 0, "255,150,255")
                op.write("%s\n" % line)

            for region in self.equivalence_regions_fwd:
                line = region.getAsBED12Line(self.ID, 0, "255,51,153")
                op.write("%s\n" % line)

        if (write_paired):
            for region_pair in self.fwd_rev_equivalence_regions_pairs:
                line = region_pair.getAsBED12Line(self.ID, 0, "0,255,0")
                op.write("%s\n" % line)

        
    def getEquivPrimerRegionPairs(self, only_with_default_primer_pair):
        if (only_with_default_primer_pair):
            return list(filter(lambda ppr: ppr.hasDefaultPrimerPairs(), self.fwd_rev_equivalence_regions_pairs))
        else:
            return self.fwd_rev_equivalence_regions_pairs


    def getEquivPrimerRegionPairsByName(self):
        return dict(map(lambda x: (x.getID(), x), self.fwd_rev_equivalence_regions_pairs))
        

    def addCompleteDesignToDatabase(self, conn, cursor):
        complete_designs_tuple = (self.chromosome, self.ID, len(self.primers_designs_sets), ' ', ' ')
        cursor.execute('INSERT INTO completedesigns VALUES (?,?,?,?,?)', complete_designs_tuple)

        primers_tuples = []
        for pds in self.primers_designs_sets:
            for ppr in pds.getPrimerPairRegions():
                primer_pairs = ppr.getPairedPrimers()
                primers_tuples.append( (self.chromosome, self.ID, pds.getExperimentNumber(),
                                        primer_pairs[0][3], primer_pairs[0][4],
                                        primer_pairs[0][10], primer_pairs[0][11],
                                        primer_pairs[0][15]) )

        cursor.executemany('INSERT INTO primers VALUES (?,?,?,?,?,?,?,?)', primers_tuples)

        conn.commit()
