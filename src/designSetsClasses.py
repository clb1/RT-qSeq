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
import SetOfPrimers
import RNAIsoform2
import GenomicPositionData
import OligoThermodynamics
import OverlappingIsoforms
import EquivPrimerTargetsRegion
import EquivPrimerTargetsRegionPair
import PPRClassCombos

import pdb

DNA_complement_table = str.maketrans("ACGTNacgtn","TGCANtgcan")

# Problems & Solution
# Problem2: Multiple matches in local transcriptome not being detected

# Problem3: Offtargeting due to decent F/R primer matching
# TODO: disallow primers to be in the same pool if they have strong and nearlby, but incomplete, matches to the same transcript. E.g. 15 nt 3' primer match good enough to be like real primer.

# Problem4: Spurious matching to RNA TE
# TODO: Add all RNA-TE to antitargets, possibly with different filtering criteria

# Problem5: primer that matches perfectly to target but not-quite-matches at 3' to other transcripts in same target locus creates products from not-quite-matching primer.
# TODO: disallow primer if it's other hybridizations in same locus are identical except for less than 3 of its last 4 3' nts


# frac_duplexed is fraction of primer that are calculated to be duplexed to its exact complementary template at the reaction temperature
# In addition to thermodynamic considerations, want to take into account primer 3' "goodness". A 3' is "better" when it has
# a balance of A+T vs C+G and when it has G or C as its 3' nucleotide
PrimerSingleton = namedtuple("PrimerSingleton", ["genome_5p_pos", "seq", "len", "Tm", "frac_duplexed", "genomic_positions", "thermo_penalty", "aux_3p_penalty"] )


class OrderedDiGraph(nx.DiGraph):
    adjlist_dict_factory = OrderedDict
    

class TimeoutException(Exception):
    '''Indicates that some logical block of code has run longer than it's alotted time.'''
    def __init__(self, message, solutions, *args):
        self.message = message
        self.solutions = solutions
        super(TimeoutException, self).__init__(message, solutions, *args)

    def hasSolutions(self):
        return len(self.solutions) > 0

    def getSolutions(self):
        return self.solutions


def selectDiversityOfPrimers(ranked_groups_of_annot_primers, max_needed, max_num_allowed_to_exceed_limit):
    indices_of_selected = [0]
    primers_positions_of_selected = [list(map(lambda p: set(p[1].genomic_positions), ranked_groups_of_annot_primers[0]))]
        
    num_allowed_to_exceed_limit = 0
    while (num_allowed_to_exceed_limit <= max_num_allowed_to_exceed_limit and len(indices_of_selected) < max_needed):
        for next_group_index, next_group in enumerate(ranked_groups_of_annot_primers[1:], 1):
            if (next_group_index in indices_of_selected):
                continue
            next_group_primers_positions = list(map(lambda p: set(p[1].genomic_positions), next_group))

            overlap_limit_exceeded = False
            for select_group_primers_positions in primers_positions_of_selected:
                truth_vals = list(map(lambda posAB: len(posAB[0] ^ posAB[1])<=6, product(next_group_primers_positions, select_group_primers_positions)))
                overlap_limit_exceeded = sum(truth_vals) > num_allowed_to_exceed_limit
                if (overlap_limit_exceeded):
                    break

            if (not overlap_limit_exceeded):
                indices_of_selected.append(next_group_index)
                primers_positions_of_selected.append( next_group_primers_positions )

                if (len(indices_of_selected) == max_needed):
                    break

        num_allowed_to_exceed_limit += 1
        
    return list(map(lambda i: ranked_groups_of_annot_primers[i], indices_of_selected))


def addsToDiversity(groups_of_annot_primers, new_annot_primers):
    adds_to_diversity = True

    primers_positions_of_new = list(map(lambda p: set(p[1].genomic_positions), new_annot_primers))
    for group_of_annot_primers in groups_of_annot_primers:
        group_primers_positions = list(map(lambda p: set(p[1].genomic_positions), group_of_annot_primers))
        truth_vals = map(lambda posAB: len(posAB[0] - posAB[1])>=3 and len(posAB[1] - posAB[0])>=3, product(group_primers_positions, primers_positions_of_new))
        adds_to_diversity = all(truth_vals)
        if (not adds_to_diversity):
            break

    return adds_to_diversity


def _expand(G, explored_nodes, explored_edges, explored_PPRs):
    """
    Expand existing solution by a process akin to BFS.

    Arguments:
    ----------
    G: networkx.Graph() instance
        full graph

    explored_nodes: set of ints
        nodes visited

    explored_edges: set of 2-tuples
        edges visited

    Returns:
    --------
    solutions: list, where each entry in turns contains two sets corresponding to explored_nodes and explored_edges
        all possible expansions of explored_nodes and explored_edges

    """
    frontier_nodes = list()
    frontier_edges = list()
    frontier_PPRs = set()
    for v in explored_nodes:
        for u in nx.neighbors(G,v):
            if (u not in explored_nodes and u[1] not in explored_PPRs):
                frontier_PPRs.add(u[1])
                frontier_nodes.append(u)
                frontier_edges.append([(u,v), (v,u)])

    explored_PPRs.update(frontier_PPRs)
    return zip([explored_nodes | frozenset([v]) for v in frontier_nodes], [explored_edges | frozenset(e) for e in frontier_edges])


def find_all_spanning_trees(G, root, G_num_diff_PPRs):
    """
    Find all spanning trees of a Graph.

    Arguments:
    ----------
    G: networkx.Graph() instance
        full graph

    Returns:
    ST: list of networkx.Graph() instances
        list of all spanning trees

    """

    # initialise solution
    explored_nodes = frozenset([root])
    explored_edges = frozenset([])
    explored_PPRs = set([root[1]])
    solutions = [(explored_nodes, explored_edges)]

    # we need to expand solutions number_of_different_PPRs_in_G-1 times
    for ii in range(G_num_diff_PPRs-1):
        # get all new solutions
        new_solutions = [_expand(G, nodes, edges, explored_PPRs) for (nodes, edges) in solutions]
        # flatten nested structure and get unique expansions
        new_solutions = set([item for sublist in new_solutions for item in sublist])
        if (len(new_solutions) > 0):
            solutions = new_solutions
        else:
            break
        
    return [nx.from_edgelist(edges) for (nodes, edges) in solutions]


def findAllSpanningTrees(G, root_node, expected_tree_size):
    all_spanning_trees = []
    
    categorized_descendants = defaultdict(list) # Categorized by PPR
    root_node_ppr = root_node[1]
    categorized_descendants[root_node_ppr] = [root_node]
    
    # TODO: currently brute force culling of incorrect solutions.
    # There is certainly a smarter and more efficient way.
    for node in nx.descendants(G, root_node):
        if (node[1] != root_node_ppr):
            categorized_descendants[node[1]].append(node)

    for spanning_tree in product(*categorized_descendants.values()):
        subG = G.subgraph(spanning_tree)
        if (nx.is_connected(subG) and len(subG)==expected_tree_size):
            all_spanning_trees.append(spanning_tree)
            
    return all_spanning_trees


def findAllSpanningTrees2(G, curr_tree, explored_pprs): # , frontier_pprs
    all_spanning_trees = []
    
    categorized_descendants = defaultdict(list) # Categorized by PPR

    #for node in [n for n in curr_tree if n[1] in frontier_pprs]:
    for node in curr_tree:
        for neigh in [m for m in G.neighbors(node) if m[1] not in explored_pprs]:
            categorized_descendants[neigh[1]].append(neigh)

    #new_frontier_pprs = set(categorized_descendants.keys())
    new_explored_pprs = explored_pprs | set(categorized_descendants.keys()) # frontier_pprs

    if (len(categorized_descendants) != 0):
        for next_tree_part in product(*categorized_descendants.values()):
            append_trees = findAllSpanningTrees2(G, next_tree_part, new_explored_pprs) # , new_frontier_pprs
            all_spanning_trees.extend( [next_tree_part+append_tree for append_tree in append_trees] )
    else:
        all_spanning_trees.append( () )

    return all_spanning_trees


def hasTooManyRepeats(nuc_seq):
    '''Fail primers with a too-high dinucleotide frequency or a trinucleotide frequency. Indicates a primer 
    in a repeat region that would have very high number of non-specific matches across the transcriptome.
    '''
    nuc_counts = Counter(nuc_seq)
    dinuc_counts = Counter([nuc_seq[i:i+2] for i in range(0,len(nuc_seq)-1,1)])
    trinuc_counts = Counter([nuc_seq[i:i+3] for i in range(0,len(nuc_seq)-2,1)])
    max_nuc_count = int(len(nuc_seq)/2)
    max_dinuc_count = int((len(nuc_seq)-1)/3)
    max_trinuc_count = int((len(nuc_seq)-2)/4)
    too_many_repeats = nuc_counts.most_common(1)[0][1]>max_nuc_count or \
                       dinuc_counts.most_common(1)[0][1]>max_dinuc_count or trinuc_counts.most_common(1)[0][1]>max_trinuc_count
    return too_many_repeats


def scanIndivPrimersAgainstLocalTranscriptome(oligo_thermo, expected_isoforms_per_primer_seq, all_equiv_region_primers,
                                              antitargets_fasta_name, local_transcriptome_fasta, ppr_unique_target_isoform_IDs, tempdir):
    '''Evaluate each primer from each set of equiv region candidates individually, without regard to a priming partner. Filter out any primers that have 
    unexpected/unwanted hybridization potential. Only return results for which frac_duplexed of primer to target > 0.01. primer and target sequences are 
    returned in 5'->3' orientation.'''
    bad_primer_seqs = set()
    filtered_equiv_region_primers = None
    primer_num_nontarget_hybs = None

    # First, filter out primers that align to antitarget mRNAs
    primers = {}
    for equiv_region_ID, equiv_region_primers in all_equiv_region_primers.items():
        for primer_num, singleton_primer in enumerate(equiv_region_primers):
            primer_id = "%s.%d" % (equiv_region_ID, primer_num)
            primers[primer_id] = singleton_primer.seq

    ids_of_antitarget_matching_primers = scanPrimersAgainstAntitargets(primers, antitargets_fasta_name, tempdir, 'B')
    print("Antitarget filtering removed %d primers" % len(ids_of_antitarget_matching_primers), file=sys.stderr)

    # Second, align to local transcriptome those primers not matching to antitargets
    primers_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    num_primers_not_bad = 0
    for primer_id, primer_seq in primers.items():
        if (primer_id not in ids_of_antitarget_matching_primers):
            num_primers_not_bad += 1
            primers_fasta.write(">%s\n%s\n" % (primer_id,primer_seq))
        else:
            bad_primer_seqs.add(primer_seq)
    primers_fasta.flush()

    if (num_primers_not_bad > 0):
        ip_blastn = tempfile.NamedTemporaryFile(mode='w+t', suffix=".blastn", dir=tempdir, delete=True)
        blastn_cmd = "blastn -task blastn-short -query %s -db %s -out %s -max_target_seqs 500000 -num_threads 3 -dust no -perc_identity 80 -evalue 1000 -qcov_hsp_perc 33 \
                      -outfmt '6 qseqid sseqid qlen length sstrand qstart qend sstart send evalue pident nident mismatch gaps qseq sseq'" % \
                          (primers_fasta.name, local_transcriptome_fasta.name, ip_blastn.name)

        try:
            blastn_output = check_output(blastn_cmd, stderr=subprocess.STDOUT, shell=True)
        except CalledProcessError as cpe:
            pdb.set_trace()
        
        ip_blastn.flush()
        ip_blastn.seek(0)

        # Compile local transcriptome matches for each primer
        primed_target_isoform_counts = defaultdict(Counter)
        primed_nontarget_isoform_counts = defaultdict(Counter)
        for line in ip_blastn:
            #query, qstrand, qlo, qhi, qrow, target, tstrand, tilo, tihi, trow, mism, opens = line.strip().split('\t')
            qseqid, sseqid, qlen, length, sstrand, qstart, qend, sstart, send, evalue, pident, nident, mismatch, gaps, qseq, sseq = line.strip().split('\t')
            equiv_region_ID, primer_num = qseqid.rsplit('.',1)
            primer_fwdrev = "fwd" if (equiv_region_ID[0] == 'F') else "rev"
            match_fwdrev = "fwd" if (sstrand=='plus') else "rev"

            if (gaps == '0'):
                query_seq, template_seq = qseq, sseq
            else:
                query_seq = qseq.replace('-','')
                template_seq = sseq.replace('-','')

            template_seq = template_seq[::-1].translate(DNA_complement_table)
            any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
            end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
            min_dG = min(any_hyb.dg, end_hyb.dg)
            frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)

            if (frac_duplexed >= 0.01):
                if ((int(qlen) < int(length) and int(qlen) - int(qend) <= 3) or (qlen == length and qseq != sseq)):
                    bad_primer_seqs.add(query_seq) # Primes inexactly but can possibly do offtarget priming
                elif (sseqid in ppr_unique_target_isoform_IDs):
                    primed_target_isoform_counts[query_seq][sseqid] += 1
                else:
                    primed_nontarget_isoform_counts[query_seq][sseqid] += 1

        ip_blastn.close()
        primers_fasta.close()

        # Identify "bad" primers that prime more than once on *any* dual-primed PPR target isoform OR that prime on an unintended *dual-primed PPR target isoform*
        # This logic should catch primers that also act oppositely (eg forward primer acting as a reverse primer)
        for primer_seq, primer_target_counts in primed_target_isoform_counts.items():
            actual_targets_primed = set(primer_target_counts.keys())
            if (primer_target_counts.most_common(1)[0][1] > 1 or len(actual_targets_primed - expected_isoforms_per_primer_seq[primer_seq])>0):
                bad_primer_seqs.add(primer_seq)
        print("Local transcriptome matching removed %d primers" % (len(bad_primer_seqs)-len(ids_of_antitarget_matching_primers),), file=sys.stderr)
                    
        # Filter out "bad" singleton primers
        filtered_equiv_region_primers = {}
        primer_num_nontarget_hybs = {}
        for equiv_region_ID in all_equiv_region_primers.keys():
            filtered_equiv_region_primers[equiv_region_ID] = set(filter(lambda sp: sp.seq not in bad_primer_seqs, all_equiv_region_primers[equiv_region_ID]))
            for primer in filtered_equiv_region_primers[equiv_region_ID]:
                primer_num_nontarget_hybs[primer] = sum(map(itemgetter(1), primed_nontarget_isoform_counts[primer.seq].most_common()))

        if (any(map(lambda v: len(v)==0, filtered_equiv_region_primers.values()))):
            filtered_equiv_region_primers = None
            
    return filtered_equiv_region_primers, primer_num_nontarget_hybs


def scanPrimersAgainstAntitargets(primers, antitargets_fasta_name, tempdir, use_case='A'):
    '''Identifies primers that match to any sequence in antitargets_fasta, in either a forward or reverse sense.
    Matching is based not on thermodynamics but on strict sequence similarity so that even low-similarity matches 
    that wouldn't thermodynamically duplex are identified. For each fwd/rev primer that matches to antitarget
    sequences, that primer's associated rev/fwd primer is returned as well.

    Assumption: Primers are named <something>_fwd or <something>_rev.
    '''

    primers_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    for primer_ID, primer_seq in primers.items():
        primers_fasta.write(">%s\n%s\n" % (primer_ID, primer_seq))
    primers_fasta.flush()
            
    ip_blastn = tempfile.NamedTemporaryFile(mode='w+t', suffix=".blastn", dir=tempdir, delete=True)
    blastn_cmd = "blastn -task blastn-short -query %s -db %s -max_target_seqs 5000 -num_threads 2 -dust no -perc_identity 80 -evalue 1000 -qcov_hsp_perc 30 \
                 -word_size 6 -outfmt '6 qseqid qlen length qend pident nident mismatch gaps qseq sseq' -out %s" % (primers_fasta.name, antitargets_fasta_name, ip_blastn.name)

    try:
        usearch_output = check_output(blastn_cmd, stderr=subprocess.STDOUT, shell=True)
    except CalledProcessError as cpe:
        pdb.set_trace()

    ip_blastn.flush()
    ip_blastn.seek(0)

    ids_of_matching_primers = set()
    if (use_case == 'A'):
        for line in ip_blastn:
            qseqid, qlen, length, qend, pident, nident, mismatch, gaps, qseq, sseq = line.strip().split('\t')
            mate_revfwd = "rev" if (qseqid[-3:] == "fwd") else "fwd"
            mate_qseqid = "%s_%s" % (qseqid[0:-4], mate_revfwd)

            dist_from_query_3p = int(qlen) - int(qend)
            nident = int(nident)
            if ((dist_from_query_3p <= 1 and nident >= 8) or (dist_from_query_3p == 2 and nident >= 10) or (dist_from_query_3p == 3 and nident > 12)):
                ids_of_matching_primers.update([qseqid,mate_qseqid])
            
    elif (use_case == 'B'):
        for line in ip_blastn:
            qseqid, qlen, length, qend, pident, nident, mismatch, gaps, qseq, sseq = line.strip().split('\t')
            dist_from_query_3p = int(qlen) - int(qend)
            nident = int(nident)
            if ((dist_from_query_3p <= 1 and nident >= 8) or (dist_from_query_3p == 2 and nident >= 10) or (dist_from_query_3p == 3 and nident > 12)):
                ids_of_matching_primers.add(qseqid)

    else:
        print("ERROR: unknown use case", file=sys.stderr, flush=True)
        sys.exit(1)

    ip_blastn.close()
    primers_fasta.close()

    return ids_of_matching_primers


def scanIndivPrimersAgainstAdapterTags(all_equiv_region_primers, tempdir):
    '''For extra level of stringency against primers forming stable hairpins with their adapter tags and against primers forming stable heteroduplexes 
    and/or extension products with their adapter tags, which are in very high concentration in a fwd/rev stranding reactions during RT-qSeq'''

    primers = {}
    primers_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    for equiv_region_ID, equiv_region_primers in all_equiv_region_primers.items():
        for primer_num, singleton_primer in enumerate(equiv_region_primers):
            primer_ID = "%s.%d" % (equiv_region_ID, primer_num)
            primers_fasta.write(">%s\n%s\n" % (primer_ID, singleton_primer.seq))
            primers[primer_ID] = singleton_primer.seq.upper()
    primers_fasta.flush()

    adapter_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    adapter_fasta.write(">Fwd_adapter\n%s\n>Rev_adapter\n%s\n" % \
                        (designParams.fwd_primer_adapter_tag[::-1].translate(DNA_complement_table), designParams.rev_primer_adapter_tag[::-1].translate(DNA_complement_table)))
    adapter_fasta.flush()

    ip_usearch = tempfile.NamedTemporaryFile(mode='w+t', suffix=".usearch", dir=tempdir, delete=True)
    usearch_cmd = ['usearch', '-search_local', adapter_fasta.name, '-db', primers_fasta.name, '-strand', 'plus', '-hspw', '4', '-evalue', '10000', '-mincols', '4',
                   '-userout', ip_usearch.name, '-userfields', 'query+qlo+qhi+ql+qstrand+target+tlo+thi+tl+tstrand+alnlen+aln+mism+ids+opens+diffs+qrow+trow']

    try:
        usearch_output = check_output(usearch_cmd, stderr=subprocess.STDOUT)
    except CalledProcessError as cpe:
        pdb.set_trace()

    ip_usearch.flush()
    ip_usearch.seek(0)

    bad_primer_seqs = set()
    for line in ip_usearch:
        query, qlo, qhi, ql, qstrand, target, tlo, thi, tl, tstrand, alnlen, aln, mism, ids, opens, diffs, qrow, trow = \
            map(lambda z: z[0](z[1]), zip((str,int,int,int,str,str,int,int,int,str,int,str,int,int,int,int,str,str),line.strip().split('\t')))

        # Alignment positions in the non-reverse complemented adapter
        #adapter_start = ql - qhi + 1
        #adapter_stop = ql - qlo + 1

        num_GC_matches = sum(list(map(lambda t: t==('G','G') or t==('C','C'), zip(qrow.upper(),trow.upper()))))

        if (query[0] == target[0]):
            # For forward/reverse primer against its forward/reverse adapter tag
            # Logic for badness, attached adapter tag (eg foward adapter tag that is attached to each forward primer)
            if ((thi>=tl-1 and alnlen >= 6 and opens == 0) or num_GC_matches >= 4 or ids >= 7) :
                bad_primer_seqs.add( primers[target] ) # The query name is the primer sequence

        elif ((thi>=tl-1 and ids >= 7 and opens == 0) or num_GC_matches >= 8 or ids >= 10):
            # For forward/reverse primer against reverse/forward adapter tag
            bad_primer_seqs.add( primers[target] ) # The query name is the primer sequence
            
    primers_fasta.close()
    adapter_fasta.close()
    ip_usearch.close()

    # Filter out "bad" singleton primers
    filtered_equiv_region_primers = {}
    for equiv_region_ID in all_equiv_region_primers.keys():
        filtered_equiv_region_primers[equiv_region_ID] = set(filter(lambda sp: sp.seq not in bad_primer_seqs, all_equiv_region_primers[equiv_region_ID]))

    if (any(map(lambda v: len(v)==0, filtered_equiv_region_primers.values()))):
        filtered_equiv_region_primers = None
            
    return filtered_equiv_region_primers


def detectPotentialPrimerInteractionOnGenome():
    blastn_cmd = ['blastn', '-query', rev_primers.fa, '-task', 'blastn-short', '-db', '/raid1/references_and_indexes/hg38/BLAST', '-max_target_seqs', '500000',
                  '-num_threads', '15', '-dust', 'no', '-perc_identity', '70', '-evalue', '100000', '-qcov_hsp_perc', '85', 
                  '-outfmt', "'6 qseqid sseqid qlen length sstrand qstart qend sstart send evalue pident nident mismatch gaps qseq sseq'", '-out', rev.blastn]

    ip = open(blastn_file, 'r')

    all_lines_by_chr = defaultdict(list)
    for line in ip:
        fields = line.strip().split("\t")
        cast_fields = list(map(lambda x: x[0](x[1]), zip((str,str,int,int,str,int,int,int,int,float,float,int,int,int,str,str),fields)))
        all_lines_by_chr[fields[1]].append(cast_fields)
    ip.close()

    all_report_lines = []
    for chromsome, chr_lines in all_lines_by_chr.items():
        chr_lines.sort(key=itemgetter(7,8))

        #report_lines = reportCloseby(chr_lines)
        report_lines = []
        for curr_line_index in range(len(chr_lines)-1):
            for next_line_index in range(curr_line_index+1,len(chr_lines)):
                if (chr_lines[curr_line_index][4] == chr_lines[next_line_index][4]):
                    continue

                if (abs(chr_lines[next_line_index][7] - chr_lines[curr_line_index][7]) > 100): # TODO: change dist
                    break
            
                if (chr_lines[curr_line_index][4] == 'plus'):
                    plus_match_pos = chr_lines[curr_line_index][7]
                    minus_match_pos = chr_lines[next_line_index][7]
                else:
                    plus_match_pos = chr_lines[next_line_index][7]
                    minus_match_pos = chr_lines[curr_line_index][7]

                if (minus_match_pos < plus_match_pos):
                    continue
            
                dist = minus_match_pos - plus_match_pos + 1
                if (dist < 100): # TODO: change dist
                    curr_line = "\t".join(map(str,chr_lines[curr_line_index]))
                    next_line = "\t".join(map(str,chr_lines[next_line_index]))
                    report_lines.append( (dist, curr_line, next_line) )
                else:
                    try:
                        assert (dist >= 0)
                    except AssertionError:
                        pdb.set_trace()
                    break

        all_report_lines.extend(report_lines)

    all_report_lines.sort(key=itemgetter(0))
    op = open(output_file, 'w')
    for (dist, curr_line, next_line) in all_report_lines:
        op.write("dist = %d\n%s\n%s\n\n" % (dist, curr_line, next_line))
    op.close()


def verifyNpartSolutionsGlobally(Npart_solutions, transcriptome_vsearch_udb, transcriptome_ref, oligo_thermo, complete_TIG_IDs_w_nonlocal_members, tempdir):
    # Accumulate all primers, align them against transcriptome, and group by query sequence
    query_primer_seqs = set()
    for ppr_combo, primer_sets in Npart_solutions:
        for primer in chain.from_iterable(map(methodcaller("getAllPrimers"), primer_sets)):
            query_primer_seqs.add(primer.seq)

    query_primers = {}
    for seq in query_primer_seqs:
        query_primers[seq] = seq
    
    all_vsearch_results, primers_w_no_perfect_match = runThoroughVsearch(query_primers, transcriptome_vsearch_udb, tempdir)

    grouped_vsearch_results = defaultdict(list)
    for query, qstrand, qlo, qhi, qrow, target, tstrand, tilo, tihi, trow, mism, opens in all_vsearch_results:
        assert (tstrand == '+'), "Template match is negative strand. How?"
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
            if (qstrand == '+'): # and '-' not in qrow[-6:] and '-' not in trow[-6:] and sum(map(lambda x:x[0]!=x[1], zip(qrow[-7:],trow[-7:]))) <= 1):
                template_seq = template_seq[::-1].translate(DNA_complement_table)
                any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                min_dG = min(any_hyb.dg, end_hyb.dg)
                frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)
                    
            else: # elif (qstrand == '-' and '-' not in qrow[0:6] and '-' not in trow[0:6] and sum(map(lambda x:x[0]!=x[1], zip(qrow[0:7],trow[0:7]))) <= 1):
                query_seq = query_seq[::-1].translate(DNA_complement_table)
                any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                min_dG = min(any_hyb.dg, end_hyb.dg)
                frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)

            if (frac_duplexed >= 0.01):
                data_tup = (match_fwdrev, query_seq, target, int(tilo), int(tihi), template_seq, opens != "0", frac_duplexed)
                grouped_vsearch_results[query].append( data_tup )


    # Try to verify each primer set for each solution part
    verified_Npart_solutions = []
    expected_num_unwanted_from_primers = 0
    for ppr_combo, primer_sets in Npart_solutions:
        ppr_combo_descriptor = "+".join(list(map(methodcaller("getID"), ppr_combo)))

        try:
            tigs_by_ID = list(chain.from_iterable(map(lambda x: x.local_target_isoform_groups_dual_primed_by_ID, ppr_combo)))
        except TypeError as te:
            pdb.set_trace()
            print(te, file=sys.stderr)
        assert (None not in tigs_by_ID), "None in tigs_by_ID"

        try:
            ntigs_by_ID = list(chain.from_iterable(map(lambda x: x.local_nontarget_isoform_groups_dual_primed_by_ID, ppr_combo)))
        except TypeError as te:
            print(te, file=sys.stderr)
            assert (None not in ntigs_by_ID), "None in ntigs_by_ID"
            
        expected_igs_tuples_by_ID = set(ntigs_by_ID + tigs_by_ID)
        target_isoform_IDs = set(list(chain.from_iterable(tigs_by_ID)))

        verified_sets_of_primers = []
        num_new_isoforms = 0
        for primer_set in primer_sets:
            all_primer_matches = []
            for erID, primer in primer_set.getAllAnnotPrimers():
                primer_fwdrev = "fwd" if (erID[0]=='F') else "rev"
                all_primer_matches.extend( list(map(lambda t: (primer_fwdrev,) + t, grouped_vsearch_results[primer.seq])) )
            ps_was_verified, new_igs_tuples_by_ID = \
                 verifyIsoformGroupsAreDistinctBySequencing2(all_primer_matches, oligo_thermo, expected_igs_tuples_by_ID, transcriptome_ref,
                                                             complete_TIG_IDs_w_nonlocal_members, target_isoform_IDs)
            if (ps_was_verified):
                primer_set.setIsoformGroupsFromCategorized(tigs_by_ID, ntigs_by_ID, new_igs_tuples_by_ID)
                verified_sets_of_primers.append(primer_set)
                num_new_isoforms += len(list(chain.from_iterable(new_igs_tuples_by_ID)))
        if (len(verified_sets_of_primers) > 0):
            verified_Npart_solutions.append( (ppr_combo, verified_sets_of_primers) )
            expected_num_unwanted_from_primers += float(num_new_isoforms)/float(len(verified_sets_of_primers))
        else:
            verified_Npart_solutions = None
            break

    return (verified_Npart_solutions, round(expected_num_unwanted_from_primers,1))


def runThoroughBlast(query_primers, transcriptome_blastn_db, tempdir):
    '''Accumulate perfect and imperfect matches of primers to transcriptome.'''
    all_blastn_results = set()

    op_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    for primer_id, primer_seq in query_primers.items():
        op_fasta.write(">%s\n%s\n" % (primer_id, primer_seq))
    op_fasta.flush()

    ip_blastn = tempfile.NamedTemporaryFile(mode='w+t', suffix=".blastn", dir=tempdir, delete=True) 
    blastn_cmd = "blastn -query %s -task blastn-short -db %s -out %s -max_target_seqs 50000 -num_threads 9 -dust no \
                  -outfmt '6 qseqid sseqid qlen length sstrand qstart qend sstart send evalue pident nident mismatch gaps qseq sseq' \
                  -word_size 7 -evalue 1000 -gapopen 2 -gapextend 2 -reward 1 -penalty -1 -perc_identity 70" % (op_fasta.name, transcriptome_blastn_db, ip_blastn.name)

    try:
        begin_time = datetime.datetime.now()
        begin_time_str = begin_time.strftime("%Y-%m-%d %H:%M:%S")
        print("\n\t\tINFO: starting blastn (%s)..." % begin_time_str, file=sys.stderr, end='', flush=True)
        blastn_output = check_output(blastn_cmd, stderr=subprocess.STDOUT, shell=True)
        end_time = datetime.datetime.now()
        end_time_str = end_time.strftime("%Y-%m-%d %H:%M:%S")
        print("done (%s)" % end_time_str, file=sys.stderr, flush=True)
        time_delta = end_time - begin_time
        seconds = time_delta.total_seconds() 
        print("\t\tINFO: runtime = %d hr, %d min, %d sec" % (int(seconds / (3600)), int(seconds / 60), seconds % 60), file=sys.stderr, flush=True)
    except CalledProcessError as cpe:
        pdb.set_trace()

    ip_blastn.flush()
    ip_blastn.seek(0)
    for line in ip_blastn:
        fields = line.strip().split('\t')
        all_blastn_results.add(tuple(fields))

    return all_blastn_results


def groupBlastnResultsByPrimerPair(all_blastn_results, oligo_thermo):
    '''Query sequence IDs must be in a specific format. See below. '''
    unverified_pp_grouped_by_ppr = defaultdict(list)
    hyb_events_per_primer = defaultdict(list)  # For identifying promiscuous binders
    frac_duplexed_cache = {}

    for (qseqid, sseqid, qlen, length, sstrand, qstart, qend, sstart, send, evalue, pident, nident, mismatch, gaps, qseq, sseq) in all_blastn_results:
        match_fwdrev = "fwd" if (sstrand=='plus') else "rev"
        query_len = int(qlen)
        dist_from_query_3p = query_len - int(qend)

        if (gaps == '0'):
            query_seq, template_seq = qseq, sseq
        else:
            query_seq = qseq.replace('-','')
            template_seq = sseq.replace('-','')

        src, primer_fwdrev = qseqid.rsplit('_',1)
        ppr_ID, pp_num = src.rsplit('.',1)
        pp_num = int(pp_num)
        frac_duplexed = 0

        cache_key = (sstrand, qseq, sseq)
        if (cache_key in frac_duplexed_cache):
            frac_duplexed, query_seq, template_seq = frac_duplexed_cache[cache_key]
            if (frac_duplexed >= 0.01):
                if (dist_from_query_3p <= 2): # Primer could potentially extend
                    data_tup = (qseqid, primer_fwdrev, match_fwdrev, query_seq, sseqid, int(sstart), int(send), template_seq, gaps != "0", dist_from_query_3p, frac_duplexed)
                    unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)].append(data_tup)
                hyb_events_per_primer[(ppr_ID, pp_num)].append((primer_fwdrev, sseqid, sstart, send, frac_duplexed))

        else:
            template_seq = template_seq[::-1].translate(DNA_complement_table)
            hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
            frac_duplexed = oligo_thermo.calcFracDuplexed(hyb.dg)
                    
            if (frac_duplexed >= 0.01):
                if (dist_from_query_3p <= 2): # Primer could potentially extend
                    data_tup = (qseqid, primer_fwdrev, match_fwdrev, query_seq, sseqid, int(sstart), int(send), template_seq, gaps != "0", dist_from_query_3p, frac_duplexed)
                    unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)].append(data_tup)
                hyb_events_per_primer[(ppr_ID, pp_num)].append((primer_fwdrev, sseqid, sstart, send, frac_duplexed))

            frac_duplexed_cache[cache_key] = (frac_duplexed, query_seq, template_seq)

    return unverified_pp_grouped_by_ppr, hyb_events_per_primer


def runThoroughBlastAndGroupResultsByPrimerPair(query_primers, transcriptome_blastn_db, oligo_thermo, tempdir):
    '''Accumulate perfect and imperfect matches of primers to transcriptome.'''
    unverified_pp_grouped_by_ppr = defaultdict(list)
    hyb_events_per_primer = defaultdict(list)  # For identifying promiscuous binders
    frac_duplexed_cache = {}

    op_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    for primer_id, primer_seq in query_primers.items():
        op_fasta.write(">%s\n%s\n" % (primer_id, primer_seq))
    op_fasta.flush()

    blastn_cmd = "blastn -query %s -task blastn-short -db %s -out - -max_target_seqs 50000 -num_threads 9 -dust no \
                  -outfmt '6 qseqid sseqid qlen length sstrand qstart qend sstart send evalue pident nident mismatch gaps qseq sseq' \
                  -word_size 7 -evalue 1000 -gapopen 2 -gapextend 2 -reward 1 -penalty -1 -perc_identity 70" % (op_fasta.name, transcriptome_blastn_db)

    begin_time = datetime.datetime.now()
    begin_time_str = begin_time.strftime("%Y-%m-%d %H:%M:%S")
    print("\n\t\tINFO: starting blastn and grouping (%s)..." % begin_time_str, file=sys.stderr, end='', flush=True)

    blastn_process = Popen(blastn_cmd, shell=True, stdout=subprocess.PIPE)

    line = blastn_process.stdout.readline().decode()
    while (not (line == '' and blastn_process.poll() is not None)):
        if (line != ''):
            (qseqid, sseqid, qlen, length, sstrand, qstart, qend, sstart, send, evalue, pident, nident, mismatch, gaps, qseq, sseq) = line.strip().split('\t')
            match_fwdrev = "fwd" if (sstrand=='plus') else "rev"
            query_len = int(qlen)
            dist_from_query_3p = query_len - int(qend)

            if (gaps == '0'):
                query_seq, template_seq = qseq, sseq
            else:
                query_seq = qseq.replace('-','')
                template_seq = sseq.replace('-','')

            src, primer_fwdrev = qseqid.rsplit('_',1)
            ppr_ID, pp_num = src.rsplit('.',1)
            pp_num = int(pp_num)
            frac_duplexed = 0

            cache_key = (sstrand, qseq, sseq)
            if (cache_key in frac_duplexed_cache):
                frac_duplexed, query_seq, template_seq = frac_duplexed_cache[cache_key]
                if (frac_duplexed >= 0.01):
                    if (dist_from_query_3p <= 2): # Primer could potentially extend
                        data_tup = (qseqid, primer_fwdrev, match_fwdrev, query_seq, sseqid, int(sstart), int(send), template_seq, gaps != "0", dist_from_query_3p, frac_duplexed)
                        unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)].append(data_tup)
                    hyb_events_per_primer[(ppr_ID, pp_num)].append((primer_fwdrev, sseqid, sstart, send, frac_duplexed))

            else:
                template_seq = template_seq[::-1].translate(DNA_complement_table)
                hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                frac_duplexed = oligo_thermo.calcFracDuplexed(hyb.dg)
                    
                if (frac_duplexed >= 0.01):
                    if (dist_from_query_3p <= 2): # Primer could potentially extend
                        data_tup = (qseqid, primer_fwdrev, match_fwdrev, query_seq, sseqid, int(sstart), int(send), template_seq, gaps != "0", dist_from_query_3p, frac_duplexed)
                        unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)].append(data_tup)
                    hyb_events_per_primer[(ppr_ID, pp_num)].append((primer_fwdrev, sseqid, sstart, send, frac_duplexed))

                frac_duplexed_cache[cache_key] = (frac_duplexed, query_seq, template_seq)
        line = blastn_process.stdout.readline().decode()

    rc = blastn_process.poll()
    assert (rc == 0), "blastn returned with exit code %d" % rc
    
    end_time = datetime.datetime.now()
    end_time_str = end_time.strftime("%Y-%m-%d %H:%M:%S")
    print("done (%s)" % end_time_str, file=sys.stderr, flush=True)
    time_delta = end_time - begin_time
    seconds = time_delta.total_seconds() 
    print("\t\tINFO: runtime = %d hr, %d min, %d sec" % (int(seconds / (3600)), int(seconds / 60), seconds % 60), file=sys.stderr, flush=True)

    return unverified_pp_grouped_by_ppr, hyb_events_per_primer


def runThoroughVsearch(query_primers, transcriptome_vsearch_udb, tempdir, max_maxrejects = 32768, num_vsearch_threads = 9):
    '''Accumulate perfect and imperfect matches of primers to transcriptome.'''
    primers_w_perfect_match = set()
    all_vsearch_results = set()
    maxrejects = 4096
    perc_w_perfect_match = 0
    while (perc_w_perfect_match < 100.0 and maxrejects <= max_maxrejects):
        op_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
        for primer_id, primer_seq in query_primers.items():
            if (primer_id not in primers_w_perfect_match):
                op_fasta.write(">%s\n%s\n" % (primer_id, primer_seq))
        op_fasta.flush()

        ip_vsearch = tempfile.NamedTemporaryFile(mode='w+t', suffix=".vsearch", dir=tempdir, delete=True) 
        vsearch_cmd = ['vsearch', '--usearch_global', None, '-db', transcriptome_vsearch_udb,
                       '--id', '0.75', '--minwordmatches', '2', '--wordlength', '7', '--dbmask', 'none', '--qmask', 'none', '--strand', 'both',
                       '--userout', None, '--userfields', 'query+qstrand+qlo+qhi+qrow+target+tstrand+tilo+tihi+trow+mism+opens',
                       '--maxaccepts', '1000000', '--maxgaps', '1', '--maxdiffs', '4', '--maxrejects', '4096', '--threads', str(num_vsearch_threads),
                       '--maxseqlength', '1000000', '--query_cov', '0.93']

        vsearch_cmd[2] = op_fasta.name
        vsearch_cmd[18] = ip_vsearch.name
        vsearch_cmd[28] = str(maxrejects)
        maxrejects *= 2

        try:
            begin_time = datetime.datetime.now()
            begin_time_str = begin_time.strftime("%Y-%m-%d %H:%M:%S")
            print("\n\t\tINFO: starting vsearch (%s)..." % begin_time_str, file=sys.stderr, end='', flush=True)
            vsearch_output = check_output(vsearch_cmd, stderr=subprocess.STDOUT)
            end_time = datetime.datetime.now()
            end_time_str = end_time.strftime("%Y-%m-%d %H:%M:%S")
            print("done (%s)" % end_time_str, file=sys.stderr, flush=True)
            time_delta = end_time - begin_time
            seconds = time_delta.total_seconds() 
            print("\t\tINFO: runtime = %d hr, %d min, %d sec" % (int(seconds / (3600)), int(seconds / 60), seconds % 60), file=sys.stderr, flush=True)
        except CalledProcessError as cpe:
            pdb.set_trace()

        ip_vsearch.flush()
        ip_vsearch.seek(0)
        for line in ip_vsearch:
            fields = line.strip().split('\t')
            all_vsearch_results.add(tuple(fields))
            if (fields[-2]=='0' and fields[-1]=='0'):
                primers_w_perfect_match.add(fields[0])

        perc_w_perfect_match = 100.0 * float(len(primers_w_perfect_match))/float(len(query_primers.keys()))
        print("\t\tINFO: matched %d of %d (%5.2f%%)" % (len(primers_w_perfect_match), len(query_primers), perc_w_perfect_match), file=sys.stderr, flush=True)

        op_fasta.close()
        ip_vsearch.close()

    primers_w_no_perfect_match = set(query_primers.keys()) - primers_w_perfect_match

    return (all_vsearch_results, primers_w_no_perfect_match)


def groupVsearchResultsByPrimerPair(all_vsearch_results, oligo_thermo):
    '''Query sequence IDs must be in a specific format. See below. '''
    unverified_pp_grouped_by_ppr = defaultdict(list)
    hyb_events_per_primer = defaultdict(list)  # For identifying promiscuous binders
    frac_duplexed_cache = {}

    for query, qstrand, qlo, qhi, qrow, target, tstrand, tilo, tihi, trow, mism, opens in all_vsearch_results:
        assert (tstrand == '+'), "Template match is negative strand. How?"
        match_fwdrev = "fwd" if (qstrand=='+') else "rev"
        query_len = int(qlo) if (qstrand=='-') else int(qhi)

        if (opens == '0'):
            query_seq, template_seq = qrow, trow
        else:
            query_seq = qrow.replace('-','')
            template_seq = trow.replace('-','')

        if (len(query_seq) == query_len): # Skip partial matches at the 5'/3' ends of mRNAs
            src, primer_fwdrev = query.rsplit('_',1)
            ppr_ID, pp_num = src.rsplit('.',1)
            pp_num = int(pp_num)
            frac_duplexed = 0

            cache_key = (qstrand, qrow, trow)
            if (cache_key in frac_duplexed_cache):
                frac_duplexed, query_seq, template_seq = frac_duplexed_cache[cache_key]
                if (frac_duplexed >= 0.01):
                    data_tup = (primer_fwdrev, match_fwdrev, query_seq, target, int(tilo), int(tihi), template_seq, opens != "0", frac_duplexed)
                    unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)].append(data_tup)
                    hyb_events_per_primer[(ppr_ID, pp_num)].append((primer_fwdrev, target, tilo, tihi, frac_duplexed))

            else:
                if (qstrand == '+'):
                    template_seq = template_seq[::-1].translate(DNA_complement_table)
                    any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                    end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                    min_dG = min(any_hyb.dg, end_hyb.dg)
                    frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)
                    
                else:
                    query_seq = query_seq[::-1].translate(DNA_complement_table)
                    any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                    end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                    min_dG = min(any_hyb.dg, end_hyb.dg)
                    frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)
                        
                if (frac_duplexed >= 0.01):
                    data_tup = (primer_fwdrev, match_fwdrev, query_seq, target, int(tilo), int(tihi), template_seq, opens != "0", frac_duplexed)
                    unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)].append(data_tup)
                    hyb_events_per_primer[(ppr_ID, pp_num)].append((primer_fwdrev, target, tilo, tihi, frac_duplexed))

                frac_duplexed_cache[cache_key] = (frac_duplexed, query_seq, template_seq)

    return unverified_pp_grouped_by_ppr, hyb_events_per_primer


def filterPrimerPairsByGlobalSpecificityAndPromiscuity(ppr_ID, igs_tuples_by_ID, query_primers, primer_pairs, unverified_pp_grouped_by_ppr, hyb_events_per_primer,
                                                       complete_TIG_IDs_w_nonlocal_members, oligo_thermo, max_promiscuous, transcriptome_ref):
    '''Adds total promiscuity counts to primer pair tuples that pass the filtering.'''
    igs_tuples_isoform_IDs = set(chain.from_iterable(igs_tuples_by_ID))
    verified_primer_pairs = []
    for pp_num, primer_pair in enumerate(primer_pairs):
        # Skip primer pair if the fwd and/or rev primer binds too many other non-target trancripts, or does not bind perfectly to any of its expected targets
        num_fwd_ontarget = len(list(filter(lambda x:x[0]=="fwd" and x[1] in igs_tuples_isoform_IDs, hyb_events_per_primer[(ppr_ID,pp_num)])))
        fwd_promiscuous = sum(list(map(itemgetter(4), filter(lambda x:x[0]=="fwd" and x[1] not in igs_tuples_isoform_IDs, hyb_events_per_primer[(ppr_ID,pp_num)]))))
        num_rev_ontarget = len(list(filter(lambda x:x[0]=="rev" and x[1] in igs_tuples_isoform_IDs, hyb_events_per_primer[(ppr_ID,pp_num)])))
        rev_promiscuous = sum(list(map(itemgetter(4), filter(lambda x:x[0]=="rev" and x[1] not in igs_tuples_isoform_IDs, hyb_events_per_primer[(ppr_ID,pp_num)]))))

        if (num_fwd_ontarget == 0 or num_rev_ontarget == 0 or fwd_promiscuous > max_promiscuous or rev_promiscuous > max_promiscuous):
            # case '== 0' indicates that the primer had so many near matches during transcriptome-wide
            # search that its exact match location couldn't be found. So, very promiscuous.
            continue

        blastn_results = unverified_pp_grouped_by_ppr[(ppr_ID, pp_num)]
        pp_was_verified, new_igs_tuples_by_ID = \
            verifyIsoformGroupsAreDistinctBySequencing_Blast(blastn_results, query_primers, oligo_thermo, igs_tuples_by_ID, transcriptome_ref, complete_TIG_IDs_w_nonlocal_members)
        if (pp_was_verified):
            aug_primer_pair = primer_pair + (fwd_promiscuous+rev_promiscuous, new_igs_tuples_by_ID)
            verified_primer_pairs.append( aug_primer_pair )
    
    return verified_primer_pairs


# Possible TODO: wrap this into its own class with other associated vsearch functionality used below
def verifyIsoformGroupsAreDistinctBySequencing_Blast(all_primer_matches, query_primers, oligo_thermo, expected_igs_tuple_by_ID, transcriptome_ref,
                                                     complete_TIG_IDs_w_nonlocal_members, disallowed_overlaps_on=set()):
    '''Verifies for a primer pair, given the transcriptome-wide matching of the fwd & rev primer in all_primer_matches.'''
    # List of tuples (expected ig, expected ig with nonlocal isoforms to form global ig, offtarget ig)
    new_igs_tuples_by_ID = ()

    fwd_primer_matches = defaultdict(list)
    rev_primer_matches = defaultdict(list)
    isoforms_grouped_by_paired_reads = defaultdict(list)
    read_len = designParams.read_len
    truth_val = True
    
    for (qseqid, primer_fwdrev, match_fwdrev, primer_seq, target, tilo, tihi, template_seq, has_gap, primer_3p_nonaln, frac_duplexed) in all_primer_matches:
        if (match_fwdrev == "fwd"):
            fwd_primer_matches[target].append( ("fwd", query_primers[qseqid], tilo, tihi, primer_3p_nonaln, frac_duplexed) )
        else:
            rev_primer_matches[target].append( ("rev", query_primers[qseqid], tilo, tihi, primer_3p_nonaln, frac_duplexed) )

    primers_product_lens = {}
    for isoform_ID in set(fwd_primer_matches.keys()) & set(rev_primer_matches.keys()):
        isoform_primer_matches = fwd_primer_matches[isoform_ID] + rev_primer_matches[isoform_ID]
        isoform_primer_matches.sort(key=itemgetter(2,3))
        amplicons_5p3p = []
        for tup_5p in filter(lambda x:x[0]=="fwd", isoform_primer_matches):
            for tup_3p in filter(lambda x:x[0]=="rev", isoform_primer_matches):
                len_primers_product = len(tup_5p[1]) + len(tup_3p[1]) + (tup_3p[3]-tup_3p[4]) - (tup_5p[3]+tup_5p[4]) - 1
                if (len_primers_product <= 10):
                    continue
                elif (len_primers_product < designParams.amplicon_min_len):
                    truth_val = False
                    break
                elif (len_primers_product <= designParams.amplicon_max_len and tup_5p[5] * tup_3p[5] > 0.01): # pair_frac_duplexed
                    inter_primer_seq = transcriptome_ref[isoform_ID][tup_5p[3]+tup_5p[4] : tup_3p[3]-1-tup_3p[4]]
                    primers_product = "%s%s%s" % (tup_5p[1], inter_primer_seq, tup_3p[1][::-1].translate(DNA_complement_table))
                    primers_product_lens[isoform_ID] = len(primers_product)
                    if (len(primers_product) != len_primers_product):
                        pdb.set_trace()
                            
                    amplicons_5p3p.append( (tup_5p[2],tup_3p[3]) )
                    paired_reads = (primers_product[0:read_len], primers_product[-read_len:])
                    isoforms_grouped_by_paired_reads[paired_reads].append( isoform_ID )

        if (not truth_val):
            break
        elif (isoform_ID in disallowed_overlaps_on):
            truth_val = all(map(lambda xy: xy[1][1]<xy[0][0] or xy[1][0]>xy[0][1], combinations(amplicons_5p3p,2)))
            if (not truth_val):
                print("\t\tFound overlapping amplicons on %s" % isoform_ID, file=sys.stderr)
                print(list(filter(lambda xy: not(xy[1][1]<xy[0][0] or xy[1][0]>xy[0][1]), combinations(amplicons_5p3p,2))), file=sys.stderr)
                break
            
    if (truth_val): # Overlapping amplicons not found on any amplicon for which overlaps are disallowed
        primed_igs_tuple_by_ID = set()
        for isoform_IDs in isoforms_grouped_by_paired_reads.values():
            isoform_IDs.sort()
            primed_igs_tuple_by_ID.add(tuple(isoform_IDs))

        truth_val = expected_igs_tuple_by_ID <= primed_igs_tuple_by_ID
       
        if (truth_val):
            # The expected isoform groups were amplified, in addition to maybe some offtarget ones
            new_igs_tuples_by_ID = tuple(sorted(primed_igs_tuple_by_ID.difference(expected_igs_tuple_by_ID)))
        else:
            # If truth_val is False, it may be because one or more of the igs in the igs tuple is/are a partial igs that has (identical) non-local isoform members.
            # Check whether the incoporation of those members would affect the truth_val.
            truth_for_each_ig = []
            primed_igs_tuple_by_ID_used = set()
            for expected_ig_by_ID in expected_igs_tuple_by_ID:
                truth_for_expected_ig = False
                if (expected_ig_by_ID in primed_igs_tuple_by_ID): 
                    # This particular ig (ie expected_ig_by_ID) was one of the igs found, so it is okay and the
                    # issues lies with one of the other igs
                    primed_igs_tuple_by_ID_used.add(expected_ig_by_ID)
                    truth_for_expected_ig = True
                elif (expected_ig_by_ID in complete_TIG_IDs_w_nonlocal_members):
                    print("In elif clause in verifyIsoformGroupsAreDistinctBySequencing2()", file=sys.stderr, flush=True)
                    if (complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID] in primed_igs_tuple_by_ID):
                        print("CHECK:", complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID], "in", primed_igs_tuple_by_ID, file=sys.stderr, flush=True)
                        #categorized_isoform_groups.append( (expected_ig_by_ID, complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID], None) )
                        primed_igs_tuple_by_ID_used.add(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID])
                        truth_for_expected_ig = True
                    elif (any(map(lambda ig: set(ig).issubset(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID]), primed_igs_tuple_by_ID))):
                        # The nonlocal members of a TIG are included because they share an identical pORF. Some of the members
                        # may not, though, have identical nucleotide sequences. This logic clause should catch such cases where
                        # the current primers amplify the local TIG isoforms and the nonlocal TIG isoforms, but in more than one grouping.
                        partial_nonlocal_TIGs = []
                        truth_for_expected_ig = True
                        for ig in filter(lambda x: x not in primed_igs_tuple_by_ID_used, primed_igs_tuple_by_ID):
                            if (set(ig).issubset(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID])):
                                partial_nonlocal_TIGs.append(ig)
                                primed_igs_tuple_by_ID_used.add(ig)
                            elif (len(set(ig) & set(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID])) > 0):
                                truth_for_expected_ig = False
                                break
                        #categorized_isoform_groups.append( (expected_ig_by_ID, partial_nonlocal_TIGs, None) )

                truth_for_each_ig.append(truth_for_expected_ig)

            truth_val = all(truth_for_each_ig)
            if (truth_val):
                new_igs_tuples_by_ID = primed_igs_tuple_by_ID.difference(primed_igs_tuple_by_ID_used)

    return (truth_val, new_igs_tuples_by_ID)


# Possible TODO: wrap this into its own class with other associated vsearch functionality used below
def verifyIsoformGroupsAreDistinctBySequencing2(all_primer_matches, oligo_thermo, expected_igs_tuple_by_ID, transcriptome_ref,
                                                complete_TIG_IDs_w_nonlocal_members, disallowed_overlaps_on=set()):
    '''Verifies for a primer pair, given the transcriptome-wide matching of the fwd & rev primer in all_primer_matches.'''
    # List of tuples (expected ig, expected ig with nonlocal isoforms to form global ig, offtarget ig)
    new_igs_tuples_by_ID = ()

    fwd_primer_matches = defaultdict(list)
    rev_primer_matches = defaultdict(list)
    isoforms_grouped_by_paired_reads = defaultdict(list)
    read_len = designParams.read_len
    truth_val = True
    
    for (primer_fwdrev, match_fwdrev, primer_seq, target, tilo, tihi, template_seq, has_gap, frac_duplexed) in all_primer_matches:
        if (match_fwdrev == "fwd"):
            fwd_primer_matches[target].append( ("fwd", primer_seq, tilo, tihi, frac_duplexed) )
        else:
            rev_primer_matches[target].append( ("rev", primer_seq, tilo, tihi, frac_duplexed) )

    primers_product_lens = {}
    for isoform_ID in set(fwd_primer_matches.keys()) & set(rev_primer_matches.keys()):
        isoform_primer_matches = fwd_primer_matches[isoform_ID] + rev_primer_matches[isoform_ID]
        isoform_primer_matches.sort(key=itemgetter(2,3))
        amplicons_5p3p = []
        for tup_5p in filter(lambda x:x[0]=="fwd", isoform_primer_matches):
            for tup_3p in filter(lambda x:x[0]=="rev", isoform_primer_matches):
                pair_frac_duplexed = tup_5p[4] * tup_3p[4]
                if (tup_3p[3]-tup_5p[3] > 10 and pair_frac_duplexed > 0.01):
                    inter_primer_seq = transcriptome_ref[isoform_ID][tup_5p[3]:tup_3p[2]-1]
                    primers_product = "%s%s%s" % (tup_5p[1], inter_primer_seq, tup_3p[1][::-1].translate(DNA_complement_table))
                    primers_product_lens[isoform_ID] = len(primers_product)
                    if (len(primers_product) < designParams.amplicon_min_len):
                        truth_val = False
                        break
                    elif (len(primers_product) <= designParams.amplicon_max_len):
                        amplicons_5p3p.append( (tup_5p[2],tup_5p[4],tup_3p[3],tup_3p[4]) )
                        paired_reads = (primers_product[0:read_len], primers_product[-read_len:])
                        isoforms_grouped_by_paired_reads[paired_reads].append( isoform_ID )

        if (not truth_val):
            break
        elif (isoform_ID in disallowed_overlaps_on):
            truth_val = all(map(lambda xy: xy[1][2]<xy[0][0] or xy[1][0]>xy[0][2], combinations(amplicons_5p3p,2)))
            if (not truth_val):
                print("\t\tFound overlapping amplicons on %s" % isoform_ID, file=sys.stderr)
                print(list(filter(lambda xy: not(xy[1][2]<xy[0][0] or xy[1][0]>xy[0][2]), combinations(amplicons_5p3p,2))), file=sys.stderr)
                break
            
    if (truth_val): # Overlapping amplicons not found on any amplicon for which overlaps are disallowed
        primed_igs_tuple_by_ID = set()
        for isoform_IDs in isoforms_grouped_by_paired_reads.values():
            isoform_IDs.sort()
            primed_igs_tuple_by_ID.add(tuple(isoform_IDs))

        truth_val = expected_igs_tuple_by_ID <= primed_igs_tuple_by_ID
       
        if (truth_val):
            # The expected isoform groups were amplified, in addition to maybe some offtarget ones
            new_igs_tuples_by_ID = tuple(sorted(primed_igs_tuple_by_ID.difference(expected_igs_tuple_by_ID)))
        else:
            # If truth_val is False, it may be because one or more of the igs in the igs tuple is/are a partial igs that has (identical) non-local isoform members.
            # Check whether the incoporation of those members would affect the truth_val.
            truth_for_each_ig = []
            primed_igs_tuple_by_ID_used = set()
            for expected_ig_by_ID in expected_igs_tuple_by_ID:
                truth_for_expected_ig = False
                if (expected_ig_by_ID in primed_igs_tuple_by_ID): 
                    # This particular ig (ie expected_ig_by_ID) was one of the igs found, so it is okay and the
                    # issues lies with one of the other igs
                    primed_igs_tuple_by_ID_used.add(expected_ig_by_ID)
                    truth_for_expected_ig = True
                elif (expected_ig_by_ID in complete_TIG_IDs_w_nonlocal_members):
                    print("In elif clause in verifyIsoformGroupsAreDistinctBySequencing2()", file=sys.stderr, flush=True)
                    if (complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID] in primed_igs_tuple_by_ID):
                        print("CHECK:", complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID], "in", primed_igs_tuple_by_ID, file=sys.stderr, flush=True)
                        #categorized_isoform_groups.append( (expected_ig_by_ID, complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID], None) )
                        primed_igs_tuple_by_ID_used.add(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID])
                        truth_for_expected_ig = True
                    elif (any(map(lambda ig: set(ig).issubset(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID]), primed_igs_tuple_by_ID))):
                        # The nonlocal members of a TIG are included because they share an identical pORF. Some of the members
                        # may not, though, have identical nucleotide sequences. This logic clause should catch such cases where
                        # the current primers amplify the local TIG isoforms and the nonlocal TIG isoforms, but in more than one grouping.
                        partial_nonlocal_TIGs = []
                        truth_for_expected_ig = True
                        for ig in filter(lambda x: x not in primed_igs_tuple_by_ID_used, primed_igs_tuple_by_ID):
                            if (set(ig).issubset(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID])):
                                partial_nonlocal_TIGs.append(ig)
                                primed_igs_tuple_by_ID_used.add(ig)
                            elif (len(set(ig) & set(complete_TIG_IDs_w_nonlocal_members[expected_ig_by_ID])) > 0):
                                truth_for_expected_ig = False
                                break
                        #categorized_isoform_groups.append( (expected_ig_by_ID, partial_nonlocal_TIGs, None) )

                truth_for_each_ig.append(truth_for_expected_ig)

            truth_val = all(truth_for_each_ig)
            if (truth_val):
                new_igs_tuples_by_ID = primed_igs_tuple_by_ID.difference(primed_igs_tuple_by_ID_used)

    return (truth_val, new_igs_tuples_by_ID)
