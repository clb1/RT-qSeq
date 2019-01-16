#!/usr/bin/env python3

import bz2
from collections import defaultdict, Counter
import pickle
import datetime
from itertools import combinations, product
import sqlite3
from operator import itemgetter
import os
import sys
import tempfile
import designParams
from subprocess import check_output, CalledProcessError

# TODO: These are included just because of localCalcFracDuplexed(), which itself is just a patch to allow parallization.
from math import exp
import primer3

from xlrd import open_workbook
from xlutils.copy import copy
from xlwt import Workbook, XFStyle, easyxf, Formula

import networkx as nx

from designSetsClasses import SetOfPrimers, PrimerSingleton, OligoThermodynamics, runThoroughVsearch
from multiprocessing import Pool

from pyfaidx import Fasta

DNA_complement_table = str.maketrans("ACGTNacgtn","TGCANtgcan")

import pdb

def collectSetsOfPrimers(database_path, params_descriptor, target_chrom_OS):
    # Collect all primers to align against transcriptome
    all_sets_of_primers = {}
    primers_to_align = set()
    primer_ids_per_primer_seq = defaultdict(set)
    all_target_isoform_IDs = set()
    loci_isoform_IDs_per_seq = defaultdict(set)
    
    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    ip = open(target_chrom_OS, 'r')
    for line in ip:
        chromosome, olap_set_ID = line.strip().split('\t')
        
        # Check whether designs for overlap set have already been completed
        selection_tuple = (chromosome, olap_set_ID, designParams.params_descriptor)
        cursor.execute("SELECT solns_pkl_file FROM ExperimentDesignMeta WHERE chrom = ? AND OS_ID = ? and params_descriptor = ?", selection_tuple)
        pkl_files = cursor.fetchall()
        if (len(pkl_files)==0):  # TEMPORARY UNTIL CompleteDesign problems fixed and all loci are made
            print("WARNING: no pickle file for %s %s" % (chromosome, olap_set_ID), file=sys.stderr, flush=True)
            continue
        if (pkl_files[0][0] == ''):
            print("WARNING: no pickle file for %s %s" % (chromosome, olap_set_ID), file=sys.stderr, flush=True)
            continue
        assert (len(pkl_files)==1), "Multiple pickle files for %s %s %s" % selection_tuple

        try:
            ppr_combo_IDs_and_pp_sets = pickle.load(open(pkl_files[0][0], 'rb'))
        except IOError as ioe:
            continue  #TODO: remove, temporary
            pdb.set_trace()
            
        for ppr_combo_IDs, pp_sets in ppr_combo_IDs_and_pp_sets:
            ppr_combo_IDs_descr = "+".join(ppr_combo_IDs)
            for pp_set_num, sop in enumerate(pp_sets):
                sop_target_isoform_IDs = sop.getGlobalTargetIsoformIDs()
                all_target_isoform_IDs.update(sop_target_isoform_IDs)
                all_sets_of_primers[(chromosome, olap_set_ID, ppr_combo_IDs_descr, pp_set_num)] = sop
                for fwd_primer in sop.getAllFwdPrimers():
                    fwd_primer_descr = (chromosome, olap_set_ID, ppr_combo_IDs_descr, pp_set_num, "fwd")
                    primer_ids_per_primer_seq[fwd_primer.seq].add(fwd_primer_descr)
                    primers_to_align.add(fwd_primer.seq)
                    loci_isoform_IDs_per_seq[fwd_primer.seq].update(sop_target_isoform_IDs)
                    
                for rev_primer in sop.getAllRevPrimers():
                    rev_primer_descr = (chromosome, olap_set_ID, ppr_combo_IDs_descr, pp_set_num, "rev")
                    primer_ids_per_primer_seq[rev_primer.seq].add(rev_primer_descr)
                    primers_to_align.add(rev_primer.seq)
                    loci_isoform_IDs_per_seq[rev_primer.seq].update(sop_target_isoform_IDs)
    conn.close()

    return (all_sets_of_primers, primers_to_align, primer_ids_per_primer_seq, all_target_isoform_IDs, loci_isoform_IDs_per_seq)


def alignPrimersAgainstTranscriptome(oligo_thermo, primers_to_align, transcriptome_vsearch_udb, tempdir, all_target_isoform_IDs, loci_isoform_IDs_per_seq):
    results_per_primer_seq = defaultdict(set)
    do_not_use = set()
    
    query_primers = dict(map(lambda x: (x,x), primers_to_align))
    all_vsearch_results, primers_w_no_perfect_match = runThoroughVsearch(query_primers, transcriptome_vsearch_udb, tempdir, max_maxrejects = 32768, num_vsearch_threads = 20)

    for query, qstrand, qlo, qhi, qrow, target, tstrand, tilo, tihi, trow, mism, opens in all_vsearch_results:
        match_fwdrev = "fwd" if (qstrand=='+') else "rev"
        query_len = int(qlo) if (qstrand=='-') else int(qhi)

        if (opens == '0'):
            query_seq, template_seq = qrow, trow
        else:
            query_seq = qrow.replace('-','')
            template_seq = trow.replace('-','')

        if (len(query_seq) == query_len): # Skip partial matches at the 5'/3' ends of mRNAs
            primer_fwdrev = "fwd" if (qstrand == '+') else "rev"
            pos_5p = int(tilo) if (qstrand == '+') else int(tihi)

            frac_duplexed = 0

            # indels in the 3' end or more than 1bp mismatch -> effectively no primer extension
            if (qstrand == '+' and '-' not in qrow[-6:] and '-' not in trow[-6:] and sum([c1!=c2 for c1,c2 in zip(qrow[-6:],trow[-6:])]) <= 2):
                template_seq = template_seq[::-1].translate(DNA_complement_table)
                any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                min_dG = min(any_hyb.dg, end_hyb.dg)
                frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)
                    
            elif (qstrand == '-' and '-' not in qrow[0:6] and '-' not in trow[0:6] and sum([c1!=c2 for c1,c2 in zip(qrow[0:6],trow[0:6])]) <= 2):
                query_seq = query_seq[::-1].translate(DNA_complement_table)
                any_hyb = oligo_thermo.calcHeterodimer(query_seq, template_seq)
                end_hyb = oligo_thermo.calcEndStability(query_seq, template_seq)
                min_dG = min(any_hyb.dg, end_hyb.dg)
                frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)

            if (frac_duplexed >= 0.01):
                results_per_primer_seq[query_seq].add((primer_fwdrev, match_fwdrev, target, pos_5p, frac_duplexed))
                if (target in loci_isoform_IDs_per_seq[query] and (opens != '0' or mism != '0')):  # TODO: Temporary until fixed in SOP creation
                    do_not_use.add(query)

    return (results_per_primer_seq, primers_w_no_perfect_match, do_not_use)


def cullPrimerSets(all_sets_of_primers, do_not_use):
    sop_to_cull = set()
    for sop_descr, sop in all_sets_of_primers.items():
        sop_primer_seqs = set()
        for primer in sop.getAllFwdPrimers():
            sop_primer_seqs.add(primer.seq)
        for primer in sop.getAllRevPrimers():
            sop_primer_seqs.add(primer.seq)

        if (len(sop_primer_seqs & do_not_use) > 0):
            sop_to_cull.add(sop_descr)

    print("INFO: culling %d primer sets" % len(sop_to_cull), file=sys.stderr, flush=True)
    for sop_descr in sop_to_cull:
        del all_sets_of_primers[sop_descr]

        
def flagPositionallyConflictingPrimers(results_per_primer_seq, primer_ids_per_primer_seq, all_target_isoform_IDs):
    primer_matches_per_isoform = defaultdict(list)
    for primer_seq in primer_ids_per_primer_seq.keys():
        assert (len(results_per_primer_seq[primer_seq]) > 0), "No vsearch results for primer sequence %s" % primer_seq
        for (primer_fwdrev, match_fwdrev, target_isoform_ID, pos_5p, frac_duplexed) in results_per_primer_seq[primer_seq]:
            for (chrom, OS_ID, ppr_combo_descr, pp_set_num, primer_fwdrev) in primer_ids_per_primer_seq[primer_seq]:
                short_primer_ID = "%s_%s" % (chrom, OS_ID)
                long_primer_ID = "%s_%s_%s_%s" % (chrom, OS_ID, ppr_combo_descr, pp_set_num)
                primer_matches_per_isoform[target_isoform_ID].append( (short_primer_ID, long_primer_ID, pos_5p, frac_duplexed, primer_seq, primer_fwdrev, match_fwdrev) )
                

    incompatible_primers = defaultdict(set)
    for isoform_ID, match_data in primer_matches_per_isoform.items():
        match_data.sort(key=itemgetter(2))

        for ml, mr in combinations(match_data,2): # left and right primer matches on transcript
            subseq_len = mr[2]-ml[2]
            if (ml[1]==mr[1]): # From same primers set
                try:
                    assert ((ml[-2] == ml[-1] and mr[-2] == mr[-1]) or (ml[-1] == "rev" and mr[-1] == "fwd") or
                            (ml[-1] == "fwd" and mr[-1] == "rev" and subseq_len > designParams.confusion_amplicon_max_len)), \
                            "Unexpected primer match orientation between primers from same primers set"
                except AssertionError:
                    print >> sys.stderr, ml
                    print >> sys.stderr, mr
                    print >> sys.stderr, "-------"
                    pdb.set_trace()
            elif (ml[0] != mr[0]): # ml[0]==mr[0] means from same overlap set but different primer set (the 'if' clause catches same primer set), 
                                   # and so are guaranteed to be used mutually exclusively anyhow.
                ml_behavior = (ml[-2],ml[-1])
                mr_behavior = (mr[-2],mr[-1])

                # A primer could act as fwd/rev even if it is designed to be rev/fwd
                if (ml_behavior == ("fwd","fwd") and mr_behavior == ("rev","rev") and subseq_len > 10 and
                    (isoform_ID in all_target_isoform_IDs or subseq_len < designParams.confusion_amplicon_max_len)):
                    incompatible_primers["fwdfwd-revrev"].add( (ml[4], mr[4]) )

                elif (ml_behavior == mr_behavior == ("fwd","fwd") and isoform_ID in all_target_isoform_IDs and
                      subseq_len < designParams.confusion_amplicon_max_len):
                    incompatible_primers["fwdfwd-fwdfwd"].add( (ml[4], mr[4]) )

                elif (ml_behavior == mr_behavior == ("rev","rev") and isoform_ID in all_target_isoform_IDs):
                    incompatible_primers["revrev-revrev"].add( (ml[4], mr[4]) )

                elif (ml_behavior == ("fwd","fwd") and mr_behavior == ("fwd","rev") and subseq_len > 10 and
                      (isoform_ID in all_target_isoform_IDs or subseq_len < designParams.confusion_amplicon_max_len)):
                    incompatible_primers["fwdfwd-fwdrev"].add( (ml[4], mr[4]) )

                elif (ml_behavior == ("rev","fwd") and mr_behavior == ("rev","rev") and subseq_len > 10 and
                      (isoform_ID in all_target_isoform_IDs or subseq_len < designParams.confusion_amplicon_max_len)):
                    incompatible_primers["revfwd-revrev"].add( (ml[4], mr[4]) )

                elif (ml_behavior == ("rev","fwd") and mr_behavior == ("fwd","rev") and subseq_len > 10 and
                      (isoform_ID in all_target_isoform_IDs or subseq_len < designParams.confusion_amplicon_max_len)):
                    incompatible_primers["revfwd-fwdrev"].add( (ml[4], mr[4]) )

                elif (ml_behavior == ("rev","rev") and mr_behavior == ("fwd","rev") and
                      (isoform_ID in all_target_isoform_IDs or subseq_len < designParams.confusion_amplicon_max_len)):
                    incompatible_primers["revrev-fwdrev"].add( (ml[4], mr[4]) )

    print("INFO: found %d positionally conflicting primer pair interactionss:" % sum(map(len, incompatible_primers.values())), file=sys.stderr, flush=True)
    for conflict_type, conflicting_pairs in incompatible_primers.items():
        print("\t%s\t%d" % (conflict_type, len(conflicting_pairs)), file=sys.stderr, flush=True)

    # Add converse, for ease-of-lookup later
    for behavior, seq_pairs in incompatible_primers.items():
        incompatible_primers[behavior] |= set(map(lambda seqs: (seqs[1],seqs[0]), seq_pairs))

    return incompatible_primers


# TODO: This function is just a patch because oligo_thermo can't be serialized to be used in parallization.
# Try using pathos package.
def localCalcFracDuplexed(seq1, seq2):
    rxn_temp_C = designParams.DNAPOL_rxn_temp
    rxn_temp_K = 273.15 + rxn_temp_C
    thermo_analysis = primer3.thermoanalysis.ThermoAnalysis(mv_conc=designParams.monovalent_salt_conc, dv_conc=designParams.divalent_salt_conc,
                                                            dntp_conc=designParams.dntp_conc, temp_c=rxn_temp_C, dna_conc=designParams.input_primer_conc_nM)

    RT = 1.98717 * rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)

    seq1seq2_any_hyb = thermo_analysis.calcHeterodimer(seq1, seq2)
    seq1seq2_end_hyb = thermo_analysis.calcEndStability(seq1, seq2)
    seq2seq1_end_hyb = thermo_analysis.calcEndStability(seq2, seq1)

    # dG units are expected to be cal/mol and NOT kcal/mol
    end_dG = sum(set(filter(lambda x: x<0, [seq1seq2_end_hyb.dg, seq2seq1_end_hyb.dg])))
    any_dG = seq1seq2_any_hyb.dg

    frac_end_duplexed = (designParams.input_primer_conc_M * exp(-end_dG/RT))/(1 + designParams.input_primer_conc_M * exp(-end_dG/RT))
    frac_any_duplexed = (designParams.input_primer_conc_M * exp(-any_dG/RT))/(1 + designParams.input_primer_conc_M * exp(-any_dG/RT))

    return (seq1, seq2, frac_end_duplexed, frac_any_duplexed)



def flagThermodynamicallyConflictingPrimers(oligo_thermo, primer_seqs):
    candidate_interacting_primers = set()
    thermo_interacting_primers = {}
    mer_size = 5

    # All vs all is intractable, so instead group sequences that share a kmer and the
    # reverse complement of the kmer and only evaluate all pairs within a group
    canonical_and_rc_mers = []
    for mer in map(''.join, product('ACGT', repeat=mer_size)):
        mer_rc = mer[::-1].translate(DNA_complement_table)
        canonical_mer, canonical_mer_rc = (mer, mer_rc) if (mer < mer_rc) else (mer_rc, mer)
        canonical_and_rc_mers.append( (canonical_mer, canonical_mer_rc) )

    mer_groups = defaultdict(set)
    for seq in primer_seqs:
        for mer in [seq[i:i+mer_size] for i in range(len(seq)-mer_size+1)]:
            mer_groups[mer].add(seq)
                
    for mer, mer_rc in canonical_and_rc_mers:
        if (mer in mer_groups and mer_rc in mer_groups):
            mer_seqs = mer_groups[mer]
            mer_rc_seqs = mer_groups[mer_rc]

            for seq1, seq2 in product(mer_seqs, mer_rc_seqs):
                if ((seq1,seq2) not in candidate_interacting_primers and (seq2,seq1) not in candidate_interacting_primers):
                    candidate_interacting_primers.add( (seq1,seq2) )
    print("%d allowed candidate primer interactions" % len(candidate_interacting_primers), file=sys.stderr, flush=True)

    with Pool(processes=19) as pool:
        #inputs = map(lambda s: s+(False,), list(candidate_interacting_primers)[0:100])
        #for (seq1,seq2,frac_duplexed) in pool.starmap(oligo_thermo.calcPrimerPairThermoDetails, inputs):
        for (seq1, seq2, frac_end_duplexed, frac_any_duplexed) in pool.starmap(localCalcFracDuplexed, candidate_interacting_primers):
            if (frac_end_duplexed > 0.01 or frac_any_duplexed > 0.01):
                thermo_interacting_primers[(seq1,seq2)] = (frac_end_duplexed, frac_any_duplexed)
                thermo_interacting_primers[(seq2,seq1)] = (frac_end_duplexed, frac_any_duplexed)
    
    return thermo_interacting_primers


def makeCompatibilityGraphForSetsOfPrimers(all_sets_of_primers, incompatible_primers_by_pos, thermo_interacting_primers):
    G = nx.Graph()
    
    incompatible_seqs_by_pos = set()
    for incompat_type in ["fwdfwd-revrev", "fwdfwd-fwdfwd", "revrev-revrev", "fwdfwd-fwdrev"]:
        incompatible_seqs_by_pos.update(incompatible_primers_by_pos[incompat_type])
    
    for (sop1_descr, sop1), (sop2_descr, sop2) in combinations(all_sets_of_primers.items(), 2):
        if (sop1_descr[0:2] != sop2_descr[0:2]):
            fwd_fwd_frac_duplexed = [0]
            fwd_rev_frac_duplexed = [0]
            rev_rev_frac_duplexed = [0]
            primers_pos_incompat = False

            for fwd_primer1, fwd_primer2 in product(sop1.getAllFwdPrimers(), sop2.getAllFwdPrimers()):
                if ((fwd_primer1.seq,fwd_primer2.seq) in incompatible_seqs_by_pos):
                    primers_pos_incompat = True
                    break

                if ((fwd_primer1.seq,fwd_primer2.seq) in thermo_interacting_primers):
                    frac_end_duplexed, frac_any_duplexed = thermo_interacting_primers[(fwd_primer1.seq,fwd_primer2.seq)]
                    fwd_fwd_frac_duplexed.append( max(frac_end_duplexed, frac_any_duplexed) )
            if (primers_pos_incompat):
                continue


            for fwd_primer1, rev_primer2 in product(sop1.getAllFwdPrimers(), sop2.getAllRevPrimers()):
                if ((fwd_primer1.seq,rev_primer2.seq) in incompatible_seqs_by_pos):
                    primers_pos_incompat = True
                    break
                
                if ((fwd_primer1.seq,rev_primer2.seq) in thermo_interacting_primers):
                    frac_end_duplexed, frac_any_duplexed = thermo_interacting_primers[(fwd_primer1.seq,rev_primer2.seq)]
                    fwd_rev_frac_duplexed.append( frac_end_duplexed )
            if (primers_pos_incompat):
                continue


            for rev_primer1, fwd_primer2 in product(sop1.getAllRevPrimers(), sop2.getAllFwdPrimers()):
                if ((rev_primer1.seq,fwd_primer2.seq) in incompatible_seqs_by_pos):
                    primers_pos_incompat = True
                    break
                
                if ((rev_primer1.seq,fwd_primer2.seq) in thermo_interacting_primers):
                    frac_end_duplexed, frac_any_duplexed = thermo_interacting_primers[(rev_primer1.seq,fwd_primer2.seq)]
                    fwd_rev_frac_duplexed.append( frac_end_duplexed  )
            if (primers_pos_incompat):
                continue


            for rev_primer1, rev_primer2 in product(sop1.getAllRevPrimers(), sop2.getAllRevPrimers()):
                if ((rev_primer1.seq,rev_primer2.seq) in incompatible_seqs_by_pos):
                    primers_pos_incompat = True
                    break
                
                if ((rev_primer1.seq,rev_primer2.seq) in thermo_interacting_primers):
                    frac_end_duplexed, frac_any_duplexed = thermo_interacting_primers[(rev_primer1.seq,rev_primer2.seq)]
                    rev_rev_frac_duplexed.append( max(frac_end_duplexed, frac_any_duplexed) )
            if (primers_pos_incompat):
                continue


            sum_fwd_fwd_frac_duplexed = sum(fwd_fwd_frac_duplexed)
            max_fwd_fwd_frac_duplexed = max(fwd_fwd_frac_duplexed)

            sum_fwd_rev_frac_duplexed = sum(fwd_rev_frac_duplexed)
            max_fwd_rev_frac_duplexed = max(fwd_rev_frac_duplexed)

            sum_rev_rev_frac_duplexed = sum(rev_rev_frac_duplexed)
            max_rev_rev_frac_duplexed = max(rev_rev_frac_duplexed)

            G.add_edge(sop1_descr, sop2_descr,
                       sum_fwd_fwd_frac_duplexed=sum_fwd_fwd_frac_duplexed, max_fwd_fwd_frac_duplexed=max_fwd_fwd_frac_duplexed,
                       sum_fwd_rev_frac_duplexed=sum_fwd_rev_frac_duplexed, max_fwd_rev_frac_duplexed=max_fwd_rev_frac_duplexed,
                       sum_rev_rev_frac_duplexed=sum_rev_rev_frac_duplexed, max_rev_rev_frac_duplexed=max_rev_rev_frac_duplexed)
    return G


def findACliqueFor(G, nodes_to_keep, index_to_node, tempdir, max_duplex_same_type = 0.01, max_duplex_opp_type = 0.5):
    subG = G.subgraph(nodes_to_keep)

    to_dimacs = {}
    from_dimacs = {}
    for dimacs_index, n in enumerate(subG.nodes(), 1):
        my_index = subG.node[n]['index']
        to_dimacs[my_index] = dimacs_index
        from_dimacs[dimacs_index] = my_index
        
    max_clique_dyn_basename = "iter_clique_build"
    output_edge_file = "%s/%s.edges" % (tempdir, max_clique_dyn_basename)
    op_edges = open(output_edge_file, 'wt')
    num_edges_written = 0
    for n1, n2, D in subG.edges(data=True):
        if (D['max_fwd_fwd_frac_duplexed'] < max_duplex_same_type and D['max_rev_rev_frac_duplexed'] < max_duplex_same_type and D['max_fwd_rev_frac_duplexed'] < max_duplex_opp_type):
            op_edges.write("e %d %d\n" % (to_dimacs[subG.node[n1]["index"]], to_dimacs[subG.node[n2]["index"]]))
            num_edges_written += 1
    op_edges.close()

    if (num_edges_written > 0):
        try:
            mcqd_output = check_output(["/usr/local/src/MaxCliquePara/maxClique", output_edge_file, "20", "1"])
            mcqd_output = str(mcqd_output,'utf-8')
        except CalledProcessError as cpe:
            pdb.set_trace()

        for line in mcqd_output.split('\n'):
            if (line.startswith("search took")):
                print("INFO: " + line.split(';')[0], file=sys.stderr, flush=True)
            elif (line.startswith("Clique")):
                summary_info, nodes_in_clique_raw = line.replace('[','').replace(']','').replace('(','').replace(')','').strip().split(':')
                summary_info = summary_info.replace("Clique", "INFO: Number of primer sets =")
                dimacs_vertices = map(int, nodes_in_clique_raw.strip().split(','))
                clique_node_indexes = list(map(lambda ind: from_dimacs[ind+1], dimacs_vertices))

                # Confirm clique
                clique_nodes = list(map(lambda x: index_to_node[x], clique_node_indexes))
                edges_exist = list(map(lambda nn: subG.has_edge(nn[0],nn[1]), combinations(clique_nodes,2)))
                try:
                    assert(all(edges_exist)), "Returned nodes do not form a clique"
                except AssertionError:
                    pdb.set_trace()

                print(summary_info, file=sys.stderr, flush=True)
            elif (line.startswith("Warning")):
                print(line, file=sys.stderr, flush=True)
    else:
        clique_node_indexes = []

    os.remove(output_edge_file)
    return clique_node_indexes


def formPoolIteratively(G, all_sets_of_primers, pool_type, tempdir):
    index_to_node = {}
    max_duplex_same_type = 0.25
    max_duplex_opp_type = 0.5

    nodes_by_num_primers = defaultdict(set)
    for node_index, n in enumerate(G.nodes(), 1):
        G.node[n]['index'] = node_index
        index_to_node[node_index] = n

        sop = all_sets_of_primers[n]
        num_fwd_primers = len(sop.getAllFwdPrimers())
        num_rev_primers = len(sop.getAllRevPrimers())
        num_primers = num_fwd_primers + num_rev_primers
        nodes_by_num_primers[num_primers].add(n)

    pdb.set_trace()
    min_primer_count = min(nodes_by_num_primers.keys())
    max_primer_count = max(nodes_by_num_primers.keys())
    #initial_primer_set =  # nodes_by_num_primers[max_primer_count] | 
    clique_node_indexes = findACliqueFor(G, nodes_by_num_primers[max_primer_count], index_to_node, tempdir, max_duplex_same_type, max_duplex_opp_type)
    
    nodes_in_clique = set(map(lambda node_index: index_to_node[node_index], clique_node_indexes))

    clique_sops = list(map(lambda sop_descr: all_sets_of_primers[sop_descr], nodes_in_clique))
    num_clique_fwd_primers, num_clique_rev_primers = 0, 0
    fwd_primer_counts, rev_primer_counts = [], []
    for sop in clique_sops:
        num_fwd_primers = len(sop.getAllFwdPrimers())
        num_rev_primers = len(sop.getAllRevPrimers())
        fwd_primer_counts.append(num_fwd_primers)
        rev_primer_counts.append(num_rev_primers)
        num_clique_fwd_primers += num_fwd_primers
        num_clique_rev_primers += num_rev_primers
    print("INFO: clique has %d fwd primers and %d rev primers" % (num_clique_fwd_primers, num_clique_rev_primers), file=sys.stderr, flush=True)
        
    remaining_primer_counts = sorted(nodes_by_num_primers.keys(), reverse=True)
    remaining_primer_counts.remove(max_primer_count)

    # Find the other primer nodes that are connected to all of the nodes in the clique
    for primer_count in remaining_primer_counts:
        completely_connected_nodes = set()
        for node in nodes_by_num_primers[primer_count]:
            if (all(map(lambda clique_node: G.has_edge(node,clique_node), nodes_in_clique))):
                completely_connected_nodes.add(node)
        print("%d of %d %d-primer nodes were compatible" % (len(completely_connected_nodes), len(nodes_by_num_primers[primer_count]), primer_count), file=sys.stderr, flush=True)

        # Find the largest clique in the connected nodes
        new_clique_node_indexes = findACliqueFor(G, completely_connected_nodes, index_to_node, tempdir, max_duplex_same_type, max_duplex_opp_type)
        new_nodes_in_clique = set(map(lambda node_index: index_to_node[node_index], new_clique_node_indexes))

        nodes_in_clique |= new_nodes_in_clique
        clique_sops = list(map(lambda sop_descr: all_sets_of_primers[sop_descr], nodes_in_clique))
        num_clique_fwd_primers, num_clique_rev_primers = 0, 0
        fwd_primer_counts, rev_primer_counts = [], []
        for sop in clique_sops:
            num_fwd_primers = len(sop.getAllFwdPrimers())
            num_rev_primers = len(sop.getAllRevPrimers())
            fwd_primer_counts.append(num_fwd_primers)
            rev_primer_counts.append(num_rev_primers)
            num_clique_fwd_primers += num_fwd_primers
            num_clique_rev_primers += num_rev_primers
        print("INFO: clique has %d fwd primers and %d rev primers" % (num_clique_fwd_primers, num_clique_rev_primers), file=sys.stderr, flush=True)

    return nodes_in_clique


def confirmIsoformGroupsCorrectlyPrimed_Local(clique_sop_descriptions, all_sets_of_primers, results_per_primer_seq, transcriptome_fasta):
    read_len = designParams.read_len
    transcriptome_ref = Fasta(transcriptome_fasta, as_raw=True, sequence_always_upper=True)

    # sop_descr = (chromosome, olap_set_ID, ppr_combo_IDs_descr, pp_set_num)
    pdb.set_trace()
    for sop_descr in clique_sop_descriptions:
        sop = all_sets_of_primers[sop_descr]
        #if ((chromosome, olap_set_ID, ppr_combo_IDs_descr, pp_set_num) != ('chr7', 'OS52', 'F40R20+F54R41+F53R50+F36R19', 0)):
        #    continue
        print(sop_descr, file=sys.stderr, flush=True)
        all_annot_primers = sop.getAllAnnotPrimers() 
        fwd_primers = filter(lambda t: t[0][0]=='F', all_annot_primers)
        rev_primers = filter(lambda t: t[0][0]=='R', all_annot_primers)

        isoforms_grouped_by_paired_reads = defaultdict(list)
        amplicons_per_paired_reads = defaultdict(set)
        for (fwd_erID, fwd_primer), (rev_erID, rev_primer) in product(fwd_primers, rev_primers):
            for (fwd_primer_fwdrev, fwd_match_fwdrev, fwd_target_isoform_ID, fwd_pos_5p, fwd_frac_duplexed) in results_per_primer_seq[fwd_primer.seq]:
                for (rev_primer_fwdrev, rev_match_fwdrev, rev_target_isoform_ID, rev_pos_5p, rev_frac_duplexed) in results_per_primer_seq[rev_primer.seq]:

                    if (fwd_target_isoform_ID==rev_target_isoform_ID and fwd_match_fwdrev=="fwd" and rev_match_fwdrev=="rev" and
                        rev_pos_5p - fwd_pos_5p > 10 and rev_pos_5p - fwd_pos_5p < designParams.confusion_amplicon_max_len):

                        #if (fwd_target_isoform_ID=='CG_60563.1.10.9'):
                        #    print(fwd_target_isoform_ID, fwd_primer.seq, rev_primer.seq, file=sys.stderr, flush=True)
                        #    print(fwd_primer_fwdrev, fwd_match_fwdrev, fwd_target_isoform_ID, fwd_pos_5p, fwd_frac_duplexed, file=sys.stderr, flush=True)
                        #    print(rev_primer_fwdrev, rev_match_fwdrev, rev_target_isoform_ID, rev_pos_5p, rev_frac_duplexed, file=sys.stderr, flush=True)
                        #    print(fwd_primer, file=sys.stderr, flush=True)
                        #    print(rev_primer, file=sys.stderr, flush=True)
                        #    print("-----", file=sys.stderr, flush=True)
                            
                        assert (fwd_primer_fwdrev=="fwd" and rev_primer_fwdrev=="rev")
                        amplicon = transcriptome_ref[fwd_target_isoform_ID][fwd_pos_5p-1:rev_pos_5p]
                        try:
                            assert (amplicon.startswith(fwd_primer.seq)), "Amplicon does not begin with forward primer"
                        except AssertionError as ae:
                            #print(fwd_primer_fwdrev, fwd_match_fwdrev, fwd_target_isoform_ID, fwd_pos_5p, fwd_frac_duplexed, file=sys.stderr, flush=True)
                            #print(amplicon, file=sys.stderr, flush=True)
                            print(ae, file=sys.stderr, flush=True)
                            #pdb.set_trace()
                            
                        rev_primer_compl = rev_primer.seq[::-1].translate(DNA_complement_table)
                        try:
                            assert (amplicon.endswith(rev_primer_compl)), "Amplicon does not end with the reverse primer"
                        except AssertionError as ae:
                            #print(rev_primer_fwdrev, rev_match_fwdrev, rev_target_isoform_ID, rev_pos_5p, rev_frac_duplexed, file=sys.stderr, flush=True)
                            #print(amplicon, file=sys.stderr, flush=True)
                            print(ae, file=sys.stderr, flush=True)
                            #pdb.set_trace()

                        read_pair = (amplicon[0:read_len], amplicon[-read_len:])
                        if (fwd_target_isoform_ID in isoforms_grouped_by_paired_reads[read_pair]):
                            pdb.set_trace()
                        isoforms_grouped_by_paired_reads[read_pair].append( fwd_target_isoform_ID )
                        amplicons_per_paired_reads[read_pair].add(amplicon)

        expected_global_tigs = set(sop.target_isoform_groups_by_ID)
        expected_nontarget_igs = set(sop.nontarget_isoform_groups_by_ID) if (sop.nontarget_isoform_groups_by_ID != None) else set()
        other_nontarget_igs = set(sop.other_isoform_groups_by_ID) if (sop.other_isoform_groups_by_ID != None) else set()
        all_nontarget_igs = expected_nontarget_igs | other_nontarget_igs

        primed_ig_tuples_by_ID = set()
        amplicons_per_ig_tuple = {}
        for read_pair, isoform_IDs in isoforms_grouped_by_paired_reads.items():
            isoform_IDs.sort()
            ig_tuple = tuple(isoform_IDs)
            primed_ig_tuples_by_ID.add(ig_tuple)
            amplicons_per_ig_tuple[ig_tuple] = amplicons_per_paired_reads[read_pair]

        try:
            assert (expected_global_tigs <= primed_ig_tuples_by_ID) # | expected_offtarget_igs 
            #print("\tHas expected global TIGs", file=sys.stderr, flush=True)
        except AssertionError:
            pdb.set_trace()
            print(primed_ig_tuples_by_ID, file=sys.stderr, flush=True)

        tig_amplicon = {}
        for global_tig in expected_global_tigs:
            if (len(amplicons_per_ig_tuple[global_tig]) == 1):
                tig_amplicon[global_tig] = list(amplicons_per_ig_tuple[global_tig])[0]
            else:
                try:
                    assert (len(set(map(lambda a: (a[0:read_len],a[-read_len:]), amplicons_per_ig_tuple[global_tig]))) == 1)
                except AssertionError as ae:
                    pdb.set_trace()
                tig_amplicon[global_tig] = list(amplicons_per_ig_tuple[global_tig])

        offtarget_ig_amplicons = {}
        offtarget_igs = primed_ig_tuples_by_ID - expected_global_tigs
        #print("\tExpected offtarget IGs: ", expected_offtarget_igs, file=sys.stderr, flush=True)
        #print("\tOfftarget IGs: ", offtarget_igs, file=sys.stderr, flush=True)
        for offtarget_ig in offtarget_igs:
            if (len(amplicons_per_ig_tuple[offtarget_ig]) == 1):
                offtarget_ig_amplicons[offtarget_ig] = list(amplicons_per_ig_tuple[offtarget_ig])[0]
            else:
                offtarget_ig_amplicons[offtarget_ig] = list(amplicons_per_ig_tuple[offtarget_ig])

    return (tig_amplicon, offtarget_ig_amplicons)


def confirmIsoformGroupsCorrectlyPrimed_Global(clique_sop_descriptions, all_sets_of_primers, results_per_primer_seq, transcriptome_fasta, all_target_isoform_IDs):
    read_len = designParams.read_len
    transcriptome_ref = Fasta(transcriptome_fasta, as_raw=True, sequence_always_upper=True)
    all_fwd_primers = []
    all_rev_primers = []
    
    for sop_descr in clique_sop_descriptions:
        sop = all_sets_of_primers[sop_descr]
        #print(sop_descr, file=sys.stderr, flush=True)
        all_annot_primers = sop.getAllAnnotPrimers() 
        all_annot_primers_aug = list(map(lambda x: x + (sop_descr,), all_annot_primers))
        all_fwd_primers += list(filter(lambda t: t[0][0]=='F', all_annot_primers_aug))
        all_rev_primers += list(filter(lambda t: t[0][0]=='R', all_annot_primers_aug))

    isoforms_grouped_by_paired_reads = defaultdict(list)
    amplicons_per_paired_reads = defaultdict(set)
    amplicons_per_isoform = defaultdict(set)
    print("Accumulating amplicons per isoform", file=sys.stderr)
    for (fwd_erID, fwd_primer, fwd_sop_descr), (rev_erID, rev_primer, rev_sop_descr) in product(all_fwd_primers, all_rev_primers):
        for (fwd_primer_fwdrev, fwd_match_fwdrev, fwd_isoform_ID, fwd_pos_5p, fwd_frac_duplexed) in results_per_primer_seq[fwd_primer.seq]:
            for (rev_primer_fwdrev, rev_match_fwdrev, rev_isoform_ID, rev_pos_5p, rev_frac_duplexed) in results_per_primer_seq[rev_primer.seq]:

                if (fwd_isoform_ID==rev_isoform_ID and fwd_match_fwdrev=="fwd" and rev_match_fwdrev=="rev" and
                    rev_pos_5p - fwd_pos_5p > 10 and rev_pos_5p - fwd_pos_5p < designParams.confusion_amplicon_max_len):
                    tup = (fwd_pos_5p, fwd_primer_fwdrev, fwd_match_fwdrev, fwd_frac_duplexed, rev_pos_5p, rev_primer_fwdrev, rev_match_fwdrev, rev_frac_duplexed)
                    amplicons_per_isoform[fwd_isoform_ID].add(tup)
                    
                    assert (fwd_primer_fwdrev=="fwd" and rev_primer_fwdrev=="rev")
                    amplicon = transcriptome_ref[fwd_isoform_ID][fwd_pos_5p-1:rev_pos_5p]
                    try:
                        assert (amplicon.startswith(fwd_primer.seq)), "Amplicon does not begin with forward primer"
                    except AssertionError as ae:
                        #print(fwd_primer_fwdrev, fwd_match_fwdrev, fwd_isoform_ID, fwd_pos_5p, fwd_frac_duplexed, file=sys.stderr, flush=True)
                        #print(amplicon, file=sys.stderr, flush=True)
                        #print(ae, file=sys.stderr, flush=True)
                        #pdb.set_trace()
                        pass
                    
                    rev_primer_compl = rev_primer.seq[::-1].translate(DNA_complement_table)
                    try:
                        assert (amplicon.endswith(rev_primer_compl)), "Amplicon does not end with the reverse primer"
                    except AssertionError as ae:
                        #print(rev_primer_fwdrev, rev_match_fwdrev, rev_isoform_ID, rev_pos_5p, rev_frac_duplexed, file=sys.stderr, flush=True)
                        #print(amplicon, file=sys.stderr, flush=True)
                        #print(ae, file=sys.stderr, flush=True)
                        #pdb.set_trace()
                        pass
                    
                    read_pair = (amplicon[0:read_len], amplicon[-read_len:])
                    isoforms_grouped_by_paired_reads[read_pair].append( fwd_isoform_ID )
                    amplicons_per_paired_reads[read_pair].add(amplicon)


    for isoform_ID, amplicons_5p3p in amplicons_per_isoform.items():
        target_descriptor = "target" if (isoform_ID in all_target_isoform_IDs) else "non-target"
        print("Checking for amplicon conflicsts on %s %s" % (target_descriptor, isoform_ID), file=sys.stderr)
        truth_val = all(map(lambda xy: xy[1][4]<xy[0][0] or xy[1][0]>xy[0][4], combinations(amplicons_5p3p,2)))
        if (not truth_val):
            pdb.set_trace()

    pdb.set_trace()

    expected_global_tigs = set(sop.target_isoform_groups_by_ID)
    expected_nontarget_igs = set(sop.nontarget_isoform_groups_by_ID) if (sop.nontarget_isoform_groups_by_ID != None) else set()
    other_nontarget_igs = set(sop.other_isoform_groups_by_ID) if (sop.other_isoform_groups_by_ID != None) else set()
    all_nontarget_igs = expected_nontarget_igs | other_nontarget_igs

    primed_ig_tuples_by_ID = set()
    amplicons_per_ig_tuple = {}
    for read_pair, isoform_IDs in isoforms_grouped_by_paired_reads.items():
        isoform_IDs.sort()
        ig_tuple = tuple(isoform_IDs)
        primed_ig_tuples_by_ID.add(ig_tuple)
        amplicons_per_ig_tuple[ig_tuple] = amplicons_per_paired_reads[read_pair]

    try:
        assert (expected_global_tigs <= primed_ig_tuples_by_ID) # | expected_offtarget_igs 
        #print("\tHas expected global TIGs", file=sys.stderr, flush=True)
    except AssertionError:
        pdb.set_trace()
        print(primed_ig_tuples_by_ID, file=sys.stderr, flush=True)

    tig_amplicon = {}
    for global_tig in expected_global_tigs:
        if (len(amplicons_per_ig_tuple[global_tig]) == 1):
            tig_amplicon[global_tig] = list(amplicons_per_ig_tuple[global_tig])[0]
        else:
            try:
                assert (len(set(map(lambda a: (a[0:read_len],a[-read_len:]), amplicons_per_ig_tuple[global_tig]))) == 1)
            except AssertionError as ae:
                pdb.set_trace()
            tig_amplicon[global_tig] = list(amplicons_per_ig_tuple[global_tig])

    offtarget_ig_amplicons = {}
    offtarget_igs = primed_ig_tuples_by_ID - expected_global_tigs
    #print("\tExpected offtarget IGs: ", expected_offtarget_igs, file=sys.stderr, flush=True)
    #print("\tOfftarget IGs: ", offtarget_igs, file=sys.stderr, flush=True)
    for offtarget_ig in offtarget_igs:
        if (len(amplicons_per_ig_tuple[offtarget_ig]) == 1):
            offtarget_ig_amplicons[offtarget_ig] = list(amplicons_per_ig_tuple[offtarget_ig])[0]
        else:
            offtarget_ig_amplicons[offtarget_ig] = list(amplicons_per_ig_tuple[offtarget_ig])

    return (tig_amplicon, offtarget_ig_amplicons)


def writePool(clique_sop_descriptions, all_sets_of_primers, pool_type, plate_basename, output_tsv):
    pool_fwd_primers = [[], [], [], []]
    pool_rev_primers = [[], [], [], []]

    sizes_of_ppr_combos_used = []
    for sop_descr in clique_sop_descriptions:
        chromosome, olap_set_ID, ppr_combo_IDs_descr, pp_set_num = sop_descr

        sop = all_sets_of_primers[sop_descr]
        all_annot_primers = sop.getAllAnnotPrimers() 
        fwd_primers = []
        rev_primers = []
        for fwd_primer in filter(lambda t: t[0][0]=='F', all_annot_primers):
            seq_name = "%s_%s_%d_%s" % ((chromosome, olap_set_ID, pp_set_num, fwd_primer[0]))
            fwd_primers.append( (seq_name, ppr_combo_IDs_descr, fwd_primer[1]) )
        for rev_primer in filter(lambda t: t[0][0]=='R', all_annot_primers):
            seq_name = "%s_%s_%d_%s" % ((chromosome, olap_set_ID, pp_set_num, rev_primer[0]))
            rev_primers.append( (seq_name, ppr_combo_IDs_descr, rev_primer[1]) )

        num_fwd = len(fwd_primers)
        num_rev = len(rev_primers)

        for pool_ind in [0,1,2,3]:
            fwd_remain = 96 - len(pool_fwd_primers[pool_ind])
            rev_remain = 96 - len(pool_rev_primers[pool_ind])
            if (fwd_remain >= num_fwd and rev_remain >= num_rev):
                sizes_of_ppr_combos_used.append( ppr_combo_IDs_descr.count('+')+1 )
                pool_fwd_primers[pool_ind].extend(fwd_primers)
                pool_rev_primers[pool_ind].extend(rev_primers)
                break

    plate_cols = list(range(1,13))
    plate_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    well_loading_order = list(map(lambda cr: "%s%d" % (cr[1],cr[0]), product(*[plate_cols,plate_rows]))) # col,row -> row,col
    scale = "500 picomole DNA Plate Oligo"
    purification = "Standard Desalting"
    norm_style = "Concentration And Quantity"
    quantity = "0.5"
    concentration = "10"
    buff = "IDTE Buffer pH 7.5 (10mH Tris-HCl, 0.1mM EDTA)"

    table_rows = []
    for pool_ind in [0,1,2,3]:
        plate_name = "%s_%s_%s" % (plate_basename, pool_ind+1, "fwd")
        fwd_primers = pool_fwd_primers[pool_ind]
        for i,well_position in enumerate(well_loading_order):
            if (i<len(fwd_primers)):
                seq_name = fwd_primers[i][0]
                seq = fwd_primers[i][2].seq
                tagged_seq = "%s%s%s" % (designParams.fwd_primer_adapter_tag, "NNNNNN", seq)
                note = fwd_primers[i][1]
                table_rows.append( (plate_name, well_position, seq_name, tagged_seq, scale, purification, norm_style, quantity, concentration, buff, note) )
            else:
                table_rows.append( (plate_name, well_position, "EMPTY", "EMPTY", "", "", "", "", "", "", "") )

        plate_name = "%s_%s_%s" % (plate_basename, pool_ind+1, "rev")
        rev_primers = pool_rev_primers[pool_ind]
        for i,well_position in enumerate(well_loading_order):
            if (i<len(rev_primers)):
                seq_name = rev_primers[i][0]
                seq = rev_primers[i][2].seq
                tagged_seq = "%s%s%s" % (designParams.rev_primer_adapter_tag, "NNNNNN", seq)
                note = rev_primers[i][1]
                table_rows.append( (plate_name, well_position, seq_name, tagged_seq, scale, purification, norm_style, quantity, concentration, buff, note) )
            else:
                table_rows.append( (plate_name, well_position, "EMPTY", "EMPTY", "", "", "", "", "", "", "") )


    with open(output_tsv, 'wt') as op:
        for row_tuple in table_rows:
            line = "\t".join(row_tuple)
            op.write("%s\n" % line)

    C = Counter(sizes_of_ppr_combos_used)
    pdb.set_trace()
    print(C.most_common(), file=sys.stderr)
    

if (__name__ == "__main__"):
    database_path, params_descriptor, target_chr_OS, tempdir_root, transcriptome_fasta, transcriptome_vsearch_udb, pool_type, plate_basename, output_tsv = sys.argv[1:]

    designParams.setParameters(params_descriptor)
    oligo_thermo = OligoThermodynamics()

    tempdir = "%s/designSets_%d" % (tempdir_root, os.getpid())
    os.mkdir(tempdir)
    print("Temporary directory is %s" % tempdir, file=sys.stderr, flush=True)

    all_sets_of_primers, primers_to_align, primer_ids_per_primer_seq, all_target_isoform_IDs, loci_isoform_IDs_per_seq = collectSetsOfPrimers(database_path, params_descriptor, target_chr_OS)

    print("INFO: read %d primers. Aligning against transcriptome..." % len(primers_to_align), file=sys.stderr, end='', flush=True)
    results_per_primer_seq, primers_w_no_perfect_match, do_not_use = alignPrimersAgainstTranscriptome(oligo_thermo, primers_to_align, transcriptome_vsearch_udb,
                                                                                                      tempdir, all_target_isoform_IDs, loci_isoform_IDs_per_seq)
    print("done", file=sys.stderr, flush=True)
    print("INFO: %d primers without perfect transcriptome match" % len(primers_w_no_perfect_match), file=sys.stderr, flush=True)
    #pickle.dump((results_per_primer_seq,primers_w_no_perfect_match, do_not_use), open('%s/results_per_primer_seq.pkl' % tempdir,'wb'))
    #results_per_primer_seq, primers_w_no_perfect_match, do_not_use = pickle.load(open('%s/results_per_primer_seq.pkl' % tempdir,'rb'))
    
    print("INFO: culling primer sets", file=sys.stderr, flush=True)
    cullPrimerSets(all_sets_of_primers, do_not_use)

    print("INFO: flagging positionally conflicting primers", file=sys.stderr, flush=True)
    incompatible_primers_by_pos = flagPositionallyConflictingPrimers(results_per_primer_seq, primer_ids_per_primer_seq, all_target_isoform_IDs)
    print("INFO: flagging thermodynamically conflicting primers", file=sys.stderr, flush=True)
    thermo_interacting_primers = flagThermodynamicallyConflictingPrimers(oligo_thermo, primers_to_align)
    #pickle.dump((incompatible_primers_by_pos, thermo_interacting_primers), open('%s/conflicted_primers.pkl' % tempdir,'wb'))
    #incompatible_primers_by_pos, thermo_interacting_primers = pickle.load(open('%s/conflicted_primers.pkl' % tempdir,'rb'))

    print("Making primer set compatibility graph", file=sys.stderr, flush=True)
    G = makeCompatibilityGraphForSetsOfPrimers(all_sets_of_primers, incompatible_primers_by_pos, thermo_interacting_primers)
    #pickle.dump(G, open('%s/G.pkl' % tempdir,'wb'))
    #G = pickle.load(open('%s/G.pkl' % tempdir,'rb'))

    print("Forming primer pool", file=sys.stderr, flush=True)
    clique_sop_descriptions = formPoolIteratively(G, all_sets_of_primers, pool_type, tempdir)
    #pickle.dump(clique_sop_descriptions, open('%s/clique_sop_descriptions.pkl' % tempdir,'wb'))
    #clique_sop_descriptions = pickle.load(open('%s/clique_sop_descriptions.pkl' % tempdir,'rb'))

    confirmIsoformGroupsCorrectlyPrimed_Local(clique_sop_descriptions, all_sets_of_primers, results_per_primer_seq, transcriptome_fasta)
    confirmIsoformGroupsCorrectlyPrimed_Global(clique_sop_descriptions, all_sets_of_primers, results_per_primer_seq, transcriptome_fasta, all_target_isoform_IDs)
    
    writePool(clique_sop_descriptions, all_sets_of_primers, pool_type, plate_basename, output_tsv)
    os.rmdir(tempdir)

    sys.exit(0)
