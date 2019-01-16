#!/usr/bin/env python
from __future__ import print_function

# NOTES
# Once graph of node=primer pair and edges=primer pairs not from same pprID, 
# Filter nodes based on
#   - either fwd or rev primer having bad self-interaction (homodimer or hairpin)
#   - fwd/rev pair having bad heterodimer interaction (ie frac_heteroduplex)
#   - either fwd or rev having high or many hybridization penalties
#
# Filter edges based on
#   - rev primers or fwd primers having bad heterodimer interactions
#   - either fwd-rev cross combination having high or many product penalties (from unwanted amplifications)
#
# To enable filtering:
# For nodes record
#   - fwd_has_hairpin, fwd_can_self_prime, fwd_frac_homoduplex
#   - rev_has_hairpin, rev_can_self_prime, rev_frac_homoduplex
#   - fwdrev_can_heteroprime, fwdrev_frac_heteroduplex
#   - fwd_hyb_scores     (fraction duplexed to intended template sites)
#   - rev_hyb_scores
#   - fwd_hyb_penalties  (fraction duplexed to unintended template sites)
#   - rev_hyb_penalties
#   - TODO: INCLUDE PRODUCT SCORES?
#
# For edges record
#   - fwd1_rev2_can_heteroprime, fwd1_rev2_frac_heteroduplex
#   - fwd2_rev1_can_heteroprime, fwd1_rev2_frac_heteroduplex
#   - fwd1_rev2_product_penalties
#   - fwd2_rev1_product_penalties
#
#
# Filtering strategy: Begin strict, and until pool can be formed, decrease strictness iteratively
#   -- NODES --
#   1) Nodes with fwd or rev hairpin always filtered
#   2) Same for each fwd & rev primer: [(consider_self_prime, max_frac_homoduplex), ...] with consider_self_prime=False for low max_frac_homoduplex values
#   3) For the [(consider_hetero_prime, max_frac_heteroduplex), ...] with consider_self_prime=False for low max_frac_homoduplex values.
#      This filtering should be more lenient, as rev primer concencration should be very low compared to fwd primer concentration. Maybe this should be
#      handled upstream in the heteroduplex calculation, with the primer concentration adjusted there. Then 2) & 3) could be handled with same logic.
#   4) For rev primer, calculate sum ontarget frac_hyb_duplex. Then calc sum offtarget frac_hyb_duplex. Then calc ratio off/on frac_hyb_duplex and filter
#      based on max allowed off/on ratio: [0.001, 0.005, 0.01, 0.05, 0.1]
#   5) For fwd primer, calculate sum ontarget frac_hyb_duplex. Then calc sum offtarget frac_hyb_duplex for cases where the corresponding rev primer wasn't
#      filtered out in 4). Apply same off/on ratio filtering.
#
#   -- EDGES --
#   6) For fwd1/rev2 & fwd2/rev1 evals, apply same logic as in 3) above.
#   7) For fwd/rev product penalties, only consider those fwd/rev for which the fwd & rev primers have not been filtered out.
#      Filter on max fwd_frac_duplex*rev_frac_duplex: [0.00001, 0.0001, 0.001, 0.01]
#
#   -- REFINEMENTS? --
#   o For 4), 5), 7), add extra "stringency factor" if isoform is rRNA.
#   o Enforce including all or none of the primers in a chr_OSID_part experiment partition.
#   o Include Universal adapters, as these could cause "1-sided" products. ADDED. NEED TO CONFIRM.
#   o Individual rev and fwd primers matching (as they should) to other isoforms in the same locus should not be penalized.
#   o Screening against transcriptome+intron needs to be able to distinguish new primer hybs and new products from those that are the same
#     as were found in just transcriptome
#   o CTRL primers don't have expected primer products because they are not in the predesign database. Thus all correctly found products are penalties,
#     and they should not be.
#   o Address "WARNING: product start/stop positions didn't match up"

from itertools import product
from operator import itemgetter, mul, xor
from collections import defaultdict
from math import exp
from subprocess import Popen, PIPE, check_call, check_output, CalledProcessError
import cPickle
import sqlite3
import sys
import os
import glob
import tempfile
import primer3
import networkx as nx
import numpy as np

import pyfaidx
from pyfaidx import Fasta

sys.path.append("/raid1/projects/CGDB/scripts")
from classMFE import MFE, MFEOptions

import designParams

import pdb

from string import maketrans, translate
DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")



def readLociFile(ribosomal_loci):
    with open(ribosomal_loci, 'r') as ip:
        loci = filter(lambda x: x != '', map(lambda x: x.strip(), ip.read().split("\n")))

    return set(loci)

        
def readPrimersFasta(primers_fasta):
    fwd_primers = {}
    rev_primers = {}
    fwd_rev_primers = {}
    
    with open(primers_fasta, 'r') as ip:
        id_and_seq_pairs = filter(lambda x: x != '', ip.read().split("\n"))

    for seq_ID, seq in zip(id_and_seq_pairs[0::2], id_and_seq_pairs[1::2]):
        primer_seq = seq.split('N')[-1] if ('N' in seq) else seq
        if (seq_ID.endswith("_5p")):
            ID = seq_ID[1:-3]
            assert (not ID in fwd_primers)
            fwd_primers[(ID, "na")] = primer_seq
        elif (seq_ID.endswith("_3p")):
            ID = seq_ID[1:-3]
            assert (not ID in rev_primers)
            rev_primers[(ID, "na")] = primer_seq
        else:
            print(seq_ID, "doesn't end with '_5p' or '_3p'", file=sys.stderr)
            sys.exit(1)

    assert (set(fwd_primers.keys()) == set(rev_primers.keys())), "Source difference for fwd and rev primers"
        
    for key in fwd_primers.keys():
        fwd_rev_primers[key] = (fwd_primers[key], rev_primers[key])

    return fwd_rev_primers


def readAndPairPrimers(fwd_primers_tsv, rev_primers_tsv, pool_name):
    fwd_primers = {}
    rev_primers = {}
    fwd_rev_primers = {}
    
    ip = open(fwd_primers_tsv, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        primer_seq = fields[2].split('N')[-1] # TODO: split on designParams.UMI
        for full_pprID in set(fields[3].strip().split(';')):
            assert (not fwd_primers.has_key(full_pprID))
            fwd_primers[(full_pprID, pool_name)] = primer_seq
    ip.close()

    ip = open(rev_primers_tsv, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        primer_seq = fields[2].split('N')[-1] # TODO: split on designParams.UMI
        for full_pprID in set(fields[3].strip().split(';')):
            assert (not rev_primers.has_key(full_pprID))
            rev_primers[(full_pprID, pool_name)] = primer_seq
    ip.close()

    assert (set(fwd_primers.keys()) == set(rev_primers.keys())), "Source difference for fwd and rev primers"
        
    for key in fwd_primers.keys():
        fwd_rev_primers[key] = (fwd_primers[key], rev_primers[key])

    return fwd_rev_primers


def getPrimersKnownTargets(fwd_rev_primers, database_path, transcriptomic_MFE):
    expected_primed_isoforms = defaultdict(set)
    expected_primer_products = defaultdict(list)
    
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    for (full_pprID, pool_name), (fwd_primer_seq, rev_primer_seq) in fwd_rev_primers.items():
        try:
            chrom, OS_ID, experiment_num, PPR_ID = full_pprID.split('_')
        except ValueError as ve:
            if (full_pprID.startswith("ERCC") or full_pprID.endswith("CTRL") or full_pprID.startswith("UDTD")):
                continue
            else:
                raise ve
                
        select_tup = (chrom, OS_ID, params_descriptor, int(experiment_num), PPR_ID, fwd_primer_seq, rev_primer_seq)
        
        cursor.execute("SELECT global_target_isoform_group,offtarget_isoform_groups from PrimerPair \
        WHERE chrom=? AND OS_ID=? AND params_descriptor=? AND part_num=? AND PPR_ID=? and fwd_primer_seq=? AND rev_primer_seq=?", select_tup)

        results = cursor.fetchall()
        assert (len(results) == 1)
        global_TIGs = results[0][0].split(';')
        offtarget_IGs = results[0][1].split(';')
            
        for ig in filter(lambda x: x != '', global_TIGs + offtarget_IGs):
            isoforms = ig.split(',')
            expected_primed_isoforms[full_pprID].update(isoforms)

        ampl_success = transcriptomic_MFE.processPrimer(fwd_primer_seq, rev_primer_seq)
        assert (ampl_success), "Primer amplification failure"
        products_per_isoform = transcriptomic_MFE.getPrimerProducts(designParams.is_directional, 1000)
        assert (len(products_per_isoform)>0), "Amplification, but no products."

        isoforms_per_product = defaultdict(list)
        for isoform_ID, products in products_per_isoform.items():
            for product in products:
                # Note: local_start and local_stop are 0-based coordinates
                local_start, local_stop, nuc_seq, is_fwd_rev_primer_product = products[0][0:4]
                if (is_fwd_rev_primer_product):
                    try:
                        assert (any(map(lambda tig: isoform_ID in tig, global_TIGs)) or any(map(lambda ig: isoform_ID in ig, offtarget_IGs))
                                or len(nuc_seq) > designParams.confusion_amplicon_max_len), \
                                "Primers yield product from isoform %s, which is not in primer pair's global TIG or offtarget IG "
                    except AssertionError as ae:
                        pdb.set_trace()
                    
                    expected_primer_products[(isoform_ID, full_pprID)].append( (local_start, local_stop) )

    conn.close()

    return (expected_primed_isoforms, expected_primer_products)


def compareRTPrimersPairwise(primers1, primers2, rxn_temp_C=57, primer_concentration=1.0*10**-8):
    '''primer_concentration is in M units, and rxn_temp_C is in Celsius'''
    primer_self_interaction = {}
    primer_pair_interaction = {}

    # Reverse complements are used to identify complementarity by identifying identical nucleotides (ie identify A/T complementary bp by A/A)
    primers1_aug = map(lambda p: (p[0], p[1], np.array(list(p[1])), np.array(list(p[1][::-1].translate(DNA_complement_table)))), primers1)
    primers2_aug = map(lambda p: (p[0], p[1], np.array(list(p[1])), np.array(list(p[1][::-1].translate(DNA_complement_table)))), primers2)

    thermo_calc = primer3.thermoanalysis.ThermoAnalysis(mv_conc=50, dv_conc=1.5, dntp_conc=0.25, temp_c=rxn_temp_C, dna_conc=50)

    rxn_temp_K = 273.15 + rxn_temp_C
    RT = 1.98717 * rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)
    # dG units are expected to be cal/mol and NOT kcal/mol
    calc_frac_duplexed = lambda dG: (primer_concentration * exp(-dG/RT))/(1 + primer_concentration * exp(-dG/RT))

    # Compile individual primer data
    for (seq_id, seq, aseq, aseq_rc) in primers1_aug + primers2_aug:
        if (not primer_self_interaction.has_key(seq)):
            seq_hairpin = thermo_calc.calcHairpin(seq)
            seq_3p_can_prime_on_self = any(map(lambda s: sum(s != aseq_rc[0:6]) <= 1, map(lambda i: aseq[i:i+6], range(len(seq)-5))))
            any_hyb = thermo_calc.calcHomodimer(seq)
            frac_duplexed_any = calc_frac_duplexed(any_hyb.dg)
            primer_self_interaction[seq] = (seq_hairpin.structure_found, seq_3p_can_prime_on_self, frac_duplexed_any)

    # Compile primer pair data
    for (seq_id1, seq1, aseq1, aseq1_rc), (seq_id2, seq2, aseq2, aseq2_rc) in product(primers1_aug, primers2_aug):
        if (primer_pair_interaction.has_key( (seq2, seq1) )):
            # Efficiency for case where primers1 is the same as primers2
            primer_pair_interaction[(seq1, seq2)] = primer_pair_interaction[(seq2, seq1)]
        elif (seq1 != seq2):
            # Check for one priming off the other. The values chosen here: RT can prime from hexamers, and fairly robustly even with a 1bp mismatch
            seq1_3p_can_prime_on_seq2 = any(map(lambda s2: sum(s2 != aseq1_rc[0:6]) <= 1, map(lambda i: aseq2[i:i+6], range(len(seq2)-5))))
            seq2_3p_can_prime_on_seq1 = any(map(lambda s1: sum(s1 != aseq2_rc[0:6]) <= 1, map(lambda i: aseq1[i:i+6], range(len(seq1)-5))))

            any_hyb = thermo_calc.calcHeterodimer(seq1, seq2)
            end1_hyb = thermo_calc.calcEndStability(seq1, seq2)
            end2_hyb = thermo_calc.calcEndStability(seq2, seq1)
            min_dG = min(any_hyb.dg, end1_hyb.dg, end2_hyb.dg)
            frac_duplexed_any = calc_frac_duplexed(min_dG)

            primer_pair_interaction[(seq1, seq2)] = (seq1_3p_can_prime_on_seq2 or seq2_3p_can_prime_on_seq1, frac_duplexed_any)

    return (primer_self_interaction, primer_pair_interaction)


def compareDNAPolPrimersPairwise(primers1, primers2, rxn_temp_C=57, primer_concentration=1.0*10**-8):
    '''primer_concentration is in M units, and rxn_temp_C is in Celsius'''
    primer_self_interaction = {}
    primer_pair_interaction = {}

    # Reverse complements are used to identify complementarity by identifying identical nucleotides (ie identify A/T complementary bp by A/A)
    primers1_aug = map(lambda p: (p[0], p[1], np.array(list(p[1])), np.array(list(p[1][::-1].translate(DNA_complement_table)))), primers1)
    primers2_aug = map(lambda p: (p[0], p[1], np.array(list(p[1])), np.array(list(p[1][::-1].translate(DNA_complement_table)))), primers2)
        
    thermo_calc = primer3.thermoanalysis.ThermoAnalysis(mv_conc=50, dv_conc=1.5, dntp_conc=0.25, temp_c=rxn_temp_C, dna_conc=50)

    rxn_temp_K = 273.15 + rxn_temp_C
    RT = 1.98717 * rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)
    calc_frac_duplexed = lambda dG: (primer_concentration * exp(-dG/RT))/(1 + primer_concentration * exp(-dG/RT))

    # Compile primer self data
    for (seq_id, seq, aseq, aseq_rc) in primers1_aug + primers2_aug:
        if (not primer_self_interaction.has_key(seq)):
            seq_hairpin = thermo_calc.calcHairpin(seq)
            seq_3p_can_prime_on_self = any(map(lambda s: sum(s != aseq_rc[0:6]) <= 1, map(lambda i: aseq[i:i+6], range(len(seq)-5))))
            any_hyb = thermo_calc.calcHomodimer(seq)
            frac_duplexed_any = calc_frac_duplexed(any_hyb.dg)
            primer_self_interaction[seq] = (seq_hairpin.structure_found, seq_3p_can_prime_on_self, frac_duplexed_any)

    # CAREFUL: The 2nd nt in each tuple is the reverse complement of what is shown in Table 1 of 
    # "The Effect of Primer-Template Mismatches on the Detection and Quantificatin of Nucleic Acids Using the 5' Nuclease Assay", Stadhouders et al
    # because of how "complementarity" is being identified below.
    #bad_terminal_mm = set([('A','T'), ('A','C'), ('G','T'), ('G','C'), ('C','G'), ('T','A'), ('T','G'),('C','A')])
    #bad_penult_mm = set([('A','T'), ('A','C'), ('G','T'), ('G','C'), ('C','G')]) 

    for (seq_id1, seq1, aseq1, aseq1_rc), (seq_id2, seq2, aseq2, aseq2_rc) in product(primers1_aug, primers2_aug):
        if (primer_pair_interaction.has_key( (seq2, seq1) )):
            # Efficiency for case where primers1 is the same as primer2
            primer_pair_interaction[(seq1, seq2)] = primer_pair_interaction[(seq2, seq1)]

        elif (seq1 != seq2):
            # Check for one priming off the other using a primer length of 6bp (to be conservative, given that Taq needs at least 8bp primer to function. See
            # "Polymerization behavior of Klenow fragment and Taq polymerase in short primer extension reactions, by Zhao)
            # Can prime if exact heptamer found or heptamer with 1bp mismatch that isn't one of "bad" mismatch cases in 

            seq1_3p_can_prime_on_seq2 = any(map(lambda s2: sum(s2 != aseq1_rc[0:6]) <= 1, map(lambda i: aseq2[i:i+6], range(len(seq2)-6))))
            seq2_3p_can_prime_on_seq1 = any(map(lambda s1: sum(s1 != aseq2_rc[0:6]) <= 1, map(lambda i: aseq1[i:i+6], range(len(seq1)-6))))

            #seq1_3p_can_prime_on_seq2 = any(map(lambda s2: sum(s2 != aseq1_rc[0:6]) <= 1 and (s2[0],aseq1_rc[0]) not in bad_terminal_mm and
            #                                    (s2[1],aseq1_rc[1]) not in bad_penult_mm, map(lambda i: aseq2[i:i+6], range(len(seq2)-6))))

            #seq2_3p_can_prime_on_seq1 = any(map(lambda s1: sum(s1 != aseq2_rc[0:6]) <= 1 and (s1[0],aseq2_rc[0]) not in bad_terminal_mm and
            #                                    (s1[1],aseq2_rc[1]) not in bad_penult_mm, map(lambda i: aseq1[i:i+6], range(len(seq1)-6))))

            # Debug 
            if (False): # and seq1_3p_can_prime_on_seq2 or seq2_3p_can_prime_on_seq1):
                print(seq1, seq2, sep="\t", file=sys.stderr)
                for s2 in map(lambda i: aseq2[i:i+6], range(len(seq2)-6)):
                    print(s2, aseq1_rc[0:6], sum(s2 != aseq1_rc[0:6]) <= 1, (s2[0],aseq1_rc[0]) not in bad_terminal_mm, (s2[1],aseq1_rc[1]) not in bad_penult_mm, file=sys.stderr)
                print("--", file=sys.stderr)
                for s1 in map(lambda i: aseq1[i:i+6], range(len(seq1)-6)):
                    print(s1, aseq2_rc[0:6], sum(s1 != aseq2_rc[0:6]) <= 1, (s1[0],aseq2_rc[0]) not in bad_terminal_mm, (s1[1],aseq2_rc[1]) not in bad_penult_mm, file=sys.stderr)
                print("--", file=sys.stderr)
                
            any_hyb = thermo_calc.calcHeterodimer(seq1, seq2)
            end1_hyb = thermo_calc.calcEndStability(seq1, seq2)
            end2_hyb = thermo_calc.calcEndStability(seq2, seq1)
            min_dG = min(any_hyb.dg, end1_hyb.dg, end2_hyb.dg)
            frac_duplexed_any = calc_frac_duplexed(min_dG)

            primer_pair_interaction[(seq1, seq2)] = (seq1_3p_can_prime_on_seq2 or seq2_3p_can_prime_on_seq1, frac_duplexed_any)

    return (primer_self_interaction, primer_pair_interaction)


def screenForUnintendedProducts(fwd_rev_primers, transcriptome_fasta, expected_primed_isoforms, expected_primer_products, temp_dir, high_stringency_loci=set()):
    '''Place 3' ends of forward and reverse primers in file and screen against transcriptome with usearch'''

    transcriptome_ref = Fasta(transcriptome_fasta, as_raw=True, sequence_always_upper=True)

    # Place the 3' ends of the fwd and rev primers in a file. Reverse complement the reverse primers so that usearch can be used in "plus" mode.
    primer_lookup = {}
    primers_file = "%s/primers.fa" % temp_dir
    primer_lookup = cPickle.load(open("%s/primer_lookup.pkl" % temp_dir, 'rb'))

    # Split the transcriptome file so that 32-bit usearch doesn't hit memory limit
    #split_cmd = ['faSplit', 'sequence', transcriptome_fasta, '25', "%s/transcriptome_split_" % temp_dir]
    #try:
    #    check_call(split_cmd)
    #except CalledProcessError as cpe:
    #    pdb.set_trace()
        
    primer_hyb_scores = defaultdict(list)
    primer_hyb_penalties = defaultdict(list)
    product_penalties = defaultdict(list)
    product_scores = defaultdict(list)

    # Match primers to templates, file-by-file
    for split_transcriptome_fasta in glob.glob("%s/transcriptome_split_*.fa" % temp_dir):
        if (not split_transcriptome_fasta.endswith("split_16.fa")):
            continue

        print("INFO: primer scanning", split_transcriptome_fasta, file=sys.stderr)
        usearch_output_file = split_transcriptome_fasta.replace(".fa", ".usearch")
        usearch_cmd = ['usearch', '-search_oligodb', split_transcriptome_fasta, '-db', primers_file, '-strand', 'plus', '-maxdiffs', '1',
                       '-threads', '6', '-target_cov', '1.0', '-userout', usearch_output_file, '-userfields', 'query+qstrand+qlo+qhi+target+tstrand+mism']
        try:
            if (not os.path.isfile(usearch_output_file)):
                check_call(usearch_cmd)
        except CalledProcessError as cpe:
            pdb.set_trace()

        ip = open(usearch_output_file, 'r')
        curr_isoform_ID = None
        curr_locus = None
        curr_isoform_matches = []
        for line in ip:
            # Example line: CG_1.1.6.1    +    469    474    chr17_OS2_1_F13R44 pool2 rev  +  1
            fields = line.strip().split('\t')

            fields[4] = tuple(fields[4].split()) # Split into full_pprID, pool_name, primer_type
            assert (fields[1] == '+')
            assert (fields[5] == '+')

            if (fields[0] != "CG_97828.2.5.7" and len(curr_isoform_matches)==0):
                continue

            if (curr_isoform_ID == None):
                curr_isoform_ID = fields[0]
                curr_locus = curr_isoform_ID.split('.',1)[0]

            if (fields[0] == curr_isoform_ID):
                curr_isoform_matches.append( fields[1:] )
            else:
                is_high_stringency_locus = curr_locus in high_stringency_loci
                calcPrimerAndProductPenalties(curr_isoform_ID, curr_isoform_matches, transcriptome_ref, primer_lookup, expected_primed_isoforms,
                                              expected_primer_products, primer_hyb_scores, primer_hyb_penalties, product_penalties, product_scores, is_high_stringency_locus)
                curr_isoform_ID = fields[0]
                curr_locus = curr_isoform_ID.split('.',1)[0]
                curr_isoform_matches = [ fields[1:] ]
        ip.close()

        # Last one
        is_high_stringency_locus = curr_locus in high_stringency_loci
        calcPrimerAndProductPenalties(curr_isoform_ID, curr_isoform_matches, transcriptome_ref, primer_lookup, expected_primed_isoforms,
                                      expected_primer_products, primer_hyb_scores, primer_hyb_penalties, product_penalties, product_scores, is_high_stringency_locus)

    return (primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties)


def calcPrimerAndProductPenalties(isoform_ID, isoform_matches, transcriptome_ref, primer_lookup, expected_primed_isoforms, 
                                  expected_primer_products, primer_hyb_scores, primer_hyb_penalties, product_penalties, product_scores, is_high_stringency_locus,
                                  ignore_product_len=1000, rxn_temp_C=57, primer_concentration=1.0*10**-8):
    '''Calculate and record the "fraction duplexed" penalty for each primer whose 3' end can (at least weakly) hybridize to the isoform template.
       Additionally, calculate and record the penalties for primer pairs that can potentially form a product.
       Do not compute penalties for primers hybridizing to their intended sites and for primer pairs that generate intended products.
       primer_concentration is in M.
    '''
    
    thermo_calc = primer3.thermoanalysis.ThermoAnalysis(mv_conc=50, dv_conc=1.5, dntp_conc=0.25, temp_c=rxn_temp_C, dna_conc=50)
    rxn_temp_K = 273.15 + rxn_temp_C
    RT = 1.98717 * rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)
    # dG units are expected to be cal/mol and NOT kcal/mol
    calc_frac_duplexed = lambda dG: (primer_concentration * exp(-dG/RT))/(1 + primer_concentration * exp(-dG/RT))
    
    frac_duplex_cutoff = 0.0001 if (is_high_stringency_locus) else 0.001

    isoform_matches.sort(key=itemgetter(1))
    unique_isoform_matches = defaultdict(set)
    section_header = "### ! %s ! ###\n" % isoform_ID if (is_high_stringency_locus) else "### %s ###\n" % isoform_ID
    for isoform_strand, isoform_start, isoform_stop, primer_ID, primer_strand, num_mismatch in isoform_matches:
        isoform_start = int(isoform_start)
        isoform_stop = int(isoform_stop)
        full_pprID, pool_name, primer_type = primer_ID
        primer_seq = primer_lookup[primer_ID]

        # Usearch is 0-based for start & stop. Slicing coordinates for pyfaidx are 0-based
        template_start = isoform_stop-len(primer_seq) if (primer_type == "fwd") else isoform_start-1
        template_stop = isoform_start+len(primer_seq)-1 if (primer_type == "rev") else isoform_stop
        if (template_start < 0 or template_stop > transcriptome_ref[isoform_ID]):
            continue
        else:
            template_seq = str(transcriptome_ref[isoform_ID][template_start:template_stop]) 

        if (primer_type == "fwd"):
            template_seq = template_seq[::-1].translate(DNA_complement_table)
        if (len(primer_seq) != len(template_seq)):
            try:
                assert (len(primer_seq) > isoform_stop or len(primer_seq) > len(transcriptome_ref[isoform_ID]) - isoform_start)
            except AssertionError:
                pdb.set_trace()
            continue
        
        any_hyb = thermo_calc.calcHeterodimer(primer_seq, template_seq)
        end_hyb = thermo_calc.calcEndStability(primer_seq, template_seq)
        min_dG = any_hyb.dg if (any_hyb.dg < end_hyb.dg) else end_hyb.dg
        frac_duplexed_any = calc_frac_duplexed(min_dG)
        frac_duplexed_any = round(frac_duplexed_any, 6)
        
        if (frac_duplexed_any >= frac_duplex_cutoff):
            print("%s%s %s %s mm\n %s\n %s" % (section_header, full_pprID, primer_type, num_mismatch, primer_seq, template_seq), file=sys.stderr)
            section_header = ''
            print(" Fraction duplexed: %6.5f" % frac_duplexed_any, file=sys.stderr)
            unique_isoform_matches[(template_start, template_stop, primer_type, primer_seq, template_seq, frac_duplexed_any)].add(full_pprID)

            # If the primer sequence is supposed to prime off of the isoform and does so at the correct location, then frac_duplex_any is score.
            # Otherwise, the frac_duplex_any is a penalty for the primer
            if (primer_type == "fwd"):
                if (any(map(lambda x: x[0]==template_start, expected_primer_products[(isoform_ID, full_pprID)]))):
                    primer_hyb_scores[primer_seq].append( (isoform_ID, frac_duplexed_any) )
                else:
                    primer_hyb_penalties[primer_seq].append( (isoform_ID, frac_duplexed_any) )
            else:
                if (any(map(lambda x: x[1]==template_stop, expected_primer_products[(isoform_ID, full_pprID)]))):
                    primer_hyb_scores[primer_seq].append( (isoform_ID, frac_duplexed_any) )
                else:
                    primer_hyb_penalties[primer_seq].append( (isoform_ID, frac_duplexed_any) )


    for tup_5p in filter(lambda x: x[2]=="fwd", unique_isoform_matches.keys()):
        # TODO: Compute a polymerase extension penalty for forward primer with bad mismatch in 3' terminal or penultimate position.
        #       An "extension" penalty (0.1 ?) to multiply with the forward primer hyb penalty.
        for tup_3p in filter(lambda t: t[2]=="rev", unique_isoform_matches.keys()):
            this_product_len = tup_3p[1] - tup_5p[0]
            if (this_product_len >= ignore_product_len or tup_3p[0] < tup_5p[0]):
                break # The remaining "rev" primers only create even longer products
            else:
                fwdrev_seq_tup = (tup_5p[3], tup_3p[3])
                common_full_pprIDs = unique_isoform_matches[tup_5p] & unique_isoform_matches[tup_3p]
                if (any(map(lambda cid: isoform_ID in expected_primed_isoforms[cid], common_full_pprIDs))): # If the current isoform is a target of the primer pair
                    # Confirm expected_product
                    product_start_stops_in_isoforms = set()
                    for list_of_ss_tups in map(lambda cid: expected_primer_products[(isoform_ID, cid)], common_full_pprIDs):
                        product_start_stops_in_isoforms.update(list_of_ss_tups)
                        
                    if ((tup_5p[0], tup_3p[1]) in product_start_stops_in_isoforms):
                        product_scores[fwdrev_seq_tup].append( (isoform_ID, tup_5p[-1], tup_3p[-1]) )
                    else:
                        print("WARNING: product start/stop positions didn't match up", file=sys.stderr)
                else:
                    # A pair of primers is legit if they are from same overlap set+partition and have same forward primer region xor same reverse primer region
                    # Example: ('chr1_OS5_1_F42R34', 'chr1_OS5_1_F44R34') shouldn't be penalties because they are supposed to have the same rev primer.
                    if (not tup_5p[2].startswith("chr") or not tup_3p[2].startswith("chr")):
                        product_penalties[fwdrev_seq_tup].append( (isoform_ID, this_product_len, tup_5p[-1], tup_3p[-1]) )
                    else:
                        fwd_primer_OSpart, fwd_primer_pprID = tup_5p[2].rsplit('_',1)
                        rev_primer_OSpart, rev_primer_pprID = tup_3p[2].rsplit('_',1)
                        fwd_primer_F, fwd_primer_R = fwd_primer_pprID[1:].split('R')
                        rev_primer_F, rev_primer_R = rev_primer_pprID[1:].split('R')
                        if (not(fwd_primer_OSpart == rev_primer_OSpart and xor(fwd_primer_F==rev_primer_F, fwd_primer_R==rev_primer_R))):
                            product_penalties[fwdrev_seq_tup].append( (isoform_ID, this_product_len, tup_5p[-1], tup_3p[-1]) )



def buildPrimerPairsGraph(fwd_rev_primers, fwd_primer_self_interaction, fwd_primer_pair_interaction, rev_primer_self_interaction, rev_primer_pair_interaction,
                          fwdrev_primer_pair_interaction, primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties):
    G = nx.Graph()
    
    # Add the nodes. A node is a fwd/rev primer pair
    node_index = 0
    for full_pprID_and_pool, (fwd_primer_seq, rev_primer_seq) in fwd_rev_primers.items():
        node_index += 1
        full_pprID, pool_name = full_pprID_and_pool
        fwd_has_hairpin, fwd_can_self_prime, fwd_frac_homoduplex = fwd_primer_self_interaction[fwd_primer_seq]
        rev_has_hairpin, rev_can_self_prime, rev_frac_homoduplex = rev_primer_self_interaction[rev_primer_seq]
        fwdrev_can_heteroprime, fwdrev_frac_heteroduplex = fwdrev_primer_pair_interaction[(fwd_primer_seq, rev_primer_seq)]

        node_data = {"full_pprID": full_pprID, "pool_name": pool_name, "index": node_index,
                     "fwd_has_hairpin": fwd_has_hairpin, "fwd_can_self_prime": fwd_can_self_prime, "fwd_frac_homoduplex": fwd_frac_homoduplex,
                     "rev_has_hairpin": rev_has_hairpin, "rev_can_self_prime": rev_can_self_prime, "rev_frac_homoduplex": rev_frac_homoduplex,
                     "fwdrev_can_heteroprime": fwdrev_can_heteroprime, "fwdrev_frac_heteroduplex": fwdrev_frac_heteroduplex,
                     "fwd_primer_hyb_scores": primer_hyb_scores[fwd_primer_seq],
                     "rev_primer_hyb_scores": primer_hyb_scores[rev_primer_seq],
                     "fwd_primer_hyb_penalties": primer_hyb_penalties[fwd_primer_seq],
                     "rev_primer_hyb_penalties": primer_hyb_penalties[rev_primer_seq],
                     "product_scores": product_scores[(fwd_primer_seq, rev_primer_seq)]}

        if (full_pprID.startswith("ERCC") or full_pprID.endswith("CTRL") or full_pprID.startswith("UDTD")):
            node_data["locus_group"] = (full_pprID, "na", "na")
        else:
            chrom, OSID, partition_num, pprID = full_pprID.split('_')
            node_data["locus_group"] = (chrom, OSID, partition_num)

        G.add_node((fwd_primer_seq, rev_primer_seq), attr_dict=node_data)


    # Add an edge between all potentially compatible primer pairs
    all_nodes = G.nodes()
    for i, n1 in enumerate(all_nodes[0:-1]):
        n1_index = G.node[n1]["index"]
        n1_locus_group = G.node[n1]["locus_group"]
        fwd1_primer_seq, rev1_primer_seq = n1
        for n2 in all_nodes[i+1:]:
            n2_index = G.node[n2]["index"]
            n2_locus_group = G.node[n2]["locus_group"]
            if (n1_locus_group[0:2] != n2_locus_group[0:2] or n1_locus_group[2] == n2_locus_group[2]):
                (fwd2_primer_seq, rev2_primer_seq) = n2
                try: # Maybe an "if statement", for speed
                    fwd1rev2_can_heteroprime, fwd1rev2_frac_heteroduplex = fwdrev_primer_pair_interaction[(fwd1_primer_seq, rev2_primer_seq)]
                    fwd2rev1_can_heteroprime, fwd2rev1_frac_heteroduplex = fwdrev_primer_pair_interaction[(fwd2_primer_seq, rev1_primer_seq)]
                except KeyError as ke:
                    pdb.set_trace()
                    if ((G.node[n1]["full_pprID"].startswith("UDTD57") and G.node[n2]["full_pprID"].startswith("UDTD75")) or
                        (G.node[n1]["full_pprID"].startswith("UDTD75") and G.node[n2]["full_pprID"].startswith("UDTD57"))):
                        assert(fwd1_primer_seq == rev2_primer_seq and fwd2_primer_seq == rev1_primer_seq)
                        has_hairpin, fwd1rev2_can_heteroprime, fwd1rev2_frac_heteroduplex = fwd_primer_self_interaction[fwd1_primer_seq] # has_hairpin not used
                        has_hairpin, fwd2rev1_can_heteroprime, fwd2rev1_frac_heteroduplex = fwd_primer_self_interaction[fwd2_primer_seq]
                    else:
                        raise ke
                        
                edge_data = {"node_indices": (n1_index, n2_index),
                             "fwd1rev2_can_heteroprime": fwd1rev2_can_heteroprime, "fwd1rev2_frac_heteroduplex": fwd1rev2_frac_heteroduplex,
                             "fwd2rev1_can_heteroprime": fwd2rev1_can_heteroprime, "fwd2rev1_frac_heteroduplex": fwd2rev1_frac_heteroduplex,
                             "fwd1rev2_product_penalties": product_penalties[(fwd1_primer_seq, rev2_primer_seq)],
                             "fwd2rev1_product_penalties": product_penalties[(fwd2_primer_seq, rev1_primer_seq)]}

                if ((G.node[n1]["full_pprID"].startswith("UDTD") or G.node[n2]["full_pprID"].startswith("UDTD")) and
                    (fwd1rev2_can_heteroprime or fwd2rev1_can_heteroprime or
                    len(edge_data["fwd1rev2_product_penalties"])>0 or len(edge_data["fwd2rev1_product_penalties"])>0)):
                    print(G.node[n1]["full_pprID"], "and", G.node[n2]["full_pprID"], file=sys.stderr)
                    print(fwd1rev2_can_heteroprime, fwd1rev2_frac_heteroduplex, file=sys.stderr)
                    print(edge_data["fwd1rev2_product_penalties"], file=sys.stderr)
                    print(fwd2rev1_can_heteroprime, fwd2rev1_frac_heteroduplex, file=sys.stderr)
                    print(edge_data["fwd2rev1_product_penalties"], file=sys.stderr)
                    
                G.add_edge(n1, n2, attr_dict=edge_data)
    return G
    


def formPool(G, max_clique_dyn_basename):
    primer_mishyb_ratio_cutoffs = [10, 10, 50, 50, 100, 100, 100]
    sum_product_penalties_cutoffs = [0.0001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.05]
    max_product_penalty_cutoffs = [1e-6, 1e-6, 1e-5, 1e-5, 1e-4, 1e-4, 1e-3]

    solutions = []
    
    for index in xrange(7):
        nodes_to_keep = []
        # Filter out nodes
        nodes_removed = []
        for n,D in G.nodes_iter(data=True):
            if (D["full_pprID"].startswith("UDTD")):
                nodes_to_keep.append(n)
                continue
                
            fwd_primer_seq, rev_primer_seq = n
            # For fwd_mishyb_ratio to be used, it needs to accumulate penalties only for where the corresponding rev primer has been kept.
            #fwd_mishyb_ratio = float(sum(map(itemgetter(1), D["fwd_primer_hyb_scores"])))/float(sum(map(itemgetter(1), D["fwd_primer_hyb_penalties"]))) 
            fwd_not_ok = D["fwd_has_hairpin"] or D["fwd_frac_homoduplex"]>0.01 or (D["fwd_can_self_prime"] and D["fwd_frac_homoduplex"]>0.001)

            try:
                rev_mishyb_ratio = float(sum(map(itemgetter(1), D["rev_primer_hyb_penalties"])))/float(sum(map(itemgetter(1), D["rev_primer_hyb_scores"])))
            except ZeroDivisionError as zde:
                if (D["full_pprID"].startswith("UDTD") or D["full_pprID"].startswith("ERCC") or D["full_pprID"].endswith("CTRL")):
                    rev_mishyb_ratio = float(sum(map(itemgetter(1), D["rev_primer_hyb_penalties"])))
                    #print >> sys.stderr, D["locus_group"][0], " ",
                else:
                    pdb.set_trace()
            rev_not_ok = D["rev_has_hairpin"] or D["rev_frac_homoduplex"]>0.01 or (D["rev_can_self_prime"] and D["rev_frac_homoduplex"]>0.001) or \
                         rev_mishyb_ratio > primer_mishyb_ratio_cutoffs[index]

            fwdrev_not_ok = D["fwdrev_frac_heteroduplex"]>0.01 or (D["fwdrev_can_heteroprime"] and D["fwdrev_frac_heteroduplex"]>0.001)

            if (D["full_pprID"].startswith("ERCC-00083") or D["full_pprID"].startswith("RBM17")):
                nodes_to_keep.append(n)
            elif (not(any([fwd_not_ok, rev_not_ok, fwdrev_not_ok]))):
                nodes_to_keep.append(n)
            else:
                nodes_removed.append(n)
                
        subG = G.subgraph(nodes_to_keep)

        # Filter out edges
        edges_to_remove = []
        for n1, n2, D in subG.edges_iter(data=True):
            fwd1rev2_not_ok = D["fwd1rev2_frac_heteroduplex"]>0.01 or (D["fwd1rev2_can_heteroprime"] and D["fwd1rev2_frac_heteroduplex"]>0.001)
            fwd2rev1_not_ok = D["fwd2rev1_frac_heteroduplex"]>0.01 or (D["fwd2rev1_can_heteroprime"] and D["fwd2rev1_frac_heteroduplex"]>0.001)
            product_penalties = map(lambda x: x[2]*x[3], D["fwd1rev2_product_penalties"]) + map(lambda x: x[2]*x[3], D["fwd2rev1_product_penalties"])
            exceeded_sum_product_penalties = sum(product_penalties) > sum_product_penalties_cutoffs[index]
            exceeded_max_product_penalties = max(product_penalties) > max_product_penalty_cutoffs[index] if (len(product_penalties)>0) else False

            if (set([subG.node[n1]["full_pprID"],subG.node[n2]["full_pprID"]]) == set(["chr7_OS10_1_F16R26", "chr16_OS2_1_F26R66"])):
                pdb.set_trace()
            if (fwd1rev2_not_ok or fwd2rev1_not_ok or exceeded_sum_product_penalties or exceeded_max_product_penalties):
                if (subG.node[n1]["full_pprID"] == "RBM17_CTRL" or subG.node[n2]["full_pprID"] == "RBM17_CTRL"):
                    #print >> sys.stderr, subG.node[n1]["full_pprID"], "\t", subG.node[n2]["full_pprID"]
                    #pdb.set_trace()
                    pass # *** TEMPORARY ***
                elif (subG.node[n1]["full_pprID"] == "ERCC-00083" or subG.node[n2]["full_pprID"] == "ERCC-00083"):
                    #print >> sys.stderr, subG.node[n1]["full_pprID"], "\t", subG.node[n2]["full_pprID"]
                    pass # *** TEMPORARY ***

                else:
                    edges_to_remove.append( (n1,n2) )

        subG.remove_edges_from(edges_to_remove)

        #  Attempt to find clique
        print("Graph %d: Removed %d nodes (explicitly) and %d edges" % (index, len(nodes_removed), len(edges_to_remove)), file=sys.stderr)
        print("Graph %d has %d nodes and %d edges" % (index, subG.number_of_nodes(), subG.size()), file=sys.stderr)
        output_edge_file = "%s/%s-%d.edges" % (temp_dir, max_clique_dyn_basename, index)
        op_edges = open(output_edge_file, 'w')
        for n1, n2 in subG.edges_iter():
            op_edges.write("e %d %d\n" % (subG.node[n1]["index"], subG.node[n2]["index"]))
        op_edges.close()

        try:
            mcqd_output = check_output(["/usr/local/src/mcqd/mcqd_custom", output_edge_file])
        except CalledProcessError as cpe:
            pdb.set_trace()

        if (mcqd_output.startswith("Maximum clique: ")):
            clq = map(int, mcqd_output.split("Maximum clique:")[1].split())
            print("INFO: Identified clique of size", len(clq), file=sys.stderr)
            if (len(clq) > 90):
                solutions.append( (index, subG, clq) )
        else:
            print("Unexpected output from mcqd_custom:\n ", mcqd_output, file=sys.stderr)
            pdb.set_trace()

        os.remove(output_edge_file)

    return solutions


def processSolution(G, solutions, output_place_picking):
    for G_index, soln_G, soln_clique in solutions:
        excluded_nodes = []
        soln_nodes = []

        for n in G.nodes():
            if (G.node[n]["index"] in soln_clique):
                soln_nodes.append(n)
            else:
                excluded_nodes.append(n)

        with open("%s-%d.txt" % (output_place_picking, G_index), 'w') as op:
            #op.write("Solution primer pairs\n")
            for n in soln_nodes:
                op.write("+ %s\n" % G.node[n]["full_pprID"])

            #op.write("Excluded primer pairs\n")
            for n in excluded_nodes:
                op.write("- %s\n" % G.node[n]["full_pprID"])


if (__name__ == "__main__"):
    params_descriptor, database_path, scratch_dir_root, transcriptome_fasta, unprocessed_transcriptome_fasta, \
        fwd_primers_1_tsv, rev_primers_1_tsv, fwd_primers_2_tsv, rev_primers_2_tsv, control_primers_fasta, \
        adapter_primers_fasta, ribosomal_loci, output_place_picking, max_clique_dyn_basename = sys.argv[1:]
    
    designParams.setParameters(params_descriptor)

    data_for_graph = False
    build_graph = True
    pick_primers = False

    temp_dir = tempfile.mkdtemp(prefix="screenForUnintendedProducts_",dir=scratch_dir_root)

    if (data_for_graph):
        high_stringency_loci = readLociFile(ribosomal_loci)

        fwd_rev_primers_1 = readAndPairPrimers(fwd_primers_1_tsv, rev_primers_1_tsv, "pool1")
        fwd_rev_primers_2 = readAndPairPrimers(fwd_primers_2_tsv, rev_primers_2_tsv, "pool2")
        fwd_rev_primers = fwd_rev_primers_1
        fwd_rev_primers.update(fwd_rev_primers_2)

        ctrl_primers = readPrimersFasta(control_primers_fasta)
        fwd_rev_primers.update(ctrl_primers)
        adapter_primers = readPrimersFasta(adapter_primers_fasta)
        fwd_rev_primers.update(adapter_primers)
        print("INFO: read", len(fwd_rev_primers), "primer pairs", file=sys.stderr)

        fwd_primers = map(lambda t: ((t[0], "fwd"), t[1][0]), fwd_rev_primers.items())
        rev_primers = map(lambda t: ((t[0], "rev"), t[1][1]), fwd_rev_primers.items())

        print("INFO: evaluating RT primers all-vs-all...", file=sys.stderr, end=" ")
        rev_primer_self_interaction, rev_primer_pair_interaction = compareRTPrimersPairwise(rev_primers, rev_primers)
        print("done", file=sys.stderr)

        print("INFO: evaluating DNAPol primers all-vs-all...", file=sys.stderr, end=" ")
        fwd_primer_self_interaction, fwd_primer_pair_interaction = compareDNAPolPrimersPairwise(fwd_primers, fwd_primers)
        print("done", file=sys.stderr)

        print("INFO: evaluating DNAPol primers vs RT primers...", file=sys.stderr, end=" ")
        fwdrev_primer_self_interaction, fwdrev_primer_pair_interaction = compareDNAPolPrimersPairwise(fwd_primers, rev_primers)
        print("done", file=sys.stderr)

        mfe_options = MFEOptions(0,1100) 
        transcriptomic_MFE = MFE(transcriptome_fasta, mfe_options)
        expected_primed_isoforms, expected_primer_products = getPrimersKnownTargets(fwd_rev_primers, database_path, transcriptomic_MFE)

        cPickle.dump((fwd_rev_primers, expected_primed_isoforms, expected_primer_products), open("fwd_rev_primers.pkl", "wb"))
        cPickle.dump((fwd_primer_self_interaction, fwd_primer_pair_interaction, rev_primer_self_interaction, rev_primer_pair_interaction,
                      fwdrev_primer_self_interaction, fwdrev_primer_pair_interaction), open("primer_interactions.pkl", "wb"))

        print("INFO: screening transcriptome for unintended products", file=sys.stderr)
        primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties = screenForUnintendedProducts(fwd_rev_primers, transcriptome_fasta,
                                                                                                                 expected_primed_isoforms, expected_primer_products,
                                                                                                                 temp_dir, high_stringency_loci)

        cPickle.dump((primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties), open("primer_and_product_penalties.pkl", 'wb'))

    if (build_graph):
        (fwd_rev_primers, expected_primed_isoforms, expected_primer_products) = cPickle.load(open("fwd_rev_primers.pkl", "rb"))
        (fwd_primer_self_interaction, fwd_primer_pair_interaction, rev_primer_self_interaction, rev_primer_pair_interaction,
         fwdrev_primer_self_interaction, fwdrev_primer_pair_interaction) = cPickle.load(open("primer_interactions.pkl", "rb"))
        primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties = cPickle.load(open("primer_and_product_penalties.pkl", 'rb'))

        # TEMP
        primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties = screenForUnintendedProducts(fwd_rev_primers, transcriptome_fasta,
                                                                                                                 expected_primed_isoforms, expected_primer_products,
                                                                                                                 temp_dir, set())
        # Needed? :fwdrev_primer_self_interaction, 
        G = buildPrimerPairsGraph(fwd_rev_primers, fwd_primer_self_interaction, fwd_primer_pair_interaction, rev_primer_self_interaction, rev_primer_pair_interaction,
                                  fwdrev_primer_pair_interaction, primer_hyb_scores, primer_hyb_penalties, product_scores, product_penalties)

        cPickle.dump(G, open("G_w_edge_data.pkl", "wb"))

    if (pick_primers):
        G = cPickle.load(open("G_w_edge_data.pkl", "rb"))
        solutions = formPool(G, max_clique_dyn_basename, temp_dir)
        processSolution(G, solutions, output_place_picking)

    os.remove(temp_dir)

    sys.exit(0)
