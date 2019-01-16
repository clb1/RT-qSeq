# Module for setting global variables for use in designSets.py and designSetsClasses.py

import sys

def setParameters(descriptor):
    global params_descriptor

    global is_directional

    global RT_rxn_temp
    global DNAPOL_rxn_temp

    global monovalent_salt_conc
    global divalent_salt_conc
    global dntp_conc
    global input_primer_conc_nM # nanomolar
    global input_primer_conc_M  # molar

    global min_primer_len
    global max_primer_len
    global opt_primer_len

    global min_primer_Tm
    global max_primer_Tm
    global opt_primer_Tm

    global read_len

    global amplicon_min_len
    global amplicon_max_len
    global amplicon_opt_len

    # PRIMER_MAX_SELF_ANY_TH
    # PRIMER_MAX_SELF_END_TH
    # PRIMER_MAX_HAIRPIN_TH
    global RT_primer_max_self_any_th
    global RT_primer_max_self_end_th
    global RT_primer_max_hairpin_th

    global DNAPol_primer_max_self_any_th
    global DNAPol_primer_max_self_end_th
    global DNAPol_primer_max_hairpin_th


    # PRIMER_PAIR_MAX_COMPL_ANY_TH
    # PRIMER_PAIR_MAX_COMPL_END_TH
    global primer_pair_max_compl_any_th
    global primer_pair_max_compl_end_th

    global max_primer_promiscuity

    global confusion_amplicon_max_len

    global max_RT_misprime_temp
    global max_DNAPol_misprime_temp

    global disallowed_primer_subseqs
    global disallowed_primer_3p_subseqs

    global fwd_primer_adapter_tag
    global rev_primer_adapter_tag

    global min_primer_adapter_hairpin_DeltaG
    global max_fwd_adapter_hairpin_Tm
    global max_rev_adapter_hairpin_Tm


    if (descriptor == "qPCR"):
        params_descriptor = descriptor

        is_directional = False

        RT_rxn_temp = 57
        DNAPOL_rxn_temp = 57

        monovalent_salt_conc = 50
        divalent_salt_conc = 1.5
        dntp_conc = 0.25
        input_primer_conc_nM = 300
        input_primer_conc_M = 3*10**-7

        min_primer_len = 18
        max_primer_len = 23
        opt_primer_len = 20

        RT_primer_max_self_any_th = RT_rxn_temp - 10
        RT_primer_max_self_end_th = RT_rxn_temp - 10
        RT_primer_max_hairpin_th = RT_rxn_temp - 10

        DNAPol_primer_max_self_any_th = DNAPOL_rxn_temp - 10
        DNAPol_primer_max_self_end_th = DNAPOL_rxn_temp - 10
        DNAPol_primer_max_hairpin_th = DNAPOL_rxn_temp - 10

        primer_pair_max_compl_any_th = DNAPOL_rxn_temp - 10
        primer_pair_max_compl_end_th = DNAPOL_rxn_temp - 10

        min_primer_Tm = 58
        max_primer_Tm = 63
        opt_primer_Tm = 60

    elif (descriptor == "dir_default" or descriptor == "default"):
        params_descriptor = descriptor

        is_directional = True

        RT_rxn_temp = 50
        DNAPOL_rxn_temp = 57

        monovalent_salt_conc = 50
        divalent_salt_conc = 1.5
        dntp_conc = 0.25
        input_primer_conc_nM = 20
        input_primer_conc_M = 2*10**-8

        min_primer_len = 18
        max_primer_len = 23
        opt_primer_len = 20

        RT_primer_max_self_any_th = RT_rxn_temp - 10
        RT_primer_max_self_end_th = RT_rxn_temp - 10
        RT_primer_max_hairpin_th = RT_rxn_temp - 10

        DNAPol_primer_max_self_any_th = DNAPOL_rxn_temp - 10
        DNAPol_primer_max_self_end_th = DNAPOL_rxn_temp - 10
        DNAPol_primer_max_hairpin_th = DNAPOL_rxn_temp - 10

        primer_pair_max_compl_any_th = DNAPOL_rxn_temp - 10
        primer_pair_max_compl_end_th = DNAPOL_rxn_temp - 10

        min_primer_Tm = 58
        max_primer_Tm = 63
        opt_primer_Tm = 60


    elif (descriptor == "dir_pe250_700" or descriptor == "PE250_700"):
        params_descriptor = descriptor

        is_directional = True

        RT_rxn_temp = 50
        DNAPOL_rxn_temp = 57

        monovalent_salt_conc = 50
        divalent_salt_conc = 1.5
        dntp_conc = 0.25
        input_primer_conc_nM = 20
        input_primer_conc_M = 2*10**-8

        min_primer_len = 18
        max_primer_len = 23
        opt_primer_len = 20

        RT_primer_max_self_any_th = RT_rxn_temp - 10
        RT_primer_max_self_end_th = RT_rxn_temp - 10
        RT_primer_max_hairpin_th = RT_rxn_temp - 10

        DNAPol_primer_max_self_any_th = DNAPOL_rxn_temp - 10
        DNAPol_primer_max_self_end_th = DNAPOL_rxn_temp - 10
        DNAPol_primer_max_hairpin_th = DNAPOL_rxn_temp - 10

        primer_pair_max_compl_any_th = DNAPOL_rxn_temp - 10
        primer_pair_max_compl_end_th = DNAPOL_rxn_temp - 10

        min_primer_Tm = 57
        max_primer_Tm = 63
        opt_primer_Tm = 60

        read_len = 250

        amplicon_min_len = 150
        amplicon_max_len = 700
        amplicon_opt_len = 350

        max_primer_promiscuity = 200
                
        confusion_amplicon_max_len = 900

        max_RT_misprime_temp = 47
        max_DNAPol_misprime_temp = 47

        # Fwd and reverse complement BspQI restriction motifs that would cause
        # all primer oligos to be degraded during OligoPool procedure.
        disallowed_primer_subseqs = [] # "GCTCTTC", "CGAGAAG"

        # Sequences that are not allowed at the 5' end of the primer sequences,
        # otherwise they would too often form a BspQI motif (or its complement) with the UMIs.
        # Based on the assumption that the UMI structure is:
        # P5: NNNNNNAT
        # P7: NNNNNCNNA
        # no disallowed 5' subseqs
        #disallowed_primer_5p_subseqs = []

        # Sequences that are not allowed at the 3' end of the primer sequences,
        # otherwise they could form a BspQI motif with an OligoPool subpool reverse primer region
        disallowed_primer_3p_subseqs = [] # "CGAGA", "GCTCT"

        # Normal Illumina adapter tags have a BspQI restriction motif:
        #
        #                      BspQI
        #                     1234567
        # P5_ADAPTER_TAG = CGACGCTCTTCCGATCT
        # P7_ADAPTER_TAG = GTGTGCTCTTCCGATC
        #
        # Since BspQI is used in the OligoPool amplification, a modified tag
        # sequence must be used in the primer design.
        #
        # The Illumina adapter tags with a G->C change at BspQI position 1
        #fwd_primer_adapter_tag = "CGACCCTCTTCCGATCT" # P5_ADAPTER_TAG = 
        #rev_primer_adapter_tag = "GTGTCCTCTTCCGATCT"  # P7_ADAPTER_TAG = 
                                     
        # These *do not* have the G->C change for BspQI.
        fwd_primer_adapter_tag = "CTACACGACGCTCTTCCGATCT"
        rev_primer_adapter_tag = "CAGACGTGTGCTCTTCCGATCT"

        min_primer_adapter_hairpin_DeltaG = -4
        max_fwd_adapter_hairpin_Tm = 47
        max_rev_adapter_hairpin_Tm = 47

    else:
        print("ERROR: unrecognized descriptor for global variable settings -> %s" % descriptor, file=sys.stderr)
        sys.exit(0)
