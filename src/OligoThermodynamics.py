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


class OligoThermodynamics(object):
    def __init__(self):
        rxn_temp_C = designParams.DNAPOL_rxn_temp
        rxn_temp_K = 273.15 + rxn_temp_C
        self.thermo_analysis = primer3.thermoanalysis.ThermoAnalysis(mv_conc=designParams.monovalent_salt_conc, dv_conc=designParams.divalent_salt_conc,
                                                                     dntp_conc=designParams.dntp_conc, temp_c=rxn_temp_C, dna_conc=designParams.input_primer_conc_nM)

        RT = 1.98717 * rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)
        # dG units are expected to be cal/mol and NOT kcal/mol
        self.calc_frac_duplexed = lambda dG: (designParams.input_primer_conc_M * exp(-dG/RT))/(1 + designParams.input_primer_conc_M * exp(-dG/RT))


    def calcHeterodimer(self, seq1, seq2):
        return self.thermo_analysis.calcHeterodimer(seq1.upper(), seq2.upper())


    def calcEndStability(self, seq1, seq2):
        return self.thermo_analysis.calcEndStability(seq1.upper(), seq2.upper())


    def calcFracDuplexed(self, dG):
        return self.calc_frac_duplexed(dG)


    def calcPrimerThermoDetails(self, primer_type, primer_seq, tagged=True):
        primer_seq = primer_seq.upper()
        
        primer_Tm = self.thermo_analysis.calcTm(primer_seq)
        primer_exact_hyb = self.thermo_analysis.calcHeterodimer(primer_seq, primer_seq[::-1].translate(DNA_complement_table))
        primer_frac_duplexed = self.calc_frac_duplexed(primer_exact_hyb.dg)

        if (primer_type == "fwd"):
            fwd_self_hairpin = self.thermo_analysis.calcHairpin(primer_seq)
            fwd_any_hyb = self.thermo_analysis.calcHomodimer(primer_seq)
            fwd_end_hyb = self.thermo_analysis.calcEndStability(primer_seq, primer_seq)
            fwd_homo_dG = sum(list(filter(lambda x: x<0, [fwd_any_hyb.dg, fwd_end_hyb.dg])))
            frac_fwd_duplexed = self.calc_frac_duplexed(fwd_homo_dG)

            if (tagged):
                fwd_tag_any_hyb = self.thermo_analysis.calcHeterodimer(primer_seq, designParams.fwd_primer_adapter_tag)
                fwd_tag_end_hyb = self.thermo_analysis.calcEndStability(primer_seq, designParams.fwd_primer_adapter_tag)
                fwd_tag_dG = sum(list(filter(lambda x: x<0, [fwd_tag_any_hyb.dg, fwd_tag_end_hyb.dg])))
                frac_fwd_tag_duplexed = self.calc_frac_duplexed(fwd_tag_dG)

                primer_has_hairpin = (fwd_self_hairpin.structure_found or frac_fwd_tag_duplexed > 0.01)
                primer_penalty = frac_fwd_duplexed + frac_fwd_tag_duplexed
            else:
                primer_has_hairpin = fwd_self_hairpin.structure_found
                primer_penalty = frac_fwd_duplexed

        elif (primer_type == "rev"):
            rev_self_hairpin = self.thermo_analysis.calcHairpin(primer_seq)
            rev_any_hyb = self.thermo_analysis.calcHomodimer(primer_seq)
            rev_end_hyb = self.thermo_analysis.calcEndStability(primer_seq, primer_seq)
            rev_homo_dG = sum(list(filter(lambda x: x<0, [rev_any_hyb.dg, rev_end_hyb.dg])))
            frac_rev_duplexed = self.calc_frac_duplexed(rev_homo_dG)

            if (tagged):
                rev_tag_any_hyb = self.thermo_analysis.calcHeterodimer(primer_seq, designParams.rev_primer_adapter_tag)
                rev_tag_end_hyb = self.thermo_analysis.calcEndStability(primer_seq, designParams.rev_primer_adapter_tag)
                rev_tag_dG = sum(list(filter(lambda x: x<0, [rev_tag_any_hyb.dg, rev_tag_end_hyb.dg])))
                frac_rev_tag_duplexed = self.calc_frac_duplexed(rev_tag_dG)

                primer_has_hairpin = (rev_self_hairpin.structure_found or frac_rev_tag_duplexed > 0.01)
                primer_penalty = frac_rev_duplexed + frac_rev_tag_duplexed
            else:
                primer_has_hairpin = rev_self_hairpin.structure_found
                primer_penalty = frac_rev_duplexed

        else:
            print("ERROR: unknown primer type -> %s" % primer_type, file=sys.stderr, flush=True)
            sys.exit(1)

        return (primer_Tm, primer_frac_duplexed, primer_has_hairpin, primer_penalty)


    def calcPrimerPairThermoDetails(self, fwd_seq, rev_seq, include_tag=True):
        fwd_seq = fwd_seq.upper()
        rev_seq = rev_seq.upper()
        
        fwdrev_any_hyb = self.thermo_analysis.calcHeterodimer(fwd_seq, rev_seq)
        fwdrev_end_hyb = self.thermo_analysis.calcEndStability(fwd_seq, rev_seq)
        revfwd_end_hyb = self.thermo_analysis.calcEndStability(rev_seq, fwd_seq)

        if (include_tag):
            # Check fwd_seq against rev adapter and rev_seq against fwd adapter
            rev_fwdtag_any_hyb = self.thermo_analysis.calcHeterodimer(rev_seq, designParams.fwd_primer_adapter_tag)
            rev_fwdtag_end_hyb = self.thermo_analysis.calcEndStability(rev_seq, designParams.fwd_primer_adapter_tag)
            min_rev_fwdtag_dG = min(rev_fwdtag_any_hyb.dg, rev_fwdtag_end_hyb.dg)

            fwd_revtag_any_hyb = self.thermo_analysis.calcHeterodimer(fwd_seq, designParams.rev_primer_adapter_tag)
            fwd_revtag_end_hyb = self.thermo_analysis.calcEndStability(fwd_seq, designParams.rev_primer_adapter_tag)
            min_fwd_revtag_dG = min(fwd_revtag_any_hyb.dg, fwd_revtag_end_hyb.dg)

            # Use sum of all that are unique and negative
            hetero_dG = sum(list(filter(lambda x: x<0, [fwdrev_any_hyb.dg, fwdrev_end_hyb.dg, revfwd_end_hyb.dg, min_rev_fwdtag_dG, min_fwd_revtag_dG])))
        else:
            hetero_dG = sum(list(filter(lambda x: x<0, [fwdrev_any_hyb.dg, fwdrev_end_hyb.dg, revfwd_end_hyb.dg])))

        frac_hetero_duplexed = self.calc_frac_duplexed(hetero_dG)

        return frac_hetero_duplexed

