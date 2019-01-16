#!/usr/bin/env python3

#import cPickle
from itertools import product
from operator import itemgetter
#import matplotlib.pyplot as plt
import sys
#from subprocess import Popen, PIPE
import designParams
from designSetsClasses import OligoThermodynamics

import pdb

def readPrimerSeqs(primers_file):
    primers = []
    ip = open(primers_file, 'r')
    for line in ip:
        line = line.strip()
        if (line != ''):
            primers.append(line)
    ip.close()
    return primers


def readPrimerPairs(primers_file):
    fwd_primers = []
    rev_primers = []
    
    ip = open(primers_file, 'r')
    for line in ip:
        if (line.startswith("ISM0001_2")):
            fields = line.split("\t")
            seq = fields[3].split('N')[-1]
            if (fields[0].endswith("fwd")):
                fwd_primers.append( (fields[2], seq) )
            else:
                rev_primers.append( (fields[2], seq) )
    ip.close()
    
    return (fwd_primers, rev_primers)


def doAllvAll2(A_primers, B_primers, name):
    relevant_results = []
    pdb.set_trace()
    for (A_id, A_seq), (B_id, B_seq) in product(A_primers, B_primers):
        any_hyb_AB = oligo_thermo.calcHeterodimer(A_seq, B_seq)
        end_hyb_AB = oligo_thermo.calcEndStability(A_seq, B_seq)
        any_hyb_BA = oligo_thermo.calcHeterodimer(B_seq, A_seq)
        end_hyb_BA = oligo_thermo.calcEndStability(B_seq, A_seq)
        min_dG = min(any_hyb_AB.dg, end_hyb_AB.dg, any_hyb_BA.dg, end_hyb_BA.dg)
        frac_duplexed = oligo_thermo.calcFracDuplexed(min_dG)
        if (frac_duplexed > 0.01):
            relevant_results.append( (A_id, B_id, A_seq, B_seq, frac_duplexed) )

    return relevant_results


def writeRelevantResults(relevant_results, output_file):
    relevant_results.sort(key=itemgetter(4), reverse=True)
    op = open(output_file, 'w')
    for tup in relevant_results:
        op.write("%s\t%s\t%s\t%s\t%3.2f\t" % tup)
    op.close()
    

def doAllvAll(A_primers, B_primers, name):
    common_primer3_params = ["PRIMER_TASK=check_primers",
                             "PRIMER_PICK_LEFT_PRIMER=1",
                             "PRIMER_PICK_INTERNAL_OLIGO=0",
                             "PRIMER_PICK_RIGHT_PRIMER=1",
                             "PRIMER_EXPLAIN_FLAG=1",
                             "PRIMER_MIN_SIZE=%d" % designParams.min_primer_len,
                             "PRIMER_MAX_SIZE=%d" % designParams.max_primer_len,
                             "PRIMER_OPT_SIZE=%d" % designParams.opt_primer_len,
                             "PRIMER_NUM_RETURN=1",
                             "PRIMER_SALT_MONOVALENT=%f" % designParams.monovalent_salt_conc,
                             "PRIMER_SALT_DIVALENT=%f" % designParams.divalent_salt_conc,
                             "PRIMER_DNTP_CONC=%f" % designParams.dntp_conc,
                             "PRIMER_DNA_CONC=%d" % designParams.input_primer_conc,
                             "PRIMER_MIN_TM=%d" % designParams.min_primer_Tm,
                             "PRIMER_MAX_TM=%d" % designParams.max_primer_Tm,
                             "PRIMER_OPT_TM=%d" % designParams.opt_primer_Tm,
                             "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/local/src/primer3-2.3.5/src/primer3_config/",
                             "=\n"]

    disallowed = set([("GGATTCAGCACCTTGGCA", "TGCCAAGGTGCTGAATCC")])

    primer3_input_record_strings = []
    counter = 0
    for A_seq, B_seq in product(A_primers, B_primers):
        counter += 1
        if ((A_seq, B_seq) not in disallowed):
            primer3_params = ["SEQUENCE_ID=%d" % counter,
                              "SEQUENCE_PRIMER=%s" % A_seq,
                              "SEQUENCE_PRIMER_REVCOMP=%s" % B_seq]
            primer3_params += common_primer3_params
            primer3_input_record_strings.append( "\n".join(primer3_params) )

    print >> sys.stderr, "Calling Primer3 with %d records..." % len(primer3_input_record_strings),
    primer3 = Popen("/usr/local/src/primer3-2.3.5/src/primer3_core", shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    primer3_input_string = "".join(primer3_input_record_strings)
    stdout_results, stderr_results = primer3.communicate(input=primer3_input_string)
    print >> sys.stderr, "done"

    with open("%s_primer3_input.txt" % name, 'w') as op:
        op.write(primer3_input_string)

    try:
        assert (primer3.returncode == 0), "Primer3_core return code was %d." % primer3.returncode
    except AssertionError:
        pdb.set_trace()

    penalties = []
    for primer3_record in stdout_results.split("=\n")[0:-1]:
        primer_search_results = dict( map(lambda x: x.strip().split("="), primer3_record.split("\n")[0:-1]) )
        if (primer_search_results.has_key("PRIMER_ERROR")):
            print >> sys.stderr, "Primer3 Error: %s" % primer_search_results["PRIMER_ERROR"]
        elif (primer_search_results["PRIMER_PAIR_NUM_RETURNED"] == "1"):
            pair_penalty = float(primer_search_results["PRIMER_PAIR_0_PENALTY"])
            A_penalty = float(primer_search_results["PRIMER_LEFT_0_PENALTY"])
            B_penalty = float(primer_search_results["PRIMER_RIGHT_0_PENALTY"])

            penalties.append( (A_penalty, B_penalty, pair_penalty) )
        elif (primer_search_results["PRIMER_PAIR_NUM_RETURNED"] == "0"):
            print >> sys.stderr, "WARNING: failure for %s" % primer_search_results["SEQUENCE_ID"]
        else:
            print >> sys.stderr, "WARNING: unexpected result for %s" % primer_search_results["SEQUENCE_ID"]
            pdb.set_trace()

    return penalties


def plotPenalties(penalties, output_basename, title_descr):

    # Histogram the penalties
    n, bins, patches = plt.hist(penalties, 50, normed=1, facecolor='green', alpha=0.75)

    plt.xlabel('Pair penalty')
    plt.ylabel('Probability')
    plt.title("%s %s" % (output_basename,title_descr))
    plt.axis([0, 15, 0, 0.5])
    plt.grid(True)

    plt.savefig("%s_%s.png" % (output_basename,title_descr))
    plt.close()


if (__name__ == "__main__"):
    param_settings_descriptor, primers_file, output_file = sys.argv[1:] # label, fwd_primers_file, rev_primers_file, output_basename

    designParams.setParameters(param_settings_descriptor)
    oligo_thermo = OligoThermodynamics()

    fwd_primers, rev_primers = readPrimerPairs(primers_file)
    relevant_results = doAllvAll2(fwd_primers, rev_primers, oligo_thermo)
    writeRelevantResults(relevant_results, output_file)

    fwd_primers = readPrimerSeqs(fwd_primers_file)
    rev_primers = readPrimerSeqs(rev_primers_file)

    penalties = doAllvAll(fwd_primers, rev_primers, "%s_fwd-rev" % label)
    plotPenalties(map(itemgetter(2), penalties), output_basename, "Fwd-Rev Pair Penalties")

    penalties = doAllvAll(fwd_primers, fwd_primers, "%s_fwd-fwd" % label)
    plotPenalties(map(itemgetter(2), penalties), output_basename, "Fwd-Fwd Pair Penalties")

    penalties = doAllvAll(rev_primers, rev_primers, "%s_rev-rev" % label)
    plotPenalties(map(itemgetter(2), penalties), output_basename, "Rev-Rev Pair Penalties")

    #cPickle.dump(penalties, open("%s_penalties.pkl" % output_basename, 'wb'))
    #penalties = cPickle.load(open("%s_penalties.pkl" % output_basename, 'rb'))

    sys.exit(0)
