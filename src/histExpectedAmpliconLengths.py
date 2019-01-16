#!/usr/bin/env python

from collections import defaultdict, Counter
from itertools import chain
from pyfaidx import Fasta
import cPickle
import sqlite3
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import designParams

import pdb

def readPrimers(primers_tsv):
    primers = defaultdict(dict)

    ip = open(primers_tsv, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        primer_seq = fields[2].split('N')[-1]
        ppr_specs = fields[3].split(';')
        for ppr_spec in ppr_specs:
            if (ppr_spec.startswith('ERCC')):
                primers[(ppr_spec,ppr_spec)][ppr_spec] = primer_seq
            else:
                elems = ppr_spec.split('_')
                chrom_OSID = (elems[0], elems[1])
                experiment_num_and_PPR_ID = (int(elems[2]), elems[3])
                assert (experiment_num_and_PPR_ID not in primers[chrom_OSID])
                primers[chrom_OSID][experiment_num_and_PPR_ID] = primer_seq
    ip.close()

    return primers


def getAmpliconLengths(params_descriptor, genome_ref, fwd_primers, rev_primers):
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    assert(set(fwd_primers.keys()) == set(rev_primers.keys())), "Mismatch between (chrom, OS_ID) sets for fwd & rev primers"

    amplicon_lens = []

    counter = 0
    for (chrom, OS_ID) in fwd_primers.keys():
        print >> sys.stderr, "%s %s" % (chrom, OS_ID)
        if (chrom.startswith('ERCC')):
            continue # For now
    counter += 1
    if (counter%10 == 0):
        print >> sys.stderr, "\r%d/%d" % (counter, len(fwd_primers.keys()))

    cursor.execute("SELECT predesign_pkl_file FROM Predesigns WHERE chrom=? AND OS_ID=? AND params_descriptor=?", (chrom, OS_ID, params_descriptor))
    predesign_pkl_file = cursor.fetchone()[0]
    olap_set = cPickle.load(open(predesign_pkl_file, 'rb'))
    olap_set.setIsoformsCoordCorrespondData()

        assert (set(fwd_primers[(chrom,OS_ID)].keys()) == set(rev_primers[(chrom,OS_ID)].keys())), \
            "Mismatch between (experiment num, PPR_ID) sets for %s, %s" % (chrom, OS_ID)

        for experiment_num, PPR_ID in fwd_primers[(chrom,OS_ID)]:
            fwd_primer_seq = fwd_primers[(chrom,OS_ID)][(experiment_num, PPR_ID)]
            rev_primer_seq = rev_primers[(chrom,OS_ID)][(experiment_num, PPR_ID)]
            select_tup = (chrom, OS_ID, params_descriptor, experiment_num, PPR_ID, fwd_primer_seq, rev_primer_seq)
            cursor.execute("SELECT local_target_isoform_group from PrimerPair WHERE chrom=? AND OS_ID=? AND params_descriptor=? AND part_num=? AND PPR_ID=? and fwd_primer_seq=? AND rev_primer_seq=?", select_tup)

            results = cursor.fetchall()
            assert (len(results) == 1)

            tigs = results[0][0].split(';')
            expected_isoform_IDs = set(chain.from_iterable(map(lambda tig: tig.split(','), tigs)))
            amplicon_lens_for_primers = set()
        for isoform in olap_set.isoforms:
                isoform_ID = isoform.getCGDBName()
                isoform.setSequence(genome_ref)
                amplicon = isoform.getPCRProduct(fwd_primer_seq, rev_primer_seq)
                if (amplicon == None):
                    assert (isoform_ID not in expected_isoform_IDs), "Primers failed to give amplicon for isoform %s" % isoform_ID
                else:
                    try:
                        assert (isoform_ID in expected_isoform_IDs or len(amplicon) > designParams.confusion_amplicon_max_len), \
                            "Primers gave unexpected amplicon for isoform %s" % isoform_ID
                    except AssertionError, ae:
                        pdb.set_trace()
                    # WARNING: assumption that different isoform/isoform group primer products are all going to have different lengths.
                    # This obscures cases where that doesn't hold. Very minor impact on histogram, if any, but not absolutely ideal.
                    # But doesn't overcount for same amplicons generated for amplicons in same isoform group.
                    amplicon_lens_for_primers.add( len(amplicon) )

            amplicon_lens.extend( list(amplicon_lens_for_primers) )
            
    olap_set.unsetIsoformsCoordCorrespondData()

    conn.close()

    return amplicon_lens

    
def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'


def writeHistograms(amplicon_lens, params_descriptor, pool_label, fig_name):
    formatter = FuncFormatter(to_percent)

    plt.hist(amplicon_lens, bins=101, normed=True)
    plt.gca().yaxis.set_major_formatter(formatter)    
    plt.title("%s, %s" % (params_descriptor, pool_label))
    plt.xlabel("Amplicon Length (nt)")
    plt.ylabel("Percent of All Amplicons")
    plt.savefig(fig_name)
    plt.gcf().clear()


if (__name__ == "__main__"):
    params_descriptor, pool_label, database_path, genome_fasta, fwd_primers_tsv, rev_primers_tsv, output_png = sys.argv[1:]

    designParams.setParameters(params_descriptor)
    genome_ref = Fasta(genome_fasta, as_raw=True, sequence_always_upper=True)

    fwd_primers = readPrimers(fwd_primers_tsv)
    rev_primers = readPrimers(rev_primers_tsv)
    
    amplicon_lens = getAmpliconLengths(params_descriptor, genome_ref, fwd_primers, rev_primers)
    
    writeHistograms(amplicon_lens, params_descriptor, pool_label, output_png)

    sys.exit(0)
