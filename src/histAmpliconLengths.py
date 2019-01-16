#!/usr/bin/env python

import cPickle
from collections import Counter, defaultdict
import sqlite3
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import pdb

def groupPrimersByExperiment(array_primers_fasta):
    ip = open(array_primers_fasta, 'r')

    fwd_primers = defaultdict(set)
    rev_primers = defaultdict(set)

    for line in ip:
        line = line.strip()
        if (line[0] == '>'):
            partition, direction, ppr_ID = line[1:].split("_")
            ppr_IDs = ppr_ID.split(';')
        elif (not ppr_ID.startswith('ERCC')):
            primer_seq = line[0:-18].split("NNNNNNNN")[1]
            for ppr_ID in ppr_IDs:
                key = (primer_seq, ppr_ID)
                if (direction == "Fwd"):
                    if (fwd_primers.has_key(key)):
                        print >> sys.stderr, "Fwd duplicate of ", key
                    fwd_primers[key].add(partition)
                else:
                    if (rev_primers.has_key(key)):
                        print >> sys.stderr, "Rev duplicate of ", key
                    rev_primers[key].add(partition)

    ip.close()

    return (fwd_primers, rev_primers)


def partitionPrimers(database_path, fwd_primers, rev_primers):
    partitioned_primers = defaultdict(list)

    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()
    cursor.execute("SELECT chrom,OS_ID,ppr_ID,fwd_primer_seq,rev_primer_seq,fwd_genomic_positions,rev_genomic_positions FROM primers")
    for chrom,OS_ID,ppr_ID,fwd_primer_seq,rev_primer_seq,fwd_genomic_positions,rev_genomic_positions in cursor.fetchall():
        fwd_key = (fwd_primer_seq, ppr_ID)
        rev_key = (rev_primer_seq, ppr_ID)
        if (fwd_primers.has_key(fwd_key)):
            try:
                assert (rev_primers.has_key(rev_key))
            except AssertionError:
                pdb.set_trace()
            fwd_partitions = fwd_primers[fwd_key]
            rev_partitions = rev_primers[rev_key]
            common_part = fwd_partitions.intersection(rev_partitions)
            assert (len(common_part) > 0)
            if (len(common_part) == 1): # Skip the > 1 case
                partition = list(common_part)[0]
                all_positions = map(int, fwd_genomic_positions.split()) + map(int, rev_genomic_positions.split())
                min_pos = min(all_positions)
                max_pos = max(all_positions)
                partitioned_primers[(chrom, OS_ID)].append( (ppr_ID, min_pos, max_pos, partition) )
        else:
            try:
                assert (not rev_primers.has_key(rev_key))
            except AssertionError:
                print >> sys.stderr, "Unmatched Rev : ", rev_key

    conn.close()

    return partitioned_primers


def getAmpliconLengthsPerPartition(partitioned_primers, database_path):
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    amplicon_lengths = defaultdict(list)
    
    counter = 0
    for (chrom, OS_ID) in partitioned_primers.keys():
        counter += 1
        if (counter%10 == 0):
            print >> sys.stderr, "\r%d/%d" % (counter, len(partitioned_primers.keys()))

        cursor.execute("SELECT predesign_pkl_file FROM predesigns WHERE chrom=? AND OS_ID=?", (chrom, OS_ID))
        predesign_pkl_file = cursor.fetchone()[0]
        olap_set = cPickle.load(open(predesign_pkl_file, 'rb'))
        olap_set.setIsoformsCoordCorrespondData()

        for ppr_ID, min_pos, max_pos, partition in partitioned_primers[(chrom, OS_ID)]:
            for isoform in olap_set.isoforms:
                if (isoform.strand == '+'):
                    amplicon_len = isoform.getSequenceLength(min_pos, max_pos)
                else:
                    amplicon_len = isoform.getSequenceLength(max_pos, min_pos)

                if (amplicon_len != None):
                    amplicon_lengths[partition].append(amplicon_len)

        olap_set.unsetIsoformsCoordCorrespondData()

    conn.close()

    return amplicon_lengths


def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'


def writeHistograms(amplicon_lengths, hist_prefix):
    formatter = FuncFormatter(to_percent)
    for partition, lengths in amplicon_lengths.items():
        lengths = filter(lambda x:x<1000, lengths)
        fig_name = "%s_%s.png" % (hist_prefix, partition)
        plt.hist(lengths, bins=101, normed=True)
        plt.gca().yaxis.set_major_formatter(formatter)    
        plt.title("Subpool %s" % partition)
        plt.xlabel("Amplicon Length (nt)")
        plt.ylabel("Percent of All Amplicons")
        plt.savefig(fig_name)
        plt.gcf().clear()


if (__name__ == "__main__"):
    database_path = "/raid1/projects/CoveringPrimersSets/PM/CGDBv1.0/v1.0/first-orders-main.db"
    predesigns_pkl_dir = "/raid1/projects/CoveringPrimersSets/PM/CGDBv1.0/v1.0/first-orders-pkl_predesign"
    array_primers_fasta = "/raid1/projects/CoveringPrimersSets/PM/CGDBv1.0/v1.0/Isommune_array_design_1.fa"
    hist_prefix = "/raid1/projects/CoveringPrimersSets/PM/CGDBv1.0/v1.0/Isommune_array_design_hist"
    
    #fwd_primers, rev_primers = groupPrimersByExperiment(array_primers_fasta)
    #partitioned_primers = partitionPrimers(database_path, fwd_primers, rev_primers)
    #amplicon_lengths = getAmpliconLengthsPerPartition(partitioned_primers, database_path)
    #cPickle.dump(amplicon_lengths, open("amplicon_lengths.pkl", "wb"))

    amplicon_lengths = cPickle.load(open("amplicon_lengths.pkl", "rb"))

    print >> sys.stderr, "INFO: writing histograms"
    writeHistograms(amplicon_lengths, hist_prefix)
    
    sys.exit(0)
