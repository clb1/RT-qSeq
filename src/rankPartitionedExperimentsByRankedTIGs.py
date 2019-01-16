#!/usr/bin/env python

from operator import itemgetter
from scipy.stats import rankdata
import sqlite3
import sys

import pdb


def readRankedTIGs(ranked_TIGs_file):
    ranked_TIGs = []
    ip = open(ranked_TIGs_file, 'r')
    for line in ip:
        ranked_TIGs.append( line.strip().split() )
    ip.close()
    return ranked_TIGs


def getRankedExperimentPartitions(database_path, params_descriptor, ranked_TIGs):
    ranked_partitions = []
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    num_found_partition = 0
    for TIG in ranked_TIGs:
        partition = set()
        for isoform in TIG:
            cursor.execute("SELECT chrom,OS_ID,part_num FROM ExperimentDesignPart WHERE params_descriptor = ? AND tigs LIKE ?", (params_descriptor,'%'+isoform+'%'))
            entries = cursor.fetchall()

            if (len(entries) == 1):
                partition.update(entries)

        if (len(partition) == 1):
            num_found_partition += 1
            TIG_str = ",".join(TIG)
            ranked_partitions.append( list(partition)[0]+(TIG_str,) )
    conn.close()

    print >> sys.stderr, "INFO: found partition for %d of %d ranked TIGs" % (num_found_partition, len(ranked_TIGs))
    return ranked_partitions


def writeRankedPartitions(ranked_partitions, output_tsv):
    op = open(output_tsv, 'w')
    op.write("#chromosome\tOS_ID\tpartition_number\tTIG\n")
    for tup in ranked_partitions:
        op.write("%s\t%s\t%d\t%s\n" % tup)
    op.close()


if (__name__ == "__main__"):
    params_descriptor, ranked_TIGs_file, database_path, output_tsv = sys.argv[1:]

    ranked_TIGs = readRankedTIGs(ranked_TIGs_file)

    ranked_partitions = getRankedExperimentPartitions(database_path, params_descriptor, ranked_TIGs)

    writeRankedPartitions(ranked_partitions, output_tsv)

    sys.exit(0)
