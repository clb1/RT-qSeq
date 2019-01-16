#!/usr/bin/env python

import pdb
import sys


def readSelectionLoci(loci_file):
    with open(loci_file, 'r') as ip:
        loci = set(ip.read().split("\n"))
    loci = map(lambda x: x.strip(), loci)
    loci = filter(lambda x: x != '', loci)
    return loci
        

def selectTIGs(selection_loci, tigs_file):
    selected_tigs = []
    utilized_loci = set()

    ip = open(tigs_file, 'r')
    for line in ip:
        tigs_loci = set(map(lambda x: x.split('.')[0], line.split()))
        if (tigs_loci.issubset(selection_loci)):
            selected_tigs.append(line)
            utilized_loci.update(tigs_loci)
    ip.close()

    print >> sys.stderr, "INFO: selected %d TIGs using %d of %d selection loci" % (len(selected_tigs), len(utilized_loci), len(selection_loci))

    return selected_tigs


def writeSelectedTIGs(selected_tigs, selected_tigs_output):
    op = open(selected_tigs_output, 'w')
    for tigs_line in selected_tigs:
        op.write("%s" % tigs_line)
    op.close()


if (__name__ == "__main__"):
    loci_file, tigs_file, selected_tigs_output = sys.argv[1:]

    loci = readSelectionLoci(loci_file)

    selected_tigs = selectTIGs(loci, tigs_file)

    writeSelectedTIGs(selected_tigs, selected_tigs_output)
    
    sys.exit(0)
