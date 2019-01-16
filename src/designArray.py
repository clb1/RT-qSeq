#!/usr/bin/env python

from collections import defaultdict
from itertools import chain
import numpy as np
from operator import itemgetter
import sqlite3
import sys
import designParams

from string import maketrans, translate
DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

import pdb

class ArrayPartition(object):
    def __init__(self, my_partition_number, required_fwd_primers, required_rev_primers, partition_rev_primer_revcomp):
        '''The required primers are lists of tuples (primer ID, primer sequence).'''
        self.my_partition_number = my_partition_number
        self.partition_rev_primer_ID = partition_rev_primer_revcomp[0]
        self.partition_rev_primer_revcomp = partition_rev_primer_revcomp[1]
        self.fwd_primers = defaultdict(list)
        self.rev_primers = defaultdict(list)

        for primerID, primer_seq in required_fwd_primers:
            self.fwd_primers[primer_seq].append(primerID)

        for primerID, primer_seq in required_rev_primers:
            self.rev_primers[primer_seq].append(primerID)


    def getPartitionNumber(self):
           return self.my_partition_number

    def getNumFwdPrimers(self):
        return len(self.fwd_primers.keys())


    def getNumRevPrimers(self):
        return len(self.rev_primers.keys())


    def getNumFwdAndRevPrimers(self):
        return len(self.fwd_primers.keys()) + len(self.rev_primers.keys())


    def numPrimersNew(self, fwd_primer_seqs, rev_primer_seqs):
        '''Returns the number of the primer sequences that are not
        already in the set of primer sequences for this partition.'''
        return sum(map(lambda x: not self.fwd_primers.has_key(x), fwd_primer_seqs)) + \
          sum(map(lambda y: not self.rev_primers.has_key(y), rev_primer_seqs))


    def addPrimers(self, fwd_primers, rev_primers):
        '''The primers are lists of tuples (primer ID, primer sequence).'''
        for primerID, primer_seq in fwd_primers:
            self.fwd_primers[primer_seq].append(primerID)

        for primerID, primer_seq in rev_primers:
            self.rev_primers[primer_seq].append(primerID)


    def getPartitionPrimers(self, fwd_adapter, rev_adapter, umi_len):
        '''Returns list of tuples (primerID, primer sequence). The 'primerID' is a composite of IDs 
        separated by a semicolon for cases where different IDs were associated with the same primer sequence'''
        sys.exit(1) # primer name needs to be more informative/specific. UMIs have now changed.
        fwd_adapter_plus_umi = fwd_adapter + 'N' * umi_len
        complete_fwd_primers = map(lambda x: ("P%d_Fwd_%s" % (self.my_partition_number, ";".join(x[1])),
                                              "%s%s%s" % (fwd_adapter_plus_umi, x[0], self.partition_rev_primer_revcomp)), self.fwd_primers.items())

        rev_adapter_plus_umi = rev_adapter + 'N' * umi_len
        complete_rev_primers = map(lambda x: ("P%d_Rev_%s" % (self.my_partition_number, ";".join(x[1])),
                                              "%s%s%s" % (rev_adapter_plus_umi, x[0], self.partition_rev_primer_revcomp)), self.rev_primers.items())

        return (complete_fwd_primers, complete_rev_primers)

        
class ArrayDesign(object):
    def __init__(self, database_file, partition_rev_primers_revcomp, required_fwd_primers, required_rev_primers, array_capacity):
        self.db = database_file
        self.num_partitions = len(partition_rev_primers_revcomp)
        self.array_partitions = []
        self.array_capacity = array_capacity
        self.prioritized_design_data = None
        for i in xrange(self.num_partitions):
            self.array_partitions.append( ArrayPartition(i+1, required_fwd_primers, required_rev_primers, partition_rev_primers_revcomp[i]) )


    def setPrioritizedOverlapSets(self, ranked_olap_sets):
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()

        self.prioritized_design_data = []
        for tup in ranked_olap_sets:
            chrom, OS_ID = tup[0:2]
            num_partitions = int(tup[6])
            cursor.execute("SELECT num_parts,ppr_combos FROM completedesigns WHERE chrom = ? AND OS_ID = ?", (chrom, OS_ID)) 
            row = cursor.fetchone()
            assert (row[0] == num_partitions)
            self.prioritized_design_data.append( (chrom, OS_ID, row[0], row[1]) )        

        conn.close()


    def prioritizeOverlapSets(self):
        sys.exit(1) # Depracated
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute("SELECT chrom,OS_ID,num_tigs,num_primerable_tigs,num_unwanted,num_parts,ppr_combos FROM completedesigns") #  WHERE timeout = Null
        all_rows = cursor.fetchall()
        conn.close()

        all_augmented_rows = []
        for row in all_rows:
            perc_primerable = round(100.0 * float(row[3])/float(row[2]), 1)
            if (row[3] > 0 and row[-1] != ' '): # TEMPORARY: row[-1]
                all_augmented_rows.append( row[0:4] + (perc_primerable,-row[4]) + row[5:] )

        all_aug_rows_sorted = sorted(all_augmented_rows, key=itemgetter(5,4,3), reverse=True) # lowest num_unwanted, then percent primerable, then number primerable

        self.prioritized_design_data = []
        for fields in all_aug_rows_sorted:
            self.prioritized_design_data.append( (fields[0], fields[1], fields[6], fields[7]) )


    def confirmPrimerDataPresent(self):
        print >> sys.stderr, "INFO: confirming the presence of all potentially-needed primers..."
        problems_found = False
        assert (self.prioritized_design_data != None and len(self.prioritized_design_data) > 0)
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()

        for chrom, OS_ID, olap_set_num_parts, ppr_IDs_all_partitions in self.prioritized_design_data:
            ppr_IDs_per_partition = map(lambda x: x.split(','), ppr_IDs_all_partitions.split(';'))
            for ppr_ID in chain.from_iterable(ppr_IDs_per_partition):
                cursor.execute("SELECT fwd_primer_seq,rev_primer_seq FROM primers WHERE chrom = ? AND OS_ID = ? AND ppr_ID = ?", (chrom, OS_ID, ppr_ID))
                row = cursor.fetchone()
                if (row == None or len(row) != 2):
                    print >> sys.stderr, "ERROR: primers problem for %s PPR %s on %s" % (OS_ID, ppr_ID, chrom)
                    problems_found = True
                else:
                    assert (row[0] != None and len(row[0]) > 10)
                    assert (row[1] != None and len(row[1]) > 10)

        conn.close()

        return problems_found


    def fillPartitions(self):
        print >> sys.stderr, "INFO: Initiating design of array"
        assert (self.prioritized_design_data != None and len(self.prioritized_design_data) > 0)

        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()

        num_primers_per_partition = np.zeros((self.num_partitions), 'i')
        for pnum in xrange(self.num_partitions):
            num_primers_per_partition[pnum] = self.array_partitions[pnum].getNumFwdAndRevPrimers()
        remaining_array_capacity = self.array_capacity - sum(num_primers_per_partition)
        assert (remaining_array_capacity >= 0)

        curr_olap_set = 0
        while (remaining_array_capacity > 1):
            chrom, OS_ID, olap_set_num_parts, ppr_IDs_all_partitions = self.prioritized_design_data[curr_olap_set]
            ppr_IDs_per_partition = map(lambda x: x.split(','), ppr_IDs_all_partitions.split(';'))

            # Set the number of usable experiment partitionings that can be used from the current overlap set
            max_num_parts = min(olap_set_num_parts, self.num_partitions)
            
            ranked_array_partitions = map(itemgetter(1), sorted(zip(num_primers_per_partition, range(self.num_partitions)), key=itemgetter(0)))

            # Add primers for up to max_num_parts partitions
            olap_set_part_index = 0
            num_parts_added = 0
            while (num_parts_added < max_num_parts and olap_set_part_index < olap_set_num_parts and remaining_array_capacity > 1):
                pnum = ranked_array_partitions.pop(0) # Select the array partition containing the next fewest total number of primers

                # Get the primers for the next overlap set partition and determine how
                # many don't already exist in the selected array partition
                new_fwd_primers, new_rev_primers = [], []
                for ppr_ID in ppr_IDs_per_partition[olap_set_part_index]:
                    cursor.execute("SELECT fwd_primer_seq,rev_primer_seq FROM primers WHERE chrom = ? AND OS_ID = ? AND ppr_ID = ?",
                                   (chrom, OS_ID, params_descriptor, ppr_ID))
                    row = cursor.fetchone()
                    new_fwd_primers.append( (ppr_ID, row[0]) )
                    new_rev_primers.append( (ppr_ID, row[1]) )
                olap_set_part_index += 1

                # Add the new primer sequences if doing so doesn't exceed the array capacity
                new_fwd_primer_seqs = set(map(itemgetter(1), new_fwd_primers))
                new_rev_primer_seqs = set(map(itemgetter(1), new_rev_primers))
                num_new_primer_seqs = self.array_partitions[pnum].numPrimersNew(new_fwd_primer_seqs, new_rev_primer_seqs)
                if (num_new_primer_seqs <= remaining_array_capacity):
                    self.array_partitions[pnum].addPrimers(new_fwd_primers, new_rev_primers)
                    num_primers_per_partition[pnum] = self.array_partitions[pnum].getNumFwdAndRevPrimers()
                    remaining_array_capacity = self.array_capacity - sum(num_primers_per_partition)
                    assert (remaining_array_capacity >= 0)
                    num_parts_added += 1
                
            curr_olap_set += 1
        conn.close()

        # Print report
        all_parts_total_primers = 0
        print >> sys.stderr, "\nPartition\tNum_Fwd\tNum_Rev\tTotal"
        print >> sys.stderr, "---------\t-------\t-------\t-----"
        for i in xrange(self.num_partitions):
            part_num = self.array_partitions[i].getPartitionNumber()
            num_fwd_primers = self.array_partitions[i].getNumFwdPrimers()
            num_rev_primers = self.array_partitions[i].getNumRevPrimers()
            total_num_primers = num_fwd_primers + num_rev_primers
            all_parts_total_primers += total_num_primers
            print >> sys.stderr, "\t%d\t%d\t%d\t%d" % (part_num, num_fwd_primers, num_rev_primers, total_num_primers)

        print >> sys.stderr, "Array total = %d, Array capacity = %d, Unused capacity = %d\n" % \
          (all_parts_total_primers, self.array_capacity, self.array_capacity - all_parts_total_primers)

        #assert (sum(primers_per_partition) == self.array_capacity), "Array not exactly filled"

        
    def writeFasta(self, P5_sequencing_adapter, P7_sequencing_adapter, output_fasta):
        op = open(output_fasta, 'w')
        for partition in self.array_partitions:
            fwd_primers, rev_primers = partition.getPartitionPrimers(P5_sequencing_adapter, P7_sequencing_adapter, 8)
            for ID, seq in fwd_primers:
                op.write(">%s\n%s\n" % (ID, seq))
            for ID, seq in rev_primers:
                op.write(">%s\n%s\n" % (ID, seq))

        op.close()


def selectPartitionReversePrimers(num_partitions_needed, partition_rev_primers_fasta):
    partition_rev_primers_revcomp = []

    ip = open(partition_rev_primers_fasta,'r')
    line = ip.readline()

    while (len(partition_rev_primers_revcomp) < num_partitions_needed):
        assert(line[0] == '>')
        rev_primer_ID = line[1:].strip()
        primer_seq = ip.readline().strip()
        primer_seq_revcomp = primer_seq[::-1].translate(DNA_complement_table)
        partition_rev_primers_revcomp.append( (rev_primer_ID, primer_seq_revcomp) )
        line = ip.readline()

    ip.close()


    return partition_rev_primers_revcomp


def readRankedOlapSets(ranked_olap_sets_tsv):
    ranked_olap_sets = []

    ip = open(ranked_olap_sets_tsv, 'r')
    for line in ip:
        if (line[0] != '#'):
            fields = line.strip().split("\t")
            ranked_olap_sets.append( fields )
    ip.close()

    return ranked_olap_sets


def readRequiredPrimers(required_primers_tsv):
    required_fwd_primers = []
    required_rev_primers = []

    ip = open(required_primers_tsv, 'r')
    for line in ip:
        if (line[0] != '#'):
            ID, fwd_primer_seq, rev_primer_seq, product = line.strip().split("\t")
            required_fwd_primers.append( (ID, fwd_primer_seq) )
            required_rev_primers.append( (ID, rev_primer_seq) )

    ip.close()

    return (required_fwd_primers, required_rev_primers)


if (__name__ == "__main__"):
    # partition_rev_primers_fasta, 
    params_descriptor, main_db, num_partitions, ranked_olap_sets_tsv, required_primers_tsv, output_fasta = sys.argv[1:]
    num_partitions = int(num_partitions)

    designParams.setParameters(params_descriptor)
    P5_sequencing_adapter = designParams.fwd_primer_adapter_tag # Correct association?
    P7_sequencing_adapter = designParams.rev_primer_adapter_tag

    ranked_olap_sets = readRankedOlapSets(ranked_olap_sets_tsv)
    partition_rev_primers_revcomp = selectPartitionReversePrimers(num_partitions, partition_rev_primers_fasta)
    required_fwd_primers, required_rev_primers = readRequiredPrimers(required_primers_tsv)

    array_design = ArrayDesign(main_db, partition_rev_primers_revcomp, required_fwd_primers, required_rev_primers, 12472)
    #array_design.prioritizeOverlapSets()
    array_design.setPrioritizedOverlapSets(ranked_olap_sets)

    problems_found = array_design.confirmPrimerDataPresent()
    if (problems_found):
        print >> sys.stderr, "ERROR: problems found with primer sequence data. Exiting without designing array."
    else:
        array_design.fillPartitions()
        array_design.writeFasta(P5_sequencing_adapter, P7_sequencing_adapter, output_fasta)

    sys.exit(0)
    
