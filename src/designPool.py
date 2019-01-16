#!/usr/bin/env python

from collections import defaultdict
from itertools import chain
import numpy as np
from operator import itemgetter
from scipy.stats import rankdata
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

        
# Rank olaps sets by num_tigs, etc criteria. Then select the first Npart set of primers that have not alread been selected for a previously
# designed pool.
class PoolDesign(object):
    def __init__(self, database_file, params_descriptor, max_pool_size):
        self.db = database_file
        self.params_descriptor = params_descriptor
        self.fwd_primers = defaultdict(list)
        self.rev_primers = defaultdict(list)
        self.max_pool_size = max_pool_size

        self.ranked_TIGs = None
        self.no_use_olap_sets = None
        self.ranked_partitions = None
        self.ranked_primer_sets = None


    def readRequiredPrimers(self, required_primers_tsv):
        required_fwd_primers = []
        required_rev_primers = []

        ip = open(required_primers_tsv, 'r')
        for line in ip:
            if (line[0] != '#'):
                ID, fwd_primer_seq, rev_primer_seq, product = line.strip().split("\t")
                required_fwd_primers.append( (ID, fwd_primer_seq) )
                required_rev_primers.append( (ID, rev_primer_seq) )
        ip.close()

        for primerID, primer_seq in required_fwd_primers:
            self.fwd_primers[primer_seq].append(primerID)

        for primerID, primer_seq in required_rev_primers:
            self.rev_primers[primer_seq].append(primerID)


    def readRankedTIGs(self, ranked_TIGs_file):
        self.ranked_TIGs = []
        ip = open(ranked_TIGs_file, 'r')
        for line in ip:
            fields = line.strip().split("\t")
            self.ranked_TIGs.append( fields[0].split() )
        ip.close()


    def readOverlapSetsToAvoid(self, no_use_olap_sets_tsv):
        self.no_use_olap_sets = set()
        with open(no_use_olap_sets_tsv, 'r') as ip:
            for line in ip:
                fields = line.strip().split()
                self.no_use_olap_sets.add( (fields[0],fields[1],int(fields[2])) )


    def rankOverlapSetPartitions(self):
        self.ranked_partitions = []
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()

        num_found_partition = 0
        for TIG in self.ranked_TIGs:
            partition = set()
            for isoform in TIG:
                cursor.execute("SELECT chrom,OS_ID,part_num FROM ExperimentDesignPart WHERE params_descriptor = ? AND tigs LIKE ?",
                               (self.params_descriptor,'%'+isoform+'%'))
                entries = cursor.fetchall()

                if (len(entries) == 1):
                    partition.update(entries)

            if (len(partition) == 1):
                num_found_partition += 1
                TIG_str = ",".join(TIG)
                self.ranked_partitions.append( list(partition)[0]+(TIG_str,) )
        conn.close()

        print >> sys.stderr, "INFO: found partition for %d of %d ranked TIGs" % (num_found_partition, len(self.ranked_TIGs))


    def rankPrimerSets(self):
        self.ranked_primer_sets = []

        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        for chrom,OS_ID,partition_num in map(itemgetter(0,1,2), self.ranked_partitions):
            if ((chrom,OS_ID,partition_num) not in self.no_use_olap_sets):
                partition_fwd_primers, partition_rev_primers = [], []
                cursor.execute("SELECT ppr_ID,fwd_primer_seq,rev_primer_seq FROM PrimerPair where chrom = ? AND OS_ID = ? and params_descriptor = ? AND part_num = ?",
                               (chrom, OS_ID, self.params_descriptor, partition_num))

                for ppr_ID, fwd_primer_seq, rev_primer_seq in cursor.fetchall():
                    partition_fwd_primers.append( (ppr_ID, fwd_primer_seq) )
                    partition_rev_primers.append( (ppr_ID, rev_primer_seq) )

                self.ranked_primer_sets.append( (chrom, OS_ID, partition_num, partition_fwd_primers, partition_rev_primers) )

        conn.close()


    def fillPool(self):
        assert (len(self.ranked_primer_sets) > 0), "No ranked designs for which to compile primer sets"

        print >> sys.stderr, "INFO: Initiating design of pool"

        fwd_primer_pool = set(self.fwd_primers.keys())
        rev_primer_pool = set(self.rev_primers.keys())

        remaining_rev_pool_capacity = self.max_pool_size - len(rev_primer_pool)
        remaining_fwd_pool_capacity = self.max_pool_size - len(fwd_primer_pool)
        
        assert (remaining_rev_pool_capacity > 0 and remaining_fwd_pool_capacity > 0)

        curr_olap_set_index = 0
        while (remaining_rev_pool_capacity > 0 and remaining_fwd_pool_capacity > 0 and curr_olap_set_index < len(self.ranked_primer_sets)):
            # Get next overlap set
            chrom, OS_ID, partition_num, partition_fwd_primers, partition_rev_primers = self.ranked_primer_sets[curr_olap_set_index]

            # Get a primer set that doesn't exceed capacity
            new_fwd_pool_size = len(fwd_primer_pool | set(map(itemgetter(1), partition_fwd_primers)))
            new_rev_pool_size = len(rev_primer_pool | set(map(itemgetter(1), partition_rev_primers)))

            pool_capacities_not_exceeded = new_fwd_pool_size <= self.max_pool_size and new_rev_pool_size <= self.max_pool_size
            if (pool_capacities_not_exceeded):
                for ppr_ID, fwd_primer_seq in partition_fwd_primers:
                    fwd_primer_pool.add(fwd_primer_seq)
                    self.fwd_primers[fwd_primer_seq].append( (chrom, OS_ID, str(partition_num), ppr_ID) )

                for ppr_ID, rev_primer_seq in partition_rev_primers:
                    rev_primer_pool.add(rev_primer_seq)
                    self.rev_primers[rev_primer_seq].append( (chrom, OS_ID, str(partition_num), ppr_ID) )

            curr_olap_set_index += 1
            remaining_rev_pool_capacity = self.max_pool_size - len(rev_primer_pool)
            remaining_fwd_pool_capacity = self.max_pool_size - len(fwd_primer_pool)

        # Print report


    def convertNumberTo96Well(self, number):
        return 'ABCDEFGH'[(number-1)%8] + "%d" % ((number-1)/8 + 1,)


    def writeFasta(self, output_fasta_basename):
        # Write the reverse primers pool
        L = []
        for rev_primer_seq, tups in self.rev_primers.items():
            merged_id = ";".join( map(lambda x: x if isinstance(x,str) else "_".join(x), tups) )
            L.append( (merged_id, rev_primer_seq) )

        L.sort(key=itemgetter(0))

        op_primers = open("%s_rev.tsv" % output_fasta_basename, 'w')
        counter = 1
        for merged_id, rev_primer_seq in L:
            well = self.convertNumberTo96Well(counter)
            simple_name = "rev%d" % counter
            tagged_rev_primer_seq = "%s%s%s" % (designParams.rev_primer_adapter_tag, "NNNNNN", rev_primer_seq)
            op_primers.write("%s\t%s\t%s\t%s\n" % (well, simple_name, tagged_rev_primer_seq, merged_id))
            counter += 1
        op_primers.close()


        # Write the forward primers pool
        L = []
        for fwd_primer_seq, tups in self.fwd_primers.items():
            merged_id = ";".join( map(lambda x: x if isinstance(x,str) else "_".join(x), tups) )
            L.append( (merged_id, fwd_primer_seq) )

        L.sort(key=itemgetter(0))

        op_primers = open("%s_fwd.tsv" % output_fasta_basename, 'w')
        counter = 1
        for merged_id, fwd_primer_seq in L:
            well = self.convertNumberTo96Well(counter)
            simple_name = "fwd%d" % counter
            tagged_fwd_primer_seq = "%s%s%s" % (designParams.fwd_primer_adapter_tag, "NNNNNN", fwd_primer_seq)
            op_primers.write("%s\t%s\t%s\t%s\n" % (well, simple_name, tagged_fwd_primer_seq, merged_id))
            counter += 1
        op_primers.close()



if (__name__ == "__main__"):
    params_descriptor, main_db, num_primer_pairs, ranked_TIGs_file, no_use_olap_sets_tsv, required_primers_tsv, output_fasta_basename = sys.argv[1:]

    designParams.setParameters(params_descriptor)

    pool_design = PoolDesign(main_db, params_descriptor, int(num_primer_pairs))
    pool_design.readRequiredPrimers(required_primers_tsv)
    pool_design.readRankedTIGs(ranked_TIGs_file)
    pool_design.readOverlapSetsToAvoid(no_use_olap_sets_tsv)
    pool_design.rankOverlapSetPartitions()
    pool_design.rankPrimerSets()

    pdb.set_trace()
    pool_design.fillPool()
    pool_design.writeFasta(output_fasta_basename)

    sys.exit(0)
    
