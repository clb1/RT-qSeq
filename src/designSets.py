#!/usr/bin/env python3
from collections import defaultdict
import datetime
from operator import itemgetter, methodcaller, or_
from itertools import chain, groupby
import networkx as nx
import gc
import stat
import sys
import os
from subprocess import check_call, CalledProcessError

import sqlite3

import pickle
from pyfaidx import Fasta
import time
    
import designParams
from designSetsClasses import *
import pdb


def setChromosomes(chromosomes_of_interest):
    chromosomes = None
    all_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
                       "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"] 

    if (chromosomes_of_interest == "all"):
        chromosomes = all_chromosomes
    else:
        assert (all(map(lambda x: x in all_chromosomes, chromosomes_of_interest.split(','))))
        chromosomes = chromosomes_of_interest.split(',')

    return chromosomes


def readTargetIsoformIDGroups(target_isoform_groups_file):
    '''Returns a list of isoform group tuples (ie sorted isoform IDs)'''
    target_isoform_ID_groups = []
    with open(target_isoform_groups_file, 'r') as ip:
        for line in ip:
            line = line.strip()
            if (line != '' and line[0] != '#'):
                isoform_IDs = line.split()
                isoform_group = sorted(isoform_IDs)
                target_isoform_ID_groups.append( tuple(isoform_group) )
    return target_isoform_ID_groups


def readSubsumedIsoformData(subsumed_isoforms_tsv):
    isoform_subsumed_by = {}
    with open(subsumed_isoforms_tsv, 'r') as ip:
        header = ip.readline()
        for line in ip:
            subsumed_isoform, subsuming_isoforms = line.strip().split("\t")
            isoform_subsumed_by[subsumed_isoform] = set(subsuming_isoforms.split())
    return isoform_subsumed_by


def readCompleteDesignTargets(chrom_OS_pairs_file):
    chr_OS_tuples = []
    ip = open(chrom_OS_pairs_file, 'r')
    for line in ip:
        if (line[0] != '#'):
            chrom, olap_set_ID = line.strip().split("\t")
            chr_OS_tuples.append( (chrom, olap_set_ID) )
    ip.close()
    return chr_OS_tuples


def createOverlapSets(chrom, target_isoform_ID_groups, CGDB_bed, tempdir, pkl_olap_sets_dir, database_path):
    all_isoforms = {}
    G = nx.Graph()

    # Create BED file and RNAIsoform objects for target isoforms
    target_isoforms_bed = "%s/%s_target_isoforms.bed"  % (tempdir, chrom)
    op_bed = open(target_isoforms_bed, 'w')
    ip_bed = open(CGDB_bed, 'r')

    flattened_target_isoform_ID_groups = set(chain.from_iterable(target_isoform_ID_groups))

    for line in ip_bed:
        bed_fields = line.strip().split("\t")
        isoform_ID = bed_fields[3]
        if (bed_fields[0] == chrom and isoform_ID in flattened_target_isoform_ID_groups): # len(bed_fields) > 0
            op_bed.write("%s" % line)
            isoform = RNAIsoform2(bed_fields[3], bed_fields[0], bed_fields[5])
            block_sizes = list(map(int, bed_fields[10].split(',')))
            block_starts = list(map(int, bed_fields[11].split(',')))
            isoform.addExons(bed_fields[0], int(bed_fields[1]), block_sizes, block_starts, bed_fields[5])
            isoform.setLengthAndExonTermini()
            isoform.setAsTarget()
            all_isoforms[bed_fields[3]] = isoform

    op_bed.close()
    ip_bed.close()
        
    print("INFO: intersecting target isoforms with all isoforms...", file=sys.stderr, flush=True, end='')
    intersect_output = "%s/intersect_output.tsv" % tempdir
    intersect_cmd = ['intersectBed', '-a', target_isoforms_bed, '-b', CGDB_bed, '-wo', '-split', '-s', '>', intersect_output]
    intersect_cmd_str = " ".join(intersect_cmd)
    try:
        return_code = check_call(intersect_cmd_str, shell=True)
        assert(return_code==0), "Return code was %d" % return_code
        print("done", file=sys.stderr, flush=True)
    except CalledProcessError as cpe:
        pdb.set_trace()
    except AssertionError as ae:
        pdb.set_trace()

    ip_intersect = open(intersect_output, 'r')
    for line in ip_intersect:
        bed_fields = line.strip().split('\t')
        isoform1_ID = bed_fields[3]
        isoform2_ID = bed_fields[15]
            
        if (isoform1_ID != isoform2_ID):
            if (isoform2_ID not in all_isoforms):
                assert (isoform1_ID in all_isoforms)
                isoform = RNAIsoform2(bed_fields[15], bed_fields[12], bed_fields[17])
                block_sizes = list(map(int, bed_fields[22].split(',')))
                block_starts = list(map(int, bed_fields[23].split(',')))
                isoform.addExons(bed_fields[12], int(bed_fields[13]), block_sizes, block_starts, bed_fields[17])
                isoform.setLengthAndExonTermini()
                isoform.setAsNontarget()
                all_isoforms[isoform2_ID] = isoform

            G.add_edge(all_isoforms[isoform1_ID], all_isoforms[isoform2_ID])

    # Cleanup
    ip_intersect.close()
    os.remove(intersect_output)
    os.remove(target_isoforms_bed)

    all_target_isoform_groups = []
    IDs_of_other_chromosome_TIG_members = {}
    for isoform_ID_group in target_isoform_ID_groups:
        if (any(map(lambda isoform_ID: isoform_ID in all_isoforms, isoform_ID_group))):
            if (all(map(lambda isoform_ID: isoform_ID in all_isoforms, isoform_ID_group))):
                isoform_group = list(map(lambda isoform_ID: all_isoforms[isoform_ID], isoform_ID_group))
                isoform_group = sorted(isoform_group, key=methodcaller("getCGDBName"))
                tig = tuple(isoform_group)
                all_target_isoform_groups.append( tig )
            else:
                # Some isoform groups are composed of non-overlapping mRNAs, possibly from different chromosomes, that are
                # grouped together because they encode sequence-identical pORFs. This step removes those isoforms in the
                # current group that are from a different chromosome.
                # IMPORTANT: What is not removed are isoforms in the isoform group that are on the same chromosome
                # but that will be in different overlapping sets of isoforms. Those cases are handled after the overlapping
                # sets of isoforms have been identified.
                same_chrom_isoform_ID_group = filter(lambda isoform_ID: isoform_ID in all_isoforms, isoform_ID_group)
                isoform_group = list(map(lambda isoform_ID: all_isoforms[isoform_ID], same_chrom_isoform_ID_group))
                isoform_group = sorted(isoform_group, key=methodcaller("getCGDBName"))
                ig_tuple = tuple(isoform_group)
                all_target_isoform_groups.append( ig_tuple )

                # Record those TIGs that are composed of identical isoforms from another chromosome, and the IDs of those isoforms.
                IDs_of_other_chrom_isoforms = list(filter(lambda isoform_ID: isoform_ID not in all_isoforms, isoform_ID_group))
                if (len(IDs_of_other_chrom_isoforms)>0):
                    assert (ig_tuple not in IDs_of_other_chromosome_TIG_members), "Should not have already encountered isoform group"
                    IDs_of_other_chromosome_TIG_members[ig_tuple] = IDs_of_other_chrom_isoforms

    # Identify the disjoinst sets of overlapping isoforms and the associated target isoform groups.
    # For historical reasons, if there are target isoforms groups within an overlap set that derive from both strands, the logic below
    # handle the isoforms from each strand separately. This now should not happen because 'intersect' above is now parameterized with
    # "same strand" set to 'True'.
    unnamed_isoform_olap_sets = []
    counter = 1

    for overlapping_isoforms in nx.connected_components(G):
        # Set the target isoform groups that pertain to this overlapping set of isoforms
        plus_strand_target_isoform_groups = set()
        minus_strand_target_isoform_groups = set()

        # These two map partial TIGs that are used to form in olap_sets to full TIG comprised of identical isoforms from outside the olap_set.
        plus_strand_complete_TIG_IDs_w_nonlocal_members = {}
        minus_strand_complete_TIG_IDs_w_nonlocal_members = {}

        for target_isoform_group in all_target_isoform_groups:
            if (set(target_isoform_group).issubset(overlapping_isoforms)):
                target_isoform_group_IDs = list(map(lambda x:x.getCGDBName(), target_isoform_group))
                target_isoform_group_IDs.sort()
                target_isoform_group_IDs = tuple(target_isoform_group_IDs)
                
                strands = set(map(lambda i: i.getStrand(), target_isoform_group))
                assert(len(strands)==1)
                group_strand = list(strands)[0]

                if (group_strand == "+"):
                    plus_strand_target_isoform_groups.add( target_isoform_group )
                    if (target_isoform_group in IDs_of_other_chromosome_TIG_members):
                        # For the TIG, record its actual membership that is composed of identical isoforms from elsewhere in the genome
                        IDs = list(map(lambda i:i.getCGDBName(), target_isoform_group)) + IDs_of_other_chromosome_TIG_members[target_isoform_group]
                        plus_strand_complete_TIG_IDs_w_nonlocal_members[target_isoform_group_IDs] = tuple(sorted(IDs))
                else:
                    minus_strand_target_isoform_groups.add( target_isoform_group )
                    if (target_isoform_group in IDs_of_other_chromosome_TIG_members):
                        # For the TIG, record its actual membership that is composed of identical isoforms from elsewhere in the genome
                        IDs = list(map(lambda i:i.getCGDBName(), target_isoform_group)) + IDs_of_other_chromosome_TIG_members[target_isoform_group]
                        minus_strand_complete_TIG_IDs_w_nonlocal_members[target_isoform_group_IDs] = tuple(sorted(IDs))

            elif (not set(target_isoform_group).isdisjoint(overlapping_isoforms)):
                assert (len(set(map(lambda i: i.getCGDBName().split('.')[0], target_isoform_group))) > 1)
                # Since mRNA isoforms have been grouped based on pORF sequences, non-overlapping isoform can be grouped. Such should be the case here.
                # Remove from the target_isoform_group those mRNA isoforms with identical pORFs that are from elsewhere on the same chromosome.
                this_olap_set_target_isoform_group = tuple(filter(lambda y: y in overlapping_isoforms, target_isoform_group))
                assert (len(this_olap_set_target_isoform_group) > 0)

                # For the this_olap_set TIG, record its actual membership that is composed of identical isoforms from elsewhere in the genome
                # (on same and/or different chromosome)
                IDs = list(map(lambda i: i.getCGDBName(), target_isoform_group))               # Includes all same-chromosome identical isoforms
                if (target_isoform_group in IDs_of_other_chromosome_TIG_members):
                    IDs += IDs_of_other_chromosome_TIG_members[target_isoform_group]     # Includes all not-same-chromosome identical chromosomes
                IDs.sort()

                strands = set(map(lambda i: i.getStrand(), this_olap_set_target_isoform_group))
                assert(len(strands)==1)
                group_strand = list(strands)[0]

                this_olap_set_target_isoform_group_IDs = list(map(lambda i:i.getCGDBName(), this_olap_set_target_isoform_group))
                this_olap_set_target_isoform_group_IDs.sort()
                this_olap_set_target_isoform_group_IDs = tuple(this_olap_set_target_isoform_group_IDs)

                if (group_strand == "+"):
                    plus_strand_target_isoform_groups.add( tuple(this_olap_set_target_isoform_group) )
                    plus_strand_complete_TIG_IDs_w_nonlocal_members[this_olap_set_target_isoform_group_IDs] = tuple(IDs)
                else:
                    minus_strand_target_isoform_groups.add( tuple(this_olap_set_target_isoform_group) )
                    minus_strand_complete_TIG_IDs_w_nonlocal_members[this_olap_set_target_isoform_group_IDs] = tuple(IDs)


        if (len(plus_strand_target_isoform_groups) > 1):
            #TIGs_w_subsuming_isoforms_IDs = addSubsumingIsoformsIDs(plus_strand_target_isoform_groups, isoform_subsumed_by, all_isoforms)
            olap_set = OverlappingIsoforms("NOT SET", overlapping_isoforms, plus_strand_target_isoform_groups, plus_strand_complete_TIG_IDs_w_nonlocal_members) # TIGs_w_subsuming_isoforms_IDs, 
            unnamed_isoform_olap_sets.append( (olap_set, olap_set.getRegionStartStop()) )
            counter += 1
                
        if (len(minus_strand_target_isoform_groups) > 1):
            #TIGs_w_subsuming_isoforms_IDs = addSubsumingIsoformsIDs(minus_strand_target_isoform_groups, isoform_subsumed_by, all_isoforms)
            olap_set = OverlappingIsoforms("NOT SET", overlapping_isoforms, minus_strand_target_isoform_groups, minus_strand_complete_TIG_IDs_w_nonlocal_members) # TIGs_w_subsuming_isoforms_IDs, 
            unnamed_isoform_olap_sets.append( (olap_set, olap_set.getRegionStartStop()) )
            counter += 1


    # Name the overlapping sets in a canonical and reproducible way
    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    unnamed_isoform_olap_sets.sort(key=itemgetter(1))
    for counter, isoform_olap_set in enumerate(map(itemgetter(0), unnamed_isoform_olap_sets), 1):
        OS_ID = "OS%d" % counter
        strand = isoform_olap_set.getStrand()
        region_start, region_stop = isoform_olap_set.getRegionStartStop()
        isoform_IDs = isoform_olap_set.getAllOverlappingIsoformsIDs()
        isoform_IDs_str = " ".join(isoform_IDs)
        isoform_olap_set.setID(OS_ID)
        pkl_file = "%s/%s_%s_olap_set.pkl" % (pkl_olap_sets_dir, chrom, OS_ID)
        cursor.execute('INSERT INTO IsoformOlapSets VALUES (?,?,?,?,?,?,?)', (chrom, OS_ID, strand, region_start, region_stop, isoform_IDs_str, pkl_file))
        print("INFO: Writing %s..." % pkl_file, file=sys.stderr, flush=True, end='')
        pickle.dump(isoform_olap_set, open(pkl_file, "wb"))
        print("done", file=sys.stderr, flush=True)

    conn.close()


def addSubsumingIsoformsIDs(target_isoform_groups, isoform_subsumed_by, all_isoforms):
    TIGs_w_subsuming_isoforms_IDs = []
    for tig in target_isoform_groups:
        tig_IDs = list(map(lambda i:i.getCGDBName(), tig))
        tig_w_subsuming_IDs = set(tig_IDs)
        for isoform_ID in filter(lambda i: i in isoform_subsumed_by, tig_IDs):
            tig_w_subsuming_IDs.update( isoform_subsumed_by[isoform_ID] )
        tig_w_subsuming_IDs = sorted(tig_w_subsuming_IDs)
        TIGs_w_subsuming_isoforms_IDs.append( tig_w_subsuming_IDs )
    return TIGs_w_subsuming_isoforms_IDs


def partitionOlapSets(database_path, num_partitions, olap_sets_to_partition_tsv, partition_files_basename):
    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    olap_sets_to_partition = set()
    with open(olap_sets_to_partition_tsv, 'r') as ip:
        for line in ip:
            chrom, OS_ID = line.strip().split("\t")
            olap_sets_to_partition.add( (chrom,OS_ID) )

    cursor.execute("SELECT chrom,OS_ID FROM IsoformOlapSets")
    all_db_olap_sets = set(cursor.fetchall())

    assert (olap_sets_to_partition.issubset(all_db_olap_sets)), "All target olap sets not in database"

    ordering = list(map(str, range(23))) + ['X','Y','M']
    all_olap_sets_annot = list(map(lambda x: (x[0], x[1], ordering.index(str(x[0][3:])), int(x[1][2:])), olap_sets_to_partition))
    all_olap_sets_annot.sort(key=itemgetter(2,3))
    ordered_olap_sets = list(map(itemgetter(0,1), all_olap_sets_annot))

    for part_num, L in enumerate(np.array_split(ordered_olap_sets, num_partitions), 1):
        op = open("%s_%d.txt" % (partition_files_basename, part_num), 'w')
        for (chrom, OS_ID) in L:
            op.write("%s\t%s\n" % (chrom, OS_ID))
        op.close()

    conn.close()


def instantiateDatabase(database_path):
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    cursor.execute('''CREATE TABLE IF NOT EXISTS IsoformOlapSets
    (chrom text NOT NULL,
    OS_ID text NOT NULL,
    strand text NOT NULL,
    genomic_start int NOT NULL,
    genomic_stop int NOT NULL,
    isoform_IDs text NOT NULL,
    olap_set_pkl_file text NOT NULL,
    PRIMARY KEY (chrom, OS_ID))
    ''') 

    cursor.execute('''CREATE TABLE IF NOT EXISTS Predesigns
    (chrom text NOT NULL,
    OS_ID text NOT NULL,
    params_descriptor text NOT NULL,
    num_target_isoform_groups integer NOT NULL,
    num_unprimerable_target_isoform_groups integer NOT NULL,
    num_primerable_target_isoform_groups integer NOT NULL,
    percent_primerable_isoform_groups float NOT NULL,
    primerable_target_isoform_groups text,
    unprimerable_target_isoform_groups text,
    num_pprs_w_primer_pairs integer NOT NULL,
    num_pprs integer NOT NULL,
    num_target_isoforms integer NOT NULL,
    num_overlapping_isoforms integer NOT NULL,
    predesign_pkl_file text NOT NULL,
    stderr_output text NOT NULL,
    PRIMARY KEY (chrom, OS_ID, params_descriptor))
    ''') 

    cursor.execute('''CREATE TABLE IF NOT EXISTS ExperimentDesignMeta
    (chrom text NOT NULL,
    OS_ID text NOT NULL,
    params_descriptor text NOT NULL,
    num_tigs int NOT NULL,
    num_primerable_tigs int NOT NULL,
    num_unwanted int NOT NULL,
    num_offtarget int NOT NULL,
    num_redundant int NOT NULL,
    expected_num_unwanted_from_primers float NOT NULL,
    num_parts int NOT NULL,
    timeout int,
    solns_pkl_file text NOT NULL,
    PRIMARY KEY (chrom, OS_ID, params_descriptor))
    ''') 

    cursor.execute('''CREATE TABLE IF NOT EXISTS ExperimentDesignPart
    (chrom text NOT NULL,
    OS_ID text NOT NULL,
    params_descriptor text NOT NULL,
    part_num int NOT NULL,
    ppr_ID_combo text NOT NULL,
    tigs text NOT NULL,
    PRIMARY KEY (chrom, OS_ID, params_descriptor, part_num))
    ''') 

    cursor.execute('''CREATE TABLE IF NOT EXISTS PrimerPair
    (chrom text NOT NULL,
    OS_ID text NOT NULL,
    params_descriptor text NOT NULL,
    part_num int NOT NULL,
    ppr_ID text NOT NULL,
    fwd_primer_seq text NOT NULL,
    target_fwd_genomic_positions text NOT NULL,
    rev_primer_seq text NOT NULL,
    target_rev_genomic_positions text NOT NULL,
    local_target_isoform_group text NOT NULL,
    global_target_isoform_group text,
    offtarget_isoform_groups text,
    primer_products_confirmed int DEFAULT 0,
    PRIMARY KEY (chrom, OS_ID, params_descriptor, ppr_ID))
    ''') 

    conn.close()


def getOlapSetsSpecs(olap_sets_to_predesign, database_path):
    olap_sets_specs = []

    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    with open(olap_sets_to_predesign, 'r') as ip:
        for line in ip:
            chrom, OS_ID = line.strip().split('\t')            
            cursor.execute("SELECT olap_set_pkl_file FROM IsoformOlapSets WHERE chrom = ? AND OS_ID = ?", (chrom, OS_ID))
            olap_set_pkl_file = cursor.fetchone()[0]
            assert (os.path.exists(olap_set_pkl_file)), "%s does not exist" % olap_set_pkl_file
            olap_sets_specs.append( (chrom, OS_ID, olap_set_pkl_file) )

    conn.close()
    
    return olap_sets_specs


def predesignPrimerSets(genome_ref, olap_sets_specs, database_path, pkl_predesigns_dir, tempdir, transcriptome_vsearch_udb, transcriptome_ref, antitargets_fasta): # transcriptomic_MFE
    oligo_thermo = OligoThermodynamics()

    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    for (chrom, OS_ID, olap_set_pkl_file) in olap_sets_specs:
        olap_set = pickle.load(open(olap_set_pkl_file,'rb'))
        olap_set_ID = olap_set.getID()
        olap_region_start, olap_region_stop = olap_set.getRegionStartStop()
        predesign_pkl_file = "%s/%s_%s_predesign.pkl" % (pkl_predesigns_dir, chrom, olap_set.getID())

        # If the predesign has already been completed for this overlap set, continue
        cursor.execute("SELECT predesign_pkl_file FROM Predesigns WHERE chrom = ? AND OS_ID = ? and params_descriptor = ?", \
                       (chrom, olap_set_ID, designParams.params_descriptor))
        all_rows = cursor.fetchall()
        if (len(all_rows) >= 1):
            assert (len(all_rows)==1), "Database has multiple entries for %s %s %s" % (chrom, olap_set_ID, designParams.params_descriptor)
            continue
        else:
            assert (not os.path.exists(predesign_pkl_file)), "%s already exists" % predesign_pkl_file

        stderr_output = []
        olap_set_descr_tuple = olap_set.getSizeDescription()
        
        stderr_output.append("\n### %s, %s:%d-%d ###" % (olap_set_ID, chrom, olap_region_start, olap_region_stop))
        print(stderr_output[-1], file=sys.stderr, flush=True)
        stderr_output.append("%d TIGs composed of %d target isoforms amongst %d overlapping isoforms" % olap_set_descr_tuple)
        print(stderr_output[-1], file=sys.stderr, flush=True)

        print("\tSetup sequence: 1 ", file=sys.stderr, end='', flush=True)
        olap_set.setIsoformsCoordCorrespondData()
        olap_set.setIsoformSequences(genome_ref)     # TODO: This could be made consistent with how same thing
        olap_set.setIsoformsSeqs(transcriptome_ref)  # is done in doCompleteDesign() below

        print("> 2 ", file=sys.stderr, end='', flush=True)
        olap_set.setPrimerGroupsPerPosition()

        print("> 3 ", file=sys.stderr, end='', flush=True)
        olap_set.buildEquivPrimerTargetsRegions(genome_ref)

        print("> 4 ", file=sys.stderr, end='', flush=True)
        olap_set.pairPrimerRegions(genome_ref)

        print("> 5 ", file=sys.stderr, end='', flush=True)
        num_pprs_w_primers, num_pprs = olap_set.designDefaultPrimersForAllPairedPrimerRegions(oligo_thermo, tempdir, transcriptome_vsearch_udb, transcriptome_ref, antitargets_fasta)
        primer_design_msg = "Could design default primers for %d of %d paired primer regions" % (num_pprs_w_primers, num_pprs)
        stderr_output.append(primer_design_msg)

        print("> 6 ", file=sys.stderr, end='', flush=True)
        design_limitations_msg, num_unprimerable_TIGs, num_TIGs, \
            primerable_isoform_groups_str, unprimerable_isoform_groups_str = olap_set.setPrimerableTargetIsoformGroups()
        num_primerable_TIGs = num_TIGs - num_unprimerable_TIGs
        perc_primerable_TIGs = round(100.0 * float(num_primerable_TIGs)/float(num_TIGs), 1)

        olap_set.unsetIsoformsCoordCorrespondData()
        olap_set.unsetIsoformsSeqs()

        print("> Done\n\t%s\n" % primer_design_msg, file=sys.stderr, flush=True)
        print(design_limitations_msg, file=sys.stderr, flush=True)
        stderr_output.append(design_limitations_msg)

        num_isoforms_in_TIGs, num_isoforms_in_olap_set = olap_set_descr_tuple[1:3]

        # TODO: add information that would allow selectAndWriteOlapSetsSubset() to not need to query stderr_output for TIG information
        predesigns_tuple = (chrom, olap_set_ID, designParams.params_descriptor, num_TIGs, num_unprimerable_TIGs, num_primerable_TIGs, perc_primerable_TIGs,
                            primerable_isoform_groups_str, unprimerable_isoform_groups_str,
                            num_pprs_w_primers, num_pprs, num_isoforms_in_TIGs, num_isoforms_in_olap_set, predesign_pkl_file, "\n".join(stderr_output))

        print("\nINFO: writing %s" % predesign_pkl_file, file=sys.stderr, flush=True)
        pickle.dump(olap_set, open(predesign_pkl_file, 'wb'))
        os.chmod(predesign_pkl_file, stat.S_IRUSR | stat.S_IRGRP)
        cursor.execute('INSERT INTO Predesigns VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', predesigns_tuple)

        olap_set = None
        num_collected = gc.collect()
        
    conn.close()


def selectAndWriteOlapSetsSubset(param_settings_descriptor, target_isoform_ID_groups, output_tsv):
    olap_sets_subset_annot = set()
    all_target_isoform_IDs = set(chain.from_iterable(target_isoform_ID_groups))

    ordering = list(map(str, range(1,23))) + ['X','Y','M']
    
    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()
    
    cursor.execute('SELECT chrom,OS_ID,isoform_IDs FROM IsoformOlapSets')
    for chrom, OS_ID, isoform_IDs_str in cursor.fetchall():
        isoform_IDs = set(isoform_IDs_str.split())
        if (len(isoform_IDs & all_target_isoform_IDs) > 0):
            olap_sets_subset_annot.add( (chrom, OS_ID, ordering.index(chrom[3:]), int(OS_ID[2:])) )
    conn.close()
    
    #print >> sys.stderr, "INFO: selected %d (with %d TIGs) of %d overlap sets" % (len(olap_sets_subset), total_num_tigs, counter)

    olap_sets_subset = sorted(olap_sets_subset_annot, key=itemgetter(2,3))
    op = open(output_tsv, 'w');
    for tup in olap_sets_subset:
        op.write("%s\t%s\n" % (tup[0], tup[1]))
    op.close()


def doCompleteDesign(chrom_OS_tuples, database_path, transcriptome_ref, transcriptome_vsearch_udb, antitargets_fasta_name, tempdir, solns_dir): # transcriptomic_MFE, 
    allowed_solve_time = 60 # In minutes
    
    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    for chrom_OS_tuple in chrom_OS_tuples:
        chrom, olap_set_ID = chrom_OS_tuple

        # Check whether designs for overlap set have already been completed
        selection_tuple = chrom_OS_tuple + (designParams.params_descriptor,)
        cursor.execute("SELECT count(*) FROM ExperimentDesignMeta WHERE chrom = ? AND OS_ID = ? and params_descriptor = ?", selection_tuple)
        num_entries = cursor.fetchone()[0]
        if (num_entries == 1):
            print("%s, %s, %s has already been computed. Skipping." % selection_tuple, file=sys.stderr, flush=True)
            continue
        else:
            assert (num_entries == 0)

        cursor.execute("SELECT predesign_pkl_file FROM Predesigns WHERE chrom = ? AND OS_ID = ? and params_descriptor = ?", selection_tuple)
        all_rows = cursor.fetchall()
        try:
            assert (len(all_rows)!=0), "Database has no entry for %s %s %s" % selection_tuple
            assert (len(all_rows)==1), "Database has multiple entries for %s %s %s" % selection_tuple
        except AssertionError as ae:
            print(str(ae), file=sys.stderr, flush=True)
            continue
        predesign_pkl_file = all_rows[0][0]
        olap_set = pickle.load(open(predesign_pkl_file, 'rb'))
        olap_set_descr_tuple = olap_set.getSizeDescription()

        if (False): # Debug BED file
            op = open("%s/%s_equiv_regions.bed" % (tempdir,olap_set.getID()), 'w')
            op.write("track type=bed description=%s name=%s itemRgb=On visibility=full\n" % (olap_set.getID(), olap_set.getID()))
            ##op_details = open("%s_equiv_regions.txt" % olap_set.getID(), 'w')
            olap_set.writeEquivalenceRegions(op, True, True) # op_details, 
            ##op_details.close()
            op.close()

        olap_region_start, olap_region_stop = olap_set.getRegionStartStop()
        print("\n### %s, %s:%d-%d ###" % (olap_set_ID, chrom, olap_region_start, olap_region_stop), file=sys.stderr, flush=True)
        print("%d TIGs composed of %d target isoforms amongst %d overlapping isoforms" % olap_set_descr_tuple, file=sys.stderr, flush=True)
        
        if (False):
            op = open("%s_equiv_regions.bed" % olap_set.getID(), 'w')
            op.write("track type=bed description=%s name=%s itemRgb=On visibility=full\n" % (olap_set.getID(), olap_set.getID()))
            op_details = open("%s_equiv_regions.txt" % olap_set.getID(), 'w')
            olap_set.writeEquivalenceRegions(op, op_details, True, True)
            op_details.close()
            op.close()

        olap_set.setIsoformsCoordCorrespondData()
        olap_set.setIsoformsSeqs(transcriptome_ref)

        begin_time = datetime.datetime.now()
        print("BEGIN: %s" % begin_time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr, flush=True)

        try:
            C = PPRClassCombos(olap_set, transcriptome_ref, transcriptome_vsearch_udb, antitargets_fasta_name) # transcriptomic_MFE, 
            C.setTemporaryDirectory(tempdir)
            C.setLocalTranscriptomeFasta()
            C.setExpireTime(begin_time, allowed_solve_time)
            ppr_combos_solution, timed_out, results_str = C.findSolution()
            timeout_value = allowed_solve_time if (timed_out) else None # TODO: Delete?
        except AssertionError as ae:
            if (str(ae) == "No primerable target isoform groups"):
                results_str = str(ae)
                ppr_combos_solution = None
                timeout_value = None
            else:
                pdb.set_trace()
                raise ae

        C.unsetLocalTranscriptomeFasta()
        olap_set.unsetIsoformsCoordCorrespondData()
        olap_set.unsetIsoformsSeqs()

        print(results_str, file=sys.stderr, flush=True)
        end_time = datetime.datetime.now()
        print("END: %s" % end_time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr, flush=True)
        time_delta = end_time - begin_time
        seconds = time_delta.total_seconds() 

        num_wanted_TIGs = olap_set.getNumTargetIsoformGroups("all")
        num_primerable_TIGs = olap_set.getNumTargetIsoformGroups("primerable")

        olap_set = None
        num_collected = gc.collect()

        continue
    
        # Record solution details
        if (ppr_combos_solution != None):
            # ppr_combos_solution[-1] is an Npart list of (ppr_combo, ppr_combo_primer_sets), where ppr_combo_primer_sets is a list of
            # lists of primer sets that solve for ppr_combo.

            solns_pkl_file = "%s/%s_%s.pkl" % (solns_dir, chrom, olap_set_ID)

            ppr_combo_IDs_and_pp_sets = []
            for ppr_combo, pp_sets in ppr_combos_solution[-1]:
                ppr_combo_IDs_and_pp_sets.append( (list(map(lambda ppr: ppr.getID(), ppr_combo)), pp_sets) )

            db_tuple = chrom_OS_tuple + (designParams.params_descriptor, num_wanted_TIGs, num_primerable_TIGs) + ppr_combos_solution[0:-1] + (timeout_value, solns_pkl_file)
            cursor.execute('INSERT INTO ExperimentDesignMeta VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', db_tuple)
            pickle.dump(ppr_combo_IDs_and_pp_sets, open(solns_pkl_file, 'wb'))
            os.chmod(solns_pkl_file, stat.S_IRUSR | stat.S_IRGRP)

        else:
            null_tuple = chrom_OS_tuple + (designParams.params_descriptor, num_wanted_TIGs, num_primerable_TIGs) + (-1, -1, -1, 0, 0, None, '')
            cursor.execute('INSERT INTO ExperimentDesignMeta VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', null_tuple)

    conn.close()


def writeGFF3(all_orig_isoforms, output_gff3): # olap_sets, CGDB_bed
    # This method probably needs to be based on reading all of CGDB and handling each isoform
    # individually in terms of whether it is written out and with whatever particular annotation (e.g. color)
    mRNA_tuples = set()
    parent_mRNAs = defaultdict(set)
    primer_region_parent_exons = defaultdict(set)
    amplicon_parent_exons = defaultdict(set)
    for isoforms in [all_orig_isoforms.values()]: # unique_isoforms
        for isoform in isoforms:
            parent_mRNA = isoform.getCGDBName()
            mRNA_tuple, exon_tuples, primer_region_tuples, amplicon_tuples = isoform.getAsGFF3(None)
            mRNA_tuples.add(mRNA_tuple)

            for exon_tuple in exon_tuples:
                parent_mRNAs[exon_tuple].add(parent_mRNA)

            for primer_region_tuple in primer_region_tuples:
                parent_exon = primer_region_tuple[-1]
                primer_region_tuple = tuple(primer_region_tuple[0:-1])
                primer_region_parent_exons[primer_region_tuple].add(parent_exon)

            for amplicon_tuple in amplicon_tuples:
                parent_exon = amplicon_tuple[-1]
                amplicon_tuple = tuple(amplicon_tuple[0:-1])
                amplicon_parent_exons[amplicon_tuple].add(parent_exon)

    # http://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
    #RGB_colors = ["240,163,255","0,117,220","153,63,0","76,0,92","0,92,49","43,206,72","255,204,153","128,128,128",
    #           "148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16",
    #            "94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5"]
    hex_colors = ["#FFFFFF", "#F0A3FF", "#0075DC", "#993F00", "#4C005C", "#005C31", "#2BCE48", "#FFCC99", "#808080"]

    op = open(output_gff3, 'w')
    op.write("##gff-version 3\n")
    op.write("##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606\n")

    for mRNA_tuple in mRNA_tuples:
        op.write("%s\tCGDB\tmRNA\t%d\t%d\t0\t%s\t.\t%s\n" % mRNA_tuple)
            
    exon_tuple_IDs = {}
    for counter, exon_tuple in enumerate(parent_mRNAs.keys(), 1):
        exon_ID = "e%d" % counter
        exon_tuple_IDs[exon_tuple] = exon_ID
        parents = ",".join(parent_mRNAs[exon_tuple])
        new_exon_tuple = tuple(list(exon_tuple) + ["ID=%s;Parent=%s;" % (exon_ID, parents)])
        op.write("%s\tCGDB\texon\t%d\t%d\t0\t%s\t.\t%s\n" % new_exon_tuple)

    #counter = 0
    #denom = len(hex_colors) - 1
    # color_index = counter%denom + 1
    for primer_region_tuple, exon_parents in primer_region_parent_exons.items():
        #parents = ",".join(list(map(lambda x: exon_tuple_IDs[x], exon_parents)))
        parents = reduce(or_, map(lambda x: parent_mRNAs[x], exon_parents))
        parents_str = ",".join(parents)
        chromosome, start, stop, strand, primer_region_ID = primer_region_tuple
        color_index = 1 if (primer_region_ID.startswith('F')) else 2
        tags = "ID=%s;Parent=%s;" % (primer_region_ID, parents_str)
        op.write("%s\tCGDB\tprimer_region\t%d\t%d\t%d\t%s\t.\t%s\n" % (chromosome, start, stop, color_index, strand, tags))
            
    sep_amplicon = defaultdict(list)
    for amplicon_tuple, exon_parents in amplicon_parent_exons.items():
        parents = reduce(or_, map(lambda x: parent_mRNAs[x], exon_parents))
        parents_str = ",".join(parents)
        chromosome, start, stop, score_val, strand, amplicon_part_ID = amplicon_tuple
        tags = "ID=%s;Parent=%s;" % (amplicon_part_ID, parents_str)
        op.write("%s\tCGDB\tamplicon\t%d\t%d\t%d\t%s\t.\t%s\n" % (chromosome, start, stop, score_val, strand, tags))

        # For displaying amplicons as fake mRNA
        amplicon_root_ID = amplicon_part_ID.split('.')[0]
        tags = "ID=sep%s;Parent=%s;" % (amplicon_part_ID, amplicon_root_ID)
        sep_amplicon[amplicon_root_ID].append( (chromosome, start, stop, score_val, strand, tags) )

    # Write the amplicons as fake mRNA
    for amplicon_root_ID, part_tuples in sep_amplicon.items():
        amplicon_start = min(list(map(itemgetter(1), part_tuples)))
        amplicon_stop = max(list(map(itemgetter(2), part_tuples)))
        score_vals = set(list(map(itemgetter(3), part_tuples)))
        assert(len(score_vals)==1)
        score_val = list(score_vals)[0]
        op.write("%s\tCGDB\tmRNA\t%d\t%d\t%d\t%s\t.\tID=%s;\n" % (chromosome, amplicon_start, amplicon_stop, score_val, strand, amplicon_root_ID))
        for part_tuple in part_tuples:
            op.write("%s\tCGDB\texon\t%d\t%d\t%d\t%s\t.\t%s\n" % part_tuple)
    op.close()


if (__name__ == "__main__"):
    action, database_path = sys.argv[1:3]
    assert (action in ["create_olap_sets", "olap_sets_subset", "partition_olap_sets", "predesign", "completedesign"]), "Unknown action '%s' specified" % action

    instantiateDatabase(database_path)

    if (action == "create_olap_sets"):
        tempdir_root, target_isoform_groups_file, CGDB_bed, pkl_olap_sets_dir = sys.argv[3:] # subsumed_isoforms_tsv, 

        designParams.setParameters("dir_default")

        tempdir = "%s/designSets_%d" % (tempdir_root, os.getpid())
        os.mkdir(tempdir)
        print("Temporary directory is %s" % tempdir, file=sys.stderr, flush=True)

        chromosomes = setChromosomes("all")
        target_isoform_ID_groups = readTargetIsoformIDGroups(target_isoform_groups_file)
        #isoform_subsumed_by = readSubsumedIsoformData(subsumed_isoforms_tsv)

        for chrom in chromosomes:
            print("INFO: creating overlap set for %s" % chrom, file=sys.stderr, flush=True)
            createOverlapSets(chrom, target_isoform_ID_groups, CGDB_bed, tempdir, pkl_olap_sets_dir, database_path)

        os.rmdir(tempdir)


    elif (action == "olap_sets_subset"):
        param_settings_descriptor, target_isoform_groups_file, output_tsv = sys.argv[2:]
        target_isoform_ID_groups = readTargetIsoformIDGroups(target_isoform_groups_file)        
        print("INFO: read %d target isoforms in %d groups" % \
                (len(list(chain.from_iterable(target_isoform_ID_groups))), len(target_isoform_ID_groups)), file=sys.stderr, flush=True)
        selectAndWriteOlapSetsSubset(param_settings_descriptor, target_isoform_ID_groups, output_tsv)


    elif (action == "partition_olap_sets"):
        num_partitions, olap_sets_to_partition_tsv, partition_files_basename = sys.argv[3:]
        partitionOlapSets(database_path, int(num_partitions), olap_sets_to_partition_tsv, partition_files_basename)

        
    elif (action == "predesign"):
        param_settings_descriptor, tempdir_root, CGDB_bed, genome_fasta, transcriptome_fasta, transcriptome_vsearch_udb, \
        olap_sets_to_predesign, antitargets, pkl_predesigns_dir = sys.argv[3:]

        designParams.setParameters(param_settings_descriptor)

        genome_ref = Fasta(genome_fasta, as_raw=True, sequence_always_upper=False)
        transcriptome_ref = Fasta(transcriptome_fasta, as_raw=True, sequence_always_upper=False)

        tempdir = "%s/designSets_%d" % (tempdir_root, os.getpid())
        os.mkdir(tempdir)
        print("Temporary directory is %s" % tempdir, file=sys.stderr, flush=True)

        antitargets_fasta = tempfile.NamedTemporaryFile(mode='wt', suffix=".fa", dir=tempdir, delete=True)
        with open(antitargets, 'r') as ip:
            for line in ip:
                isoform_ID = line.strip()
                if (isoform_ID != ''):
                    antitargets_fasta.write(">%s\n%s\n" % (isoform_ID, str(transcriptome_ref[isoform_ID])))
        antitargets_fasta.flush()
        try:
            makeblastdb_output = check_output(['makeblastdb', '-dbtype', 'nucl', '-in', antitargets_fasta.name], stderr=subprocess.STDOUT)
        except CalledProcessError as cpe:
            pdb.set_trace()

        olap_sets_specs = getOlapSetsSpecs(olap_sets_to_predesign, database_path)

        predesignPrimerSets(genome_ref, olap_sets_specs, database_path, pkl_predesigns_dir, tempdir,
                            transcriptome_vsearch_udb, transcriptome_ref, antitargets_fasta.name) # transcriptomic_MFE

        antitargets_fasta.close()
        os.rmdir(tempdir)


    elif (action == "completedesign"):
        param_settings_descriptor, tempdir_root, chrom_OS_pairs_file, antitargets, transcriptome_fasta, transcriptome_vsearch_udb, solns_dir = sys.argv[3:]

        designParams.setParameters(param_settings_descriptor)

        tempdir = "%s/designSets_%d" % (tempdir_root, os.getpid())
        os.mkdir(tempdir)
        print("Temporary directory is %s" % tempdir, file=sys.stderr, flush=True)

        #mfe_options = MFEOptions()
        #transcriptomic_MFE = MFE(transcriptome_fasta, mfe_options)
        transcriptome_ref = Fasta(transcriptome_fasta, as_raw=True, sequence_always_upper=True)

        antitargets_fasta = tempfile.NamedTemporaryFile(mode='wt', suffix=".fa", dir=tempdir, delete=True)
        with open(antitargets, 'r') as ip:
            for line in ip:
                isoform_ID = line.strip()
                if (isoform_ID != ''):
                    antitargets_fasta.write(">%s\n%s\n" % (isoform_ID, str(transcriptome_ref[isoform_ID])))
        antitargets_fasta.flush()
        try:
            check_call(['makeblastdb', '-dbtype', 'nucl', '-in', antitargets_fasta.name], stderr=subprocess.STDOUT)
        except CalledProcessError as cpe:
            pdb.set_trace()

        chrom_OS_tuples = readCompleteDesignTargets(chrom_OS_pairs_file)
        doCompleteDesign(chrom_OS_tuples, database_path, transcriptome_ref, transcriptome_vsearch_udb, antitargets_fasta.name, tempdir, solns_dir) # transcriptomic_MFE, 

        antitargets_fasta.close()
        os.rmdir(tempdir)

    sys.exit(0)
