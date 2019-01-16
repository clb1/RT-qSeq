#!/usr/bin/env python

import cPickle
from operator import itemgetter
import sqlite3
import sys

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range

import pdb


def getOverlapSetIsoformsAndPrimerSets(database_path, chromosome, OS_ID):
    chrom_OS_tuple = (chromosome, OS_ID)

    conn = sqlite3.connect(database_path, isolation_level=None)
    cursor = conn.cursor()

    # Confirm that a design has been registered for the overlap set
    cursor.execute("SELECT count(*) FROM completedesigns WHERE chrom = ? AND OS_ID = ?", chrom_OS_tuple)
    num_entries = cursor.fetchone()[0]
    assert (num_entries == 1)

    cursor.execute("SELECT predesign_pkl_file FROM predesigns WHERE chrom = ? AND OS_ID = ?", chrom_OS_tuple)
    all_rows = cursor.fetchall()
    assert (len(all_rows)!=0), "Database has no entry for %s %s" % chrom_OS_tuple
    assert (len(all_rows)==1), "Database has multiple entries for %s %s" % chrom_OS_tuple
    predesign_pkl_file = all_rows[0][0]
    olap_set = cPickle.load(open(predesign_pkl_file, 'rb'))
    
    cursor.execute("SELECT ppr_combos FROM completedesigns WHERE chrom = ? AND OS_ID = ?", chrom_OS_tuple)
    row = cursor.fetchone()
    pprs_per_experiment = row[0].split(";")

    ppr_primer_data = {}
    cursor.execute("SELECT ppr_ID,fwd_primer_seq,fwd_genomic_positions,rev_primer_seq,rev_genomic_positions FROM primers WHERE chrom = ? AND OS_ID = ?", chrom_OS_tuple)
    all_rows = cursor.fetchall()
    for ppr_ID, fwd_primer_seq, fwd_primer_genomic_pos_str, rev_primer_seq, rev_primer_genomic_pos_str in all_rows:
        fwd_primer_genomic_pos = set(map(int, fwd_primer_genomic_pos_str.split(' ')))
        rev_primer_genomic_pos = set(map(int, rev_primer_genomic_pos_str.split(' ')))
        ppr_primer_data[ppr_ID] = (fwd_primer_seq, fwd_primer_genomic_pos, rev_primer_seq, rev_primer_genomic_pos)

    conn.close()

    return (olap_set, pprs_per_experiment, ppr_primer_data)


def compressIntrons(olap_set):
    olap_set.setIsoformsCoordCorrespondData()
        
    # For now, just create an uncompressed mapping
    region_start, region_stop = olap_set.getRegionStartStop()
    region_len = region_stop - region_start + 1
    faux_genome_positions = dict(zip(range(region_start, region_stop+1), range(1,region_len+1)))

    return faux_genome_positions


def addPrimeringAnnotation(faux_genome_positions, olap_set, pprs_per_experiment, ppr_primer_data, style, output_dir):
    fgp = faux_genome_positions

    all_isoforms = olap_set.getAllOverlappingIsoforms()

    pdb.set_trace()

    # Create a figure for each experiment.
    for experiment_num, ppr_IDs in enumerate(pprs_per_experiment,1):
        gt_mRNAs_to_draw = []

        # Add primer to each isoform
        for isoform in all_isoforms:
            seqid = isoform.getChromosome()
            isoform_ID = isoform.getCGDBName()
            isoform_strand = isoform.getStrand()
            isoform_start, isoform_stop = isoform.getStartStop()
            gt_mRNA = FeatureNode.create_new(seqid, "mRNA", fgp[isoform_start], fgp[isoform_stop], isoform_strand)
            gt_mRNA.set_source(isoform_ID)

            print >> sys.stderr, isoform_ID

            # Add the introns. Exons are added below
            for intron_start, intron_stop in isoform.getIntronTuples():
                gt_intron = FeatureNode.create_new(seqid, "intron", fgp[intron_start], fgp[intron_stop], isoform_strand)
                gt_mRNA.add_child(gt_intron)

            fwd_symbol, rev_symbol = ('+','-') if (isoform_strand == '+') else ('-','+')
            genomic_positions_per_exon = isoform.getExonsGenomicPositions()

            primer_starts = set()
            for ppr_ID in ppr_IDs.split(','):
                print >> sys.stderr, "\t%s" % ppr_ID
                fwd_primer, fwd_primer_genomic_pos, rev_primer, rev_primer_genomic_pos = ppr_primer_data[ppr_ID]
                if (isoform_strand == '+'):
                    print >> sys.stderr, "\t\tprimer region: %d-%d" % (min(fwd_primer_genomic_pos), max(rev_primer_genomic_pos))
                    if (isoform.hasPositions(fwd_primer_genomic_pos) and isoform.areAllSequential(fwd_primer_genomic_pos)):
                        fwd_primer_start_pos = isoform.getIsoformCoord(min(fwd_primer_genomic_pos))
                        primer_starts.add( ("Fwd", fwd_primer_start_pos) )
                    if (isoform.hasPositions(rev_primer_genomic_pos) and isoform.areAllSequential(rev_primer_genomic_pos)):
                        rev_primer_start_pos = isoform.getIsoformCoord(max(rev_primer_genomic_pos))
                        primer_starts.add( ("Rev", rev_primer_start_pos) )
                else:
                    print >> sys.stderr, "\t\tprimer region: %d-%d" % (min(rev_primer_genomic_pos),max(fwd_primer_genomic_pos))
                    if (isoform.hasPositions(fwd_primer_genomic_pos) and isoform.areAllSequential(fwd_primer_genomic_pos)):
                        fwd_primer_start_pos = isoform.getIsoformCoord(max(fwd_primer_genomic_pos))
                        primer_starts.add( ("Fwd", fwd_primer_start_pos) )
                    if (isoform.hasPositions(rev_primer_genomic_pos) and isoform.areAllSequential(rev_primer_genomic_pos)):
                        rev_primer_start_pos = isoform.getIsoformCoord(min(rev_primer_genomic_pos))
                        primer_starts.add( ("Rev", rev_primer_start_pos) )

            primer_starts = sorted(primer_starts, key=itemgetter(0))
            #if (len(primer_starts) > 0):
            #    pdb.set_trace()

            # "Chunk" the isoform into ranges defined by singleton Rev primers or Fwd/Rev primer pairs
            only_rev = []
            fwd_and_rev = []
            if (len(primer_starts) > 1):
                li = 0 if (primer_starts[0][0]=="Rev") else primer_starts[0][1]
                for i, (primer_type, primer_start) in enumerate(primer_starts):
                    if (primer_type == "Rev"):
                        if (i==0 or primer_starts[i-1][0] == "Rev"):
                            only_rev.append( (li, primer_start) )
                        else:
                            fwd_and_rev.append( (li, primer_start) )
                        li = primer_start+1
                    else:
                        li = primer_start

            # Determine the "state" of each position in the isoform. Form the ranges for each separate state "chunk"
            rev_chunks = []
            fwd_chunks = []
            for li, ri in only_rev:
                rev_chunk_genomic_coords = isoform.convertToGenomicCoords( range(li,ri+1) )
                rev_chunks.append( set(rev_chunk_genomic_coords) )
            for li, ri in fwd_and_rev:
                sep_dist = ri - li + 1
                fwd_chunk_genomic_coords = isoform.convertToGenomicCoords( range(li, min(li+150,ri+1)) )
                fwd_chunks.append( set(fwd_chunk_genomic_coords) )
                rev_chunk_genomic_coords = isoform.convertToGenomicCoords( range(max(ri-150, li), ri+1) )
                rev_chunks.append( set(rev_chunk_genomic_coords) )

            # For each exon, annotate the portions of it that are in different states
            for exon_positions in genomic_positions_per_exon:
                exon_start = min(exon_positions)
                exon_stop = max(exon_positions)
                gt_exon = FeatureNode.create_new(seqid, "exon", fgp[exon_start], fgp[exon_stop], isoform_strand)

                for fwd_chunk in fwd_chunks:
                    fwd_positions = fwd_chunk & exon_positions
                    if (len(fwd_positions) > 0):
                        gt_fwd_primer = FeatureNode.create_new(seqid, "fwd_primer", fgp[min(fwd_positions)], fgp[max(fwd_positions)], fwd_symbol)
                        gt_exon.add_child(gt_fwd_primer)

                for rev_chunk in rev_chunks:
                    rev_positions = rev_chunk & exon_positions
                    if (len(rev_positions) > 0):
                        gt_rev_primer = FeatureNode.create_new(seqid, "rev_primer", fgp[min(rev_positions)], fgp[max(rev_positions)], rev_symbol)
                        gt_exon.add_child(gt_rev_primer)

                gt_mRNA.add_child(gt_exon)

            gt_mRNAs_to_draw.append(gt_mRNA)

        # Draw
        range_max = max(fgp.values())
        png_file = "%s/experiment_%d.png" % (output_dir, experiment_num)
        diagram = Diagram.from_array(gt_mRNAs_to_draw, Range(1, range_max), style)
        layout = Layout(diagram, range_max, style)
        height = layout.get_height()
        canvas = CanvasCairoFile(style, range_max, height)
        layout.sketch(canvas)
        canvas.to_file(png_file)


    olap_set.unsetIsoformsCoordCorrespondData()


if (__name__ == "__main__"):
    database_path, style_file, chromosome, OS_ID, output_dir = sys.argv[1:]

    style = Style()
    style.load_file(style_file)

    olap_set, pprs_per_experiment, ppr_primer_data = getOverlapSetIsoformsAndPrimerSets(database_path, chromosome, OS_ID)
    faux_genome_positions = compressIntrons(olap_set)

    addPrimeringAnnotation(faux_genome_positions, olap_set, pprs_per_experiment, ppr_primer_data, style, output_dir)

    sys.exit(0)
