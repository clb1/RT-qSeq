#!/usr/bin/env python

from operator import itemgetter
from collections import defaultdict, Counter
from itertools import chain, groupby
from pyfaidx import Fasta
import cPickle
import sqlite3
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

sys.path.append("/raid1/projects/CGDB/scripts")
from classMFE import MFE, MFEOptions

from string import maketrans, translate
DNA_complement_table = maketrans("ACGTNacgtn","TGCANtgcan")

import designSetsClasses
import designParams

import pdb

def readPrimers(primers_tsv, ID_suffix):
    primers = {}

    ip = open(primers_tsv, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        primer_seq = fields[2].split('N')[-1] # TODO: split on designParams.UMI
        for full_pprID in set(fields[3].strip().split(';')):
            primer_ID = "%s%s" % (full_pprID, ID_suffix)
            assert (not primers.has_key(full_pprID)), "primers dictionary already has key %s" % full_pprID
            primers[full_pprID] = (primer_ID, primer_seq)

    ip.close()

    return primers

# 1) Process all primer pair and get all amplicons
# 2) Confirm global_TIGs and offtarget_igs recapitulated
#    Record mRNA start,stop position in each isoform for each primer pair-isoform combination
# 3) Read CGDB bed and grab bed line for all isoforms found in 1)
# 4) For each isoform, use strand info plus local mRNA start,stop to extract BED for primers product
# 
def getPrimerProducts(params_descriptor, transcriptomic_MFE, fwd_primers, rev_primers): # genome_ref, 
    primer_products = []
    amplicon_lens = []
    isoforms_w_product = set()
    
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    try:
        assert(set(fwd_primers.keys()) == set(rev_primers.keys())), "Mismatch between (chrom, OS_ID) sets for fwd & rev primers"
    except AssertionError as ae:
        pdb.set_trace()
    common_keys = set(fwd_primers.keys()) & set(rev_primers.keys())

    for full_pprID in common_keys: #fwd_primers.keys():
        print >> sys.stderr, full_pprID
        
        if (full_pprID.startswith('ERCC')):
            # TODO: handle later. Need to add amplicon to primer_products
            continue

        fwd_primer_ID, fwd_primer_seq = fwd_primers[full_pprID]
        rev_primer_ID, rev_primer_seq = rev_primers[full_pprID]        
        
        try:
            chrom, OS_ID, experiment_num, PPR_ID = full_pprID.split('_')
        except ValueError, ve:
            pdb.set_trace()
        select_tup = (chrom, OS_ID, params_descriptor, int(experiment_num), PPR_ID, fwd_primer_seq, rev_primer_seq)
        
        cursor.execute("SELECT local_target_isoform_group,global_target_isoform_group,offtarget_isoform_groups from PrimerPair \
        WHERE chrom=? AND OS_ID=? AND params_descriptor=? AND part_num=? AND PPR_ID=? and fwd_primer_seq=? AND rev_primer_seq=?", select_tup)

        results = cursor.fetchall()
        assert (len(results) == 1)
        local_TIGs = results[0][0].split(';')
        global_TIGs = results[0][1].split(';')
        offtarget_IGs = results[0][2].split(';')

    ampl_success = transcriptomic_MFE.processPrimer(fwd_primer_seq, rev_primer_seq)
    assert (ampl_success), "Primer amplification failure"
    products_per_isoform = transcriptomic_MFE.getPrimerProducts(designParams.is_directional)

        isoforms_per_product = defaultdict(list)
        for isoform_ID, products in products_per_isoform.items():
            isoforms_w_product.add(isoform_ID)
            if (len(products) > 1):
                print >> sys.stderr, "WARNING: multiple %s products for primer pair" % isoform_ID
                pdb.set_trace()
                
            # Note: local_start and local_stop are 0-based coordinates
            local_start, local_stop, nuc_seq, is_fwd_rev_primer_product = products[0][0:4]
            assert (is_fwd_rev_primer_product), "Product not based on fwd and rev primer combination"

            try:
                assert (any(map(lambda tig: isoform_ID in tig, global_TIGs)) or any(map(lambda ig: isoform_ID in ig, offtarget_IGs))
                        or len(nuc_seq) > designParams.confusion_amplicon_max_len), \
                        "Primers yield product from isoform %s, which is not in primer pair's global TIG or offtarget IG "
            except AssertionError, ae:
                pdb.set_trace()
                    
            isoforms_per_product[nuc_seq].append( (isoform_ID, local_start, local_stop) )

        for product, list_of_isoform_details in isoforms_per_product.items():
            isoform_IDs = map(itemgetter(0), list_of_isoform_details)
            seq_ID = "%s_%s_%s_%s|%s" % (chrom, OS_ID, experiment_num, PPR_ID, ",".join(isoform_IDs))
            amplicon_lens.append( len(product) )
            primer_products.append( (seq_ID, product, list_of_isoform_details) )

    conn.close()

    print >> sys.stderr, "INFO: %d amplicon lengths, %d products" % (len(amplicon_lens), len(primer_products))
    
    return amplicon_lens, primer_products, isoforms_w_product


def getIsoformBED(isoforms_w_product, CGDB_Bed_file):
    isoform_beds = {}
    ip = open(CGDB_Bed_file, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        if (fields[3] in isoforms_w_product):
            isoform_beds[fields[3]] = fields
    ip.close()
    return isoform_beds


def createProductBEDLines(primer_products, isoform_beds, genome_ref):
    product_beds = {}

    for seq_ID, product, list_of_isoform_details in primer_products:
        BED_for_product = None
        for isoform_ID, local_start, local_stop in list_of_isoform_details:
            isoform_bed = isoform_beds[isoform_ID]
            chromosome = isoform_bed[0]
            start = int(isoform_bed[1]) # BED, so 0-based
            strand = isoform_bed[5]
            block_sizes = map(int, isoform_bed[10].split(','))
            block_starts = map(int, isoform_bed[11].split(','))

            isoform_genomic_positions = [] # In 1-based coordinates
            for block_start, block_size in zip(block_starts,block_sizes):
                isoform_genomic_positions.extend( range(start+block_start+1, start+block_start+block_size+1) )

            if (strand == '-'):
                isoform_genomic_positions = isoform_genomic_positions[::-1]
                product_genomic_positions = isoform_genomic_positions[local_start:local_stop]
                product_genomic_positions = product_genomic_positions[::-1]
                #rev_local_start = len(isoform_genomic_positions) - local_stop - 1   # Results are 1-based. Remember,
                #rev_local_stop = len(isoform_genomic_positions) - local_start - 1   # local_start and local_stop are in 0-based coordinates.
                #product_genomic_positions = isoform_genomic_positions[rev_local_start-1:rev_local_stop]
            else:
                assert (strand == '+')
                product_genomic_positions = isoform_genomic_positions[local_start:local_stop]

            try:
                assert (len(product_genomic_positions) == len(product)), "Extracted number of genomic positions for product differs from product length"
                assert (len(isoform_genomic_positions) == sum(block_sizes)), "Different number of genomic and transcriptomic positions for isoform"
            except AssertionError, ae:
                print >> sys.stderr, ae.message
                pdb.set_trace()
                
            # Extract genomic positions of product and create the products BED
            product_start = product_genomic_positions[0]
            product_stop = product_genomic_positions[-1]
            product_nuc_seq = ""
            product_block_sizes = []
            product_block_starts = []
            nuc_seq_chunks = []
        for k, g in groupby(enumerate(product_genomic_positions), lambda (i,x):i-x):
        group = map(itemgetter(1), g)
        product_block_starts.append( group[0] - product_start )
        product_block_sizes.append( group[-1] - group[0] + 1 )
                nuc_seq_chunk = genome_ref[chromosome][group[0]-1:group[-1]] # For pyfaidx, starts are 1-based and ends are 0-based
                assert (len(nuc_seq_chunk) == product_block_sizes[-1]), "Nucleotide sequence of block doesn't match size of block"
                nuc_seq_chunks.append(nuc_seq_chunk)
        product_nuc_seq += nuc_seq_chunk
        product_nuc_seq = str(product_nuc_seq)
        
        if (strand == '-'):
        product_nuc_seq = product_nuc_seq[::-1].translate(DNA_complement_table)
            
            try:
                assert (product_nuc_seq == product), "Genome-extracted nucleotide sequence does not match amplicon nucleotide sequence"
            except AssertionError, ae:
                print >> sys.stderr, ae.message
                pdb.set_trace()
                
            bed_product_start = product_start-1
            bed_product_stop = bed_product_start + product_block_starts[-1] + product_block_sizes[-1]
            BED_tup = (chromosome, bed_product_start, bed_product_stop, strand, len(product_block_starts), ','.join(map(str,product_block_sizes)), ','.join(map(str,product_block_starts)))
            if (BED_for_product == None):
                BED_for_product = [BED_tup]
            else:
                BED_for_product.append(BED_tup)

        product_beds[seq_ID] = BED_for_product

    return product_beds
    

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

    plt.hist(amplicon_lens, bins=101, normed=False)
    #plt.gca().yaxis.set_major_formatter(formatter)    
    plt.title("%s, %s" % (params_descriptor, pool_label))
    plt.xlabel("Amplicon Length (bp)")
    plt.ylabel("Number of Amplicons")
    plt.savefig(fig_name)
    plt.gcf().clear()


def writePrimers(fwd_primers, rev_primers, output_primers_fasta):
    assert(set(fwd_primers.keys()) == set(rev_primers.keys())), "Mismatch between (chrom, OS_ID) sets for fwd & rev primers"
    op = open(output_primers_fasta, 'w')
    for full_pprID in fwd_primers.keys():
        fwd_primer_ID, fwd_primer_seq = fwd_primers[full_pprID]
        rev_primer_ID, rev_primer_seq = rev_primers[full_pprID]        
        op.write(">%s\n%s\n" % (fwd_primer_ID, fwd_primer_seq))
        op.write(">%s\n%s\n" % (rev_primer_ID, rev_primer_seq))
    op.close()


def writeProducts(primer_products, output_fasta_file):
    op = open(output_fasta_file, 'w')
    for seqID, seq in map(itemgetter(0,1), primer_products):
        op.write(">%s\n%s\n" % (seqID, seq))
    op.close()


def writeBed(product_beds, output_bed):
    itemRgb = "0,0,0"
    op = open(output_bed, 'w')
    for product_ID, bed_tuples in product_beds.items():
        product_ID = product_ID.replace(' ', '_')
        for chromosome, start, stop, strand, num_blocks, block_sizes, block_starts in bed_tuples:
        bed12_line = "%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s" % \
                 (chromosome, start, stop, product_ID, 1, strand, start, start, itemRgb, num_blocks, block_sizes, block_starts)
        op.write("%s\n" % bed12_line)
    op.close()


if (__name__ == "__main__"):
    params_descriptor, pool_label, database_path, genome_fasta, transcriptome_fasta, CGDB_bed_file, fwd_primers_tsv, rev_primers_tsv, output_basename = sys.argv[1:] 

    # Complete the output files' names
                                           
    output_primer_products_fasta = "%s_expected_products.fa" % output_basename
    output_primers_fasta = "%s.fa" % output_basename
    output_bed = "%s_expected_products.bed" % output_basename
    output_hist_fig = "%s_expected_products_lengths_hist.png" % output_basename

    designParams.setParameters(params_descriptor)
    genome_ref = Fasta(genome_fasta, as_raw=True, sequence_always_upper=True)

    mfe_options = MFEOptions()
    transcriptomic_MFE = MFE(transcriptome_fasta, mfe_options)

    fwd_primers = readPrimers(fwd_primers_tsv, "_5p")
    rev_primers = readPrimers(rev_primers_tsv, "_3p")
    amplicon_lens, primer_products, isoforms_w_products = getPrimerProducts(params_descriptor, transcriptomic_MFE, fwd_primers, rev_primers)
    isoform_beds = getIsoformBED(isoforms_w_products, CGDB_bed_file)

    #cPickle.dump((fwd_primers, rev_primers, amplicon_lens, primer_products, isoforms_w_products, isoform_beds), open("compile_debug.pkl", 'wb'))
    #(fwd_primers, rev_primers, amplicon_lens, primer_products, isoforms_w_products, isoform_beds) = cPickle.load(open("compile_debug.pkl", 'rb'))

    product_beds = createProductBEDLines(primer_products, isoform_beds, genome_ref)
    
    writeHistograms(amplicon_lens, params_descriptor, pool_label, output_hist_fig)
    writePrimers(fwd_primers, rev_primers, output_primers_fasta)
    writeProducts(primer_products, output_primer_products_fasta)
    writeBed(product_beds, output_bed)

    sys.exit(0)
