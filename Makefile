SHELL = /bin/bash

TARGET_SET = PM
CGDB_VERSION = v1.0

ROOT_DIR = /raid1/projects/CoveringPrimersSets
SCRIPTS = ${ROOT_DIR}/scripts
CGDB_BED = /raid1/projects/CGDB/models/CGDB${CGDB_VERSION}.bed
SCRATCH_DIR = /raid1/projects/scratch
GENOME_FASTA = /raid1/references_and_indexes/hg38/hg38.fa

TARGET_SET_ROOT_DIR = ${ROOT_DIR}/target_sets/${TARGET_SET}

OLAP_SETS_PKL_DIR = ${TARGET_SET_ROOT_DIR}/pkl_olap_sets/${CGDB_VERSION}
PREDESIGN_PKL_DIR = ${TARGET_SET_ROOT_DIR}/pkl_predesign/${CGDB_VERSION}

P5_ADAPTER_TAG = CGACGCTCTTCCGATCT
P7_ADAPTER_TAG = GTGTGCTCTTCCGATC

%_covering_set_primers.gff3 %_covering_set_primers.tsv: target_isoforms-%.txt ${CGDB_BED}
    ${SCRIPTS}/designSets.py ${SCRATCH_DIR} ${P5_ADAPTER_TAG} ${P7_ADAPTER_TAG} $^ ${GENOME_FASTA} $*_covering_set_primers.gff3 $*_covering_set_primers.tsv
