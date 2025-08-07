#!/bin/zsh
#script that accepts a bam and outputs stats from ANGSD
#-------------------------------------------------------
directory="/path/to/bams"
run_name=""
bam_list=""
ref="/path/to/ref_fasta"
win_size="50000"
win_step="5000"

cd ${directory}
angsd -P 8 -bam ${bam_list} -doSaf 1 -out ${run_name} -anc ${ref} -GL 2
realSFS ${run_name}.saf.idx -P 24 -fold 1 > ${run_name}.sfs
realSFS saf2theta ${run_name}.saf.idx -outname ${run_name} -sfs ${run_name}.sfs -fold 1
thetaStat do_stat ${run_name}.thetas.idx -win ${win_size} -step ${win_step} -outnames ${run_name}.thetasWindow.gz
