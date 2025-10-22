#!/bin/zsh
#script that accepts a bam and outputs stats from ANGSD
#reference genome used was UCI_Dpse_MV25 D.pse reference genome which was submitted to RefSeq on March 3rd 2020. It can be accessed at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009870125.1/.
#ANGSD version used was 0.941
#SweepFinder2 was accessed at https://degiorgiogroup.fau.edu/sf2.html on Jan 12th 2025.
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
