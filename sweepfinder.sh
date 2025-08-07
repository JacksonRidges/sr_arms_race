#!/bin/sh
set -e

directory="/path/to/bams/"
run_name=""
threads="8"
cd ${directory}

#use angsd to generate .maf file to estimate allele frequency. Used D.pse MV2-25 reference genome as the ancestral genome.
angsd \
    -GL 2 \
    -out ${run_name} \
    -nThreads ${threads} \
    -doMajorMinor 5 \
    -doMaf 1 \
    -SNP_pval 1e-6 \
    -nInd 14 \
    -bam SRbams.filelist \
    -anc dpse_X.fasta

#python script to convert .maf file from angsd into SF2 format
gzip -d -c ${run_name}.mafs.gz | \
maf_to_sf.py > ${run_name}.aff

sweepfinder2 -s 10000 ${run_name}.aff ${run_name}.SF

rm ${run_name}.arg
