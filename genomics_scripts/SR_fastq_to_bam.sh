#!/bin/zsh

# script to take fastq files for the sequencing of wild D.pse SR and ST to bam files
#requires indexed ref fastas
#trimmomatic version: 0.39
#bwa version: 0.7.19
#samtools version: 1.22.1
#reference genome used was UCI_Dpse_MV25 D.pse reference genome which was submitted to RefSeq on March 3rd 2020. It can be accessed at https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009870125.1/.
#---------------------------------------------------------------------------------


filename="/path/to/fastq"
ref="/path/to/reference_fastq"

#trimmomatic
java -jar ~/genomics/trimmomatic-0.39.jar PE -threads 1 -phred33 ${filename}_r1.fastq.gz ${filename}_r2.fastq.gz ${filename}_R1_paired.fastq.gz ${filename}_R1_unpaired.fastq.gz ${filename}_R2_paired.fastq.gz ${filename}_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

#bwa
nice bwa mem -t 8 ${ref} ${filename}_R1_paired.fastq.gz ${filename}_R2_paired.fastq.gz > ${filename}_paired.sam
nice bwa mem -t 8 ${ref} ${filename}_R1_unpaired.fastq.gz ${filename}_R2_unpaired.fastq.gz > ${filename}_singles.sam

#post bwa to final bam
samtools view -@ 4 -b ${filename}_paired.sam > ${filename}_paired.bam
samtools view -@ 4 -b ${filename}_singles.sam > ${filename}_singles.bam

samtools sort -@ 4 -o ${filename}_paired_sorted.bam ${filename}_paired.bam
samtools sort -@ 4 -o ${filename}_singles_sorted.bam ${filename}_singles.bam

samtools merge -@ 4 ${filename}.bam ${filename}_singles_sorted.bam ${filename}_paired_sorted.bam

rm ${filename}_R1_paired.fastq.gz ${filename}_R1_unpaired.fastq.gz ${filename}_R2_paired.fastq.gz ${filename}_R2_unpaired.fastq.gz ${filename}_paired.sam ${filename}_singles.sam ${filename}_paired.bam ${filename}_singles.bam ${filename}_paired_sorted.bam ${filename}_singles_sorted.bam

#index BAM file
samtools index "${filename}.bam"
