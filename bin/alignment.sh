#! /bin/bash
#$ -cwd
#$ -r n
#$ -S /bin/bash
#$ -P gentech-rnd.p
#$ -l h_rt=4:00:00
#$ -pe pe_slots 8
#$ -w e
#$ -j y
#$ -R y

##########################################################################################
# Join reads using FLASH
##########################################################################################
module load jgibio

## inserted
module load flash
###################

reads=reads/$samp.fastq.gz
[[ ! -e "$reads" ]] && echo "$reads not found" && exit 1

## changed
## flash300 --> flash
flash -r 300 -f 324 -s 30 -o ${reads%%.*} -t 8 -I $reads


##########################################################################################
# Map reads using BWA
##########################################################################################
joined=reads/$samp.extendedFrags.fastq
[[ ! -e "$joined" ]] && echo "$joined not found" && exit 1
module load bwa/0.7.4
bwa mem -t 8 -B 2 ~/DataAnalysis/jgi_arac_analysis/ref/reference.fasta $joined > aligned/$samp.sam

module load picard
picard SortSam I=aligned/$samp.sam O=aligned/$samp.bam SO=coordinate CREATE_INDEX=true
cd aligned && ln -s $samp.bai $samp.bam.bai && cd ..

##########################################################################################
# readsam.py
##########################################################################################

### inserted
#module load biopython
########################
#
#bam=aligned/$samp.bam
#[[ ! -e "$bam" ]] && echo "$bam NOT FOUND" && exit 1
#python readsam.py $bam ref/reference.fasta > aligned/$samp.varseqs.txt
