#!/bin/sh

# usage: sweep_script.sh manifest.tsv reference.fasta

VAR_TSV=$1
VAR_REF=$2

mkdir QMS/
mkdir QMS/TRIMMING/
mkdir QMS/ALIGNING/
mkdir QMS/EXPECTED/
mkdir QMS/RECALIBRATING/
mkdir QMS/LOFREQ/

# copy reference genome to current working directory
cp "$VAR_REF" QMS/ALIGNING/reference.fasta

# build bowtie2 database for reference 
bowtie2-build QMS/ALIGNING/reference.fasta QMS/ALIGNING/reference

# read QC
while read i j k 
do

    # copy and rename raw FASTQ files
    cp "$j" QMS/TRIMMING/"$i"_R1.fastq.gz
    cp "$k" QMS/TRIMMING/"$i"_R2.fastq.gz

    # trim and filter
    trim_galore \
        --paired \
        QMS/TRIMMING/"$i"_R1.fastq.gz \
        QMS/TRIMMING/"$i"_R2.fastq.gz \
        -o QMS/TRIMMING/ \
        -q 30 \
        --phred33 \
        --fastqc \
        --length 50

    # remove copied and renamed raw FASTQ files
    rm -f QMS/TRIMMING/"$i"_R1.fastq.gz
    rm -f QMS/TRIMMING/"$i"_R2.fastq.gz   

done < "$VAR_TSV"

while read i j k
do

    # align reads to reference
    bowtie2 \
        -p 4 \
        --no-mixed \
        -x QMS/ALIGNING/reference \
        -1 QMS/TRIMMING/"$i"_R1_val_1.fq.gz \
        -2 QMS/TRIMMING/"$i"_R2_val_2.fq.gz | \
            samtools sort -@ 4 -o QMS/ALIGNING/"$i"_sorted.bam

done < "$VAR_TSV"

while read i j k
do

    # identify SNVs with >1% prevalence to give to GATK later as 'expected SNVs' 
    bcftools mpileup \
        -Ou \
        -d 1000 \
        -f QMS/ALIGNING/reference.fasta QMS/ALIGNING/"$i"_sorted.bam | \
            bcftools call -mv -Ov | \
                    bcftools view -q 0.01 -o QMS/EXPECTED/"$i"_bcftools.vcf

done < "$VAR_TSV"

# base quality recalibration with GATK
while read i j k
do

    # add reads to a read group
    picard AddOrReplaceReadGroups \
        I=QMS/ALIGNING/"$i"_sorted.bam \
        O=QMS/RECALIBRATING/"$i"_grouped.bam \
        RGID=1 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM="$i"

    # index feature file
    gatk IndexFeatureFile -I QMS/EXPECTED/"$i"_bcftools.vcf

    # create FASTA dictionary file
    gatk CreateSequenceDictionary -R QMS/ALIGNING/reference.fasta

    # base recalibration
    gatk BaseRecalibrator \
        -I QMS/RECALIBRATING/"$i"_grouped.bam \
        -R QMS/ALIGNING/reference.fasta \
        --known-sites QMS/EXPECTED/"$i"_bcftools.vcf \
        -O QMS/RECALIBRATING/"$i"_recal_data.table

    # apply BQSR
    gatk ApplyBQSR \
        -I QMS/RECALIBRATING/"$i"_grouped.bam \
        -R QMS/ALIGNING/reference.fasta \
        --bqsr-recal-file QMS/RECALIBRATING/"$i"_recal_data.table \
        -O QMS/RECALIBRATING/"$i"_recal.bam

    # sort recalibrated BAM files
    samtools sort \
        -@ 8 \
        -o QMS/RECALIBRATING/"$i"_recal_sorted.bam \
        QMS/RECALIBRATING/"$i"_recal.bam

done < "$VAR_TSV"

while read i j k 
do

    # lofreq indelqual
    lofreq indelqual \
        --dindel \
        --ref QMS/ALIGNING/reference.fasta \
        -o QMS/LOFREQ/"$i"_dindel.bam \
        QMS/RECALIBRATING/"$i"_recal_sorted.bam

    # lofreq call variants 
    lofreq call \
        -f QMS/ALIGNING/reference.fasta \
        -o QMS/LOFREQ/"$i"_variants.vcf \
        --sig 0.05 \
        --call-indels QMS/LOFREQ/"$i"_dindel.bam

done < "$VAR_TSV"

while read i j k 
do

    rm -f QMS/ALIGNING/"$i"_sorted.bam \
        QMS/LOFREQ/"$i"_dindel.bam \
        QMS/RECALIBRATING/"$i"*

done < "$VAT_TSV"
