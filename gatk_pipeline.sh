#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# This script is for demonstration purposes only

# if false
# then
# download data
# wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
# wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


# echo "Run Prep files..."

################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################



# download reference files
# wget -P ~/Desktop/demo/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# gunzip ~/Desktop/demo/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
# samtools faidx ~/Desktop/demo/supporting_files/hg38/hg38.fa


# ref dict - .dict file before running haplotype caller
# gatk CreateSequenceDictionary R=~/Desktop/demo/supporting_files/hg38/hg38.fa O=~/Desktop/demo/supporting_files/hg38/hg38.dict


# download known sites files for BQSR from GATK resource bundle
# wget -P ~/Desktop/demo/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
# wget -P ~/Desktop/demo/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx



###################################################### VARIANT CALLING STEPS ####################################################################
# fi

set -e


eval "$(conda shell.bash hook)"
conda activate gatk

## kindly specify the pathname for the path variable
path=<path/to/experiment_name>
cd $path

# mkdir -p supporting_files data results aligned_reads

# directories/variables (kindly specify the pathname for the ref variable)
ref=${path}/supporting_files/reference.fasta
# known_sites="${path}/supporting_reads/DF0114.merged.vcf"
aligned_reads=${path}/aligned_reads/${base}
reads=${path}/reads/trimmed_paired
data=${path}/data/${base}
results=${path}/results/${base}
sam=${aligned_reads}/${base}.sam
sort_bam=${aligned_reads}/${base}_sorted_dedup_reads.bam
# sort_bqsr_bam="${aligned_reads}/${base}_sorted_dedup_bqsr_reads.bam"
recal_bam=${aligned_reads}/${base}_recal_reads.bam
post_recal_table=${data}/${base}_post_recal_data.table
recal_plot=${data}/${base}_recalibration_plots.pdf


# -------------------------------------------
echo "DATA PRE-PROCESSING FOR VARIANT CALLING"
# -------------------------------------------

## Reference prep
#samtools faidx <reference.fasta>
#gatk CreateSequenceDictionary -R <reference.fasta> -O ${path}/supporting_files/ref.dict

# -------------------
# STEP 1: QC - Run fastqc 
# -------------------
# echo "STEP 1: QC - Run fastqc"
# fastqc ${reads}/*.fastq.gz -o ${reads}/raw_quality
# Trimming required afterwards followed by 

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------
echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index $ref

# BWA alignment (indicate the name portion of all the trimmed reads that are consistent across all samples eg. _trimmed.R1.fastq.gz)
for fq1 in ${reads}/*_trimmed_forward_reads.fastq.gz
do
    base=$(basename ${fq1} _trimmed_forward_reads.fastq.gz)
    
    # known_sites="${path}/supporting_reads/DF0114.merged.vcf"
    aligned_reads=${path}/aligned_reads/${base}
    reads=${path}/reads/trimmed_paired
    fq2=${reads}/${base}_trimmed_reverse_reads.fastq.gz
    results=${path}/results/${base}
    data=${path}/data/${base}
    sam=${aligned_reads}/${base}.sam
    sort_bam=${aligned_reads}/${base}_sorted_dedup_reads.bam
    # sort_bqsr_bam="${aligned_reads}/${base}_sorted_dedup_bqsr_reads.bam"
    raw_vcf=${results}/${base}_raw_variants.vcf
    raw_snps=${results}/${base}_raw_snps.vcf
    raw_indels=${results}/${base}_raw_indels.vcf
    filter_snps1=${results}/${base}_filtered_snps.vcf
    filter_indels1=${results}/${base}_filtered_indels.vcf
    bqsr_snps=${results}/${base}_bqsr_snps.vcf
    bqsr_indels=${results}/${base}_bqsr_indels.vcf
    aln_metrics=${aligned_reads}/${base}_alignment_metrics.txt
    ins_metrics=${aligned_reads}/${base}_insert_size_metrics.txt
    ins_hist=${aligned_reads}/${base}_insert_size_histogram.pdf
    recal_table=${data}/${base}_recal_data.table
    recal_bam=${aligned_reads}/${base}_recal_reads.bam
    post_recal_table=${data}/${base}_post_recal_data.table
    recal_plot=${data}/${base}_recalibration_plots.pdf
    raw_recal_vcf=${results}/${base}_raw_recal_variants.vcf
    raw_recal_snps=${results}/${base}_raw_recal_snps.vcf
    raw_recal_indels=${results}/${base}_raw_recal_indels.vcf
    filter_snps2=${results}/${base}_filtered_snps_final.vcf
    filter_indels2=${results}/${base}_filtered_indels_final.vcf
    filter_snps3=${results}/${base}.filtered_snps_final2.vcf
    filter_indels3=${results}/${base}.filtered_indels_final2.vcf
    
    ## Create sample-specific directories for each result directory
    mkdir -p $aligned_reads $results $data
    
    ## Edit the Read Group name according to the specification of the sequencing device used
    bwa mem -t 12 -R "@RG\tID:${base}\tPL:Illumina\tPM:NextSeq2k\tSM:${base}" $ref $fq1 $fq2 > $sam

# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------
    echo "STEP 3: Mark Duplicates and Sort - GATK4"
    gatk MarkDuplicatesSpark -I $sam -O $sort_bam
    rm -r $sam

# -----------------------------------------------
# STEP 4: Collect Alignment & Insert Size Metrics
# -----------------------------------------------
    echo "STEP 4: Collect Alignment & Insert Size Metrics"
    gatk CollectAlignmentSummaryMetrics R=$ref I=$sort_bam O=$aln_metrics
    gatk CollectInsertSizeMetrics INPUT=$sort_bam OUTPUT=$ins_metrics HISTOGRAM_FILE=$ins_hist
   ### multiqc can be used to generate summaries on the metrics above (as long as they are .txt/.html files)

# ---------------------------------------------
# STEP 5: Call Variants - gatk HaplotypeCaller
# ---------------------------------------------
    echo "STEP 5: Call Variants - gatk HaplotypeCaller"
    gatk HaplotypeCaller -R $ref -I $sort_bam -O $raw_vcf

# -----------------------------
# STEP 6: Extract SNPs & INDELS
# -----------------------------
    echo "STEP 6: Extract SNPs & INDELS"
    gatk SelectVariants -R $ref -V $raw_vcf --select-type SNP -O $raw_snps
    gatk SelectVariants -R $ref -V $raw_vcf --select-type INDEL -O $raw_indels

# ----------------------------------------
# STEP 7: Filter SNPs and INDELs
# ----------------------------------------
    echo "STEP 7: Filter SNPs and INDELs"
    gatk VariantFiltration -R $ref -V $raw_snps -O $filter_snps1 -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    gatk VariantFiltration -R $ref -V $raw_indels -O $filter_indels1 -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

# ---------------------------------
# STEP 8: Exclude Filtered Variants
# ---------------------------------
    echo "STEP 8: Exclude Poor Quality Filtered Variants"
    gatk SelectVariants --exclude-filtered -V $filter_snps1 -O $bqsr_snps
    gatk SelectVariants --exclude-filtered -V $filter_indels1 -O $bqsr_indels

# ------------------------------------------
# STEP 9: Base Quality Score Recalibration 1
# ------------------------------------------
    echo "STEP 9: Base Quality Score Recalibration 1"
# 1. build the model
    gatk BaseRecalibrator -I $sort_bam -R $ref --known-sites $bqsr_snps --known-sites $bqsr_indels -O $recal_table --maximum-cycle-value 500000
# 2. Apply the model to adjust the base quality scores
    gatk ApplyBQSR -I $sort_bam -R $ref --bqsr $recal_table -O $recal_bam

# ------------------------------------------
# STEP 10: Base Quality Score Recalibration 2
# ------------------------------------------
    echo "STEP 10: Base Quality Score Recalibration 2"
# 1. build the model
    gatk BaseRecalibrator -I $recal_bam -R $ref --known-sites $bqsr_snps --known-sites $bqsr_indels -O $post_recal_table --maximum-cycle-value 500000

# ---------------------------
# STEP 11: Analyze Covariates
# ---------------------------
    echo "STEP 11: Analyze Covariates"
    gatk AnalyzeCovariates -before $recal_table -after $post_recal_table -plots $recal_plot

# ---------------------------------------------
# STEP 12: Call Variants - gatk HaplotypeCaller
# ---------------------------------------------
    echo "STEP 12: Call Variants - gatk HaplotypeCaller"
    gatk HaplotypeCaller -R $ref -I $recal_bam -O $raw_recal_vcf

# -----------------------------
# STEP 13: Extract SNPs & INDELS
# -----------------------------
    echo "STEP 13: Extract SNPs & INDELS"
    gatk SelectVariants -R $ref -V $raw_recal_vcf --select-type SNP -O $raw_recal_snps
    gatk SelectVariants -R $ref -V $raw_recal_vcf --select-type INDEL -O $raw_recal_indels

# -------------------------------
# STEP 14: Filter SNPs and INDELs
# -------------------------------
    echo "STEP 14: Filter SNPs and INDELs"
    gatk VariantFiltration -R $ref -V $raw_recal_snps -O $filter_snps2 -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    gatk VariantFiltration -R $ref -V $raw_recal_indels -O $filter_indels2 -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

    gatk SelectVariants --exclude-filtered -V $filter_snps2 -O $filter_snps3
    gatk SelectVariants --exclude-filtered -V $filter_indels2 -O $filter_indels3

done

conda deactivate


