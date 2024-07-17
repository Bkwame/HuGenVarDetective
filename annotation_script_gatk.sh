#!bin/bash


set -e

# do this to get conda working
# source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate gatk


path=/home/ngslab/Documents/adu_projects/diana/mmp9
fq1=/home/ngslab/Documents/adu_projects/diana/mmp9/reads/trimmed_paired/*_R1.paired.fastq.gz
base=$(basename ${fq1} _R1.paired.fastq.gz)
results=${path}/results/${base}
filter_snps2=${results}/${base}_filtered_snps_final.vcf
filter_indels2=${results}/${base}_filtered_indels_final.vcf
indels_vcf=/home/ngslab/Documents/adu_projects/diana/mmp9/results/indels_vcf
passed_indels=/home/ngslab/Documents/adu_projects/diana/mmp9/results/indels_vcf/passed_indels
snps_vcf=/home/ngslab/Documents/adu_projects/diana/mmp9/results/snps_vcf
passed_snps=/home/ngslab/Documents/adu_projects/diana/mmp9/results//snps_vcf/passed_snps
filter_snps3=${results}/${base}.filtered_snps_final2.vcf
filter_indels3=${results}/${base}.filtered_indels_final2.vcf

## move vcfs into their appropriate folders for annotation (snp/indel)
cd $results/..
mkdir -p $snps_vcf $indels_vcf

for folder in $results
do
    cd $folder
    gatk SelectVariants --exclude-filtered -V $filter_snps2 -O $filter_snps3
    cp ${folder}_filtered_snps_final.vcf $snps_vcf
    gatk SelectVariants --exclude-filtered -V $filter_indels2 -O $filter_indels3
    cp ${folder}_filtered_indels_final.vcf $indels_vcf
    cd ..
done

conda deactivate

## filter out poor quality variants
## for indels
cd $indels_vcf
mkdir -p $passed_indels

for file in ${indels_vcf}/*.vcf.gz
do
    base=$(basename ${file} _filtered_indels_final.vcf.gz)
    bcftools filter --exclude 'FILTER="QD_filter" || FILTER="DP_filter" || INFO/DP < 30' -o ${base}.pass.indels.vcf.gz $file
    bcftools index -c ${base}.pass.indels.vcf.gz
    mv ${base}.pass.indels.vcf.gz ${base}.pass.indels.vcf.gz.csi passed_indels/
done

##for snps
cd $snps_vcf
mkdir -p passed_snps

for file in ${snps_vcf}/*.vcf.gz
do
    base=$(basename ${file} _filtered_snps_final.vcf.gz)
    bcftools filter --exclude 'FILTER="FS_filter" || FILTER="SOR_filter" || FILTER="QD_filter" || FILTER="LowQual" || INFO/DP < 25' -o ${base}.pass.snps.vcf.gz $file 
    #bcftools filter --exclude 'FILTER="QD_filter" || FILTER="DP_filter" || INFO/DP < 25' -o ${base}.pass.variants.vcf.gz $file
    bcftools index -c ${base}.pass.snps.vcf.gz
    mv ${base}.pass.snps.vcf.gz ${base}.pass.snps.vcf.gz.csi passed_snps/
done

## Merge snps and indels 
cd $passed_snps
bcftools merge -Oz -o all.mmp9.merged.snps.vcf.gz *.pass.snps.vcf.gz

cd $passed_indels
bcftools merge -Oz -o all.mmp9.merged.indels.vcf.gz *.pass.indels.vcf.gz





















