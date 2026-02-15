#!/bin/bash

# make_beagle_refs.sh
# This script processes 1000 Genomes Phase 3 VCF files to create reference panels for Beagle imputation.

# Download 1000 Genomes Phase 3 VCF files 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/*

# Load bcftools
module load bcftools

# Set path to hg38 reference genome
ref_file="/projects/standard/aventeic/balay011/references/reference_genome/GRCh38.primary_assembly.genome.fa"

# To add 'chr' prefix, make file with old and new chromosome names
rm -f chr_names.txt
for CHR in {1..22} X; do 
    echo "${CHR} chr${CHR}" >> chr_names.txt
done

# Add 'chr' prefix -> 
# remove variants with < 3 minor allele count -> 
# subset SNPs and INDELs -> 
# left-align indels -> 
# set ID to CHR_POS_REF_ALT -> 
# remove duplicates -> 
# remove variants with missing genotypes -> 
# compress and index
for CHR in {1..22} X; do
    echo "Processing Chr ${CHR}..."
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -Ou | \
    bcftools view -e 'INFO/AC<3' -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f $ref_file -Ou | \
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou | \
    bcftools norm -d any -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz
done

for CHR in X; do
    echo "Processing Chr ${CHR}..."
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -Ou | \
    bcftools view -e 'INFO/AC<3' -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f $ref_file -Ou | \
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou | \
    bcftools norm -d any -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz
done

# Convert haploid genotypes to homozygous diploid on chrX (i.e. 0 or 0/0 to 0|0)
echo "chrX 1 156040895 M 2" > ploidy.txt
bcftools +fixploidy \
    1000GP_chrX.vcf.gz -Ov -- -p ploidy.txt | \
    sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Oz -o "1000GP_chr23.vcf.gz"
mv 1000GP_chr23.vcf.gz 1000GP_chrX.vcf.gz

# Delete original VCF files and intermediate files
rm -f ALL.chr*.vcf.gz* chr_names.txt ploidy.txt

# Convert vcf.gz to bref
wget https://faculty.washington.edu/browning/beagle/bref.27Jan18.7e1.jar
for CHR in {1..22} X; do
    echo "Creating bref for Chr ${CHR}..."
    java -jar bref.27Jan18.7e1.jar \
        1000GP_chr${CHR}.vcf.gz > 1000GP_chr${CHR}.bref
done

echo "Done."
