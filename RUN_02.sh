## REPOSITORY (RQ534361)
## PROJECT NAME - KA
## INVITAE PAIRED END RUNS 
## FASTQ FILES
#SQ6981_S1_L001_R1_001.fastq.gz
#SQ6981_S1_L002_R1_001.fastq.gz 
#SQ6981_S1_L003_R1_001.fastq.gz 
#SQ6981_S1_L004_R1_001.fastq.gz
#SQ6981_S1_L001_R2_001.fastq.gz
#SQ6981_S1_L002_R2_001.fastq.gz 
#SQ6981_S1_L003_R2_001.fastq.gz 
#SQ6981_S1_L004_R2_001.fastq.gz

#Computer: this script is configured to run on the ubuntu-deploy server.
#java -version
#openjdk version "1.8.0_191"
#OpenJDK Runtime Environment (build 1.8.0_191-8u191-b12-2ubuntu0.18.04.1-b12)
#OpenJDK 64-Bit Server VM (build 25.191-b12, mixed mode)

#Strategies
#https://software.broadinstitute.org/gatk/best-practices/

#working directory
#/media/alexander/Elements/RQ534361-KA/Data

# Start script 
clear

##########################################################
# Add title and other info to display here.
##########################################################

echo . 

##########################################################
# Data integrity check - sha256, done in another file.
##########################################################

# See other script.

##########################################################
# Combine fastq files 
##########################################################

#COMBINE RUNS INTO LARGER FASTQ
cat SQ6981_S1_L001_R1_001.fastq.gz SQ6981_S1_L002_R1_001.fastq.gz SQ6981_S1_L003_R1_001.fastq.gz SQ6981_S1_L004_R1_001.fastq.gz > SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
cat SQ6981_S1_L001_R2_001.fastq.gz SQ6981_S1_L002_R2_001.fastq.gz SQ6981_S1_L003_R2_001.fastq.gz SQ6981_S1_L004_R2_001.fastq.gz > SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

##########################################################
# FASTQ Quality Checking
##########################################################

#RUN FASTQC
#Version 0.11.5
fastqc SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
fastqc SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

##########################################################
# Read Trimming - Trimmomatic or other?
##########################################################
#Quality is very good in this dataset.

##########################################################
# Mapping
##########################################################

#Mapping
#Reference downloaded Illumina iGenomes hg19
#Reference Directory: /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
#/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

#Index our hg19 genome fasta (take a while.)
#Version: 0.7.17-r1188
bwa index genome.fa

#Reads Directory: /media/alexander/Elements/RQ534361-KA/Data
#/media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
#/media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

#BWA-Mem mapping, use -M flag to make this alignment Picard friendly (Useful for GATK)
#Version: 0.7.17-r1188
#bwa mem -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_SQ6981.sam
bwa mem -M -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_M_SQ6981.sam

#SAM to BAM
#Version: 1.7 (using htslib 1.7-2)
samtools view -bS aligned_SQ6981.sam > aligned_SQ6981.bam

##########################################################
# MarkDuplicates & Sorting
##########################################################

#Sort our BAM file
#Version: 1.7 (using htslib 1.7-2)
#http://www.htslib.org/doc/samtools-1.1.html
samtools sort aligned_SQ6981.bam > aligned_sorted_SQ6981.bam
#NameSort the BAM file
samtools sort -n -o aligned_namesorted_SQ6981.bam aligned_SQ6981.bam

#Index our sorted BAM
#Version: 1.7 (using htslib 1.7-2)
#samtools index aligned_sorted_SQ6981.bam
samtools index aligned_namesorted_SQ6981.bam

#BEDTools - Generate a bedgraph of coverage
#Version:   v2.26.0
bedtools genomecov -ibam aligned_sorted_SQ6981.bam -bg > aligned_sorted_SQ6981.bedgraph
bedtools genomecov -ibam aligned_namesorted_SQ6981.bam -bg > aligned_namesorted_SQ6981.bedgraph

#Alignment flagstats
#Version: 1.7 (using htslib 1.7-2)
samtools flagstat aligned_sorted_SQ6981.bam > flagstat_aligned_sorted_SQ6981.txt
samtools flagstat aligned_namesorted_SQ6981.bam > flagstat_aligned_namesorted_SQ6981.txt

# -- MarkDuplicates - samtools
#Version: 1.7 (using htslib 1.7-2)
#http://www.htslib.org/doc/samtools.html
# The first sort can be omitted if the file is already name ordered
#samtools sort -n -o aligned_namesorted_SQ6981.bam aligned_SQ6981.bam
# Add ms and MC tags for markdup to use later
samtools fixmate -m aligned_namesorted_SQ6981.bam aligned_fixmate_namesorted_SQ6981.bam
# Markdup needs position order
samtools sort -o aligned_positionsort_SQ6981.bam aligned_fixmate_namesorted_SQ6981.bam
# Finally mark duplicates
samtools markdup aligned_positionsort_SQ6981.bam aligned_markdup_SQ6981.bam

# -- MarkDuplicates - Picard

##########################################################
# mpileup - makes genotype likelihood calls, does not do variant calling
##########################################################

#samtools mpileup
#Version: 1.7 (using htslib 1.7-2)
#http://samtools.sourceforge.net/mpileup.shtml
samtools mpileup -E -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_sorted_SQ6981.bam > samtools_aligned_sorted_SQ6981.mpileup
#MarkDuplicates Version  
samtools mpileup -E -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam > samtools_aligned_markdup_SQ6981.mpileup

#bcftools mpileup
#Version: 1.7 (using htslib 1.7-2)
#https://samtools.github.io/bcftools/bcftools.html#mpileup
#https://www.biostars.org/p/335121/
bcftools mpileup -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_sorted_SQ6981.bam > bcftools_aligned_sorted_SQ6981.mpileup 
bcftools mpileup --threads 6 -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam > bcftools_aligned_markdup_SQ6981.mpileup 
bcftools mpileup --threads 6 -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | bcftools call --threads 6 -Ov -mv > bcftools_aligned_markdup_SQ6981.mpileup.vcf
##########################################################
# Variant Calling
##########################################################

#FreeBayes 
#version:  v1.0.0
#Input: as BAM
#https://github.com/ekg/freebayes
#freebayes -f ref.fa aln.bam > var.vcf
freebayes -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_sorted_SQ6981.bam > Freebayes_aligned_sorted_SQ6981.vcf
freebayes -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Analysis/aligned_markdup_SQ6981.bam > Freebayes_aligned_markdup_SQ6981.bam.vcf
#freebayes -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -@ /media/alexander/Elements/RQ534361-KA/Analysis/samtools_aligned_markdup_SQ6981.mpileup > Freebayes_samtools_aligned_markdup_SQ6981.mpileup.vcf

#bcftools call - https://samtools.github.io/bcftools/bcftools.html
#Version: 1.7 (using htslib 1.7-2)
#Input: bcftools mpileup
bcftools call -v -m SQ8992_L001_R1_R2_hg38_sorted.mpileup > SQ8992_L001_R1_R2_hg38_sorted.mpileup_variants.vcf
bcftools call -v -m samtools_aligned_markdup_SQ6981.mpileup > samtools_aligned_markdup_SQ6981.mpileup_variants.vcf
bcftools call -v -m bcftools_aligned_markdup_SQ6981.mpileup > bcftools_aligned_markdup_SQ6981.mpileup_variants.vcf

#VarScan - http://varscan.sourceforge.net/using-varscan.html
#v2.4.3 seems to not like samtools mpileup, using bcftools
#Input: bcftools mpileup

#varscan mpileup2snp samtools_aligned_markdup_SQ6981.mpileup --min-coverage 30 --output-vcf 1 > varscan_samtools_aligned_markdup_SQ6981.mpileup.vcf
#varscan mpileup2snp bcftools_aligned_markdup_SQ6981.mpileup --min-coverage 30 --output-vcf 1 > varscan_bcftools_aligned_markdup_SQ6981.mpileup.vcf
#bcftools mpileup --threads 6 -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2indel > varscan_indel_bcftools_aligned_markdup_SQ6981.mpileup.vcf 
bcftools mpileup -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2snp > varscan_snp_bcftools_aligned_markdup_SQ6981.mpileup.vcf

#varscan mpileup2indel samtools_aligned_markdup_SQ6981.mpileup --min-coverage 30
#varscan mpileup2indel bcftools_aligned_markdup_SQ6981.mpileup --min-coverage 30
bcftools mpileup -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2indel > varscan_indel_bcftools_aligned_markdup_SQ6981.mpileup.vcf

#VarScan - Combine the varscan vcf's

#Naive Variant Caller

## --- GATK --- ##
#USING GATK 4.1.1.0 - https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/
#On ubuntu deploy - /home/alexander/Downloads/gatk-4.1.1.0
#Input: as BAM

#http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile
#Picard - https://github.com/broadinstitute/picard/releases/tag/2.20.0
#Version 2.20.0
java -jar picard.jar ValidateSamFile I=/media/alexander/Elements/RQ534361-KA/Analysis/aligned_markdup_SQ6981.bam MODE=SUMMARY
#https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_AddOrReplaceReadGroups.php
java -jar picard.jar AddOrReplaceReadGroups I=/media/alexander/Elements/RQ534361-KA/Analysis/aligned_markdup_SQ6981.bam O=/media/alexander/Elements/RQ534361-KA/Analysis/AddOrReplaceReadGroups_aligned_markdup_SQ6981.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

#GATK seems to only run on Java8 - https://www.digitalocean.com/community/tutorials/how-to-install-java-with-apt-on-ubuntu-18-04
./gatk --java-options "-Xmx4g" HaplotypeCaller -R /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I /media/alexander/Elements/RQ534361-KA/Analysis/aligned_markdup_SQ6981.bam -O GATK_output.vcf -bamout GATK_bamout.bam
./gatk --java-options "-Xmx4g" HaplotypeCaller -R /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I /media/alexander/Elements/RQ534361-KA/Analysis/AddOrReplaceReadGroups_aligned_markdup_SQ6981.bam -O GATK_output.vcf -bamout GATK_bamout.bam
./gatk --java-options "-Xmx4g" HaplotypeCaller -R /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I /media/alexander/Elements/RQ534361-KA/Analysis/AddOrReplaceReadGroups_aligned_markdup_SQ6981.bam -O GATK_output.vcf.gz -bamout GATK_bamout.bam

##########################################################
# Filtering
##########################################################
#Version: 1.7 (using htslib 1.7-2)
bcftools filter -s LowQual -e '%QUAL<20 || DP>100' > var.flt.vcf

##########################################################
# VCF Intersects
##########################################################

##########################################################
# Annotation
##########################################################

##########################################################
# Annotation - Filtering
##########################################################

##########################################################
# Plots
##########################################################

##########################################################
# End of Script.
##########################################################


