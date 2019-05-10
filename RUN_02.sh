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
# FASTQ to uBAM
##########################################################
java -jar /home/alexander/Downloads/picard.jar FastqToSam F1=SQ6981_S1_L00X_MASTER_R1_001.fastq.gz O=SQ6981_S1_L00X_MASTER_R1_001.unmapped.bam SM=SQ6981_R1
java -jar /home/alexander/Downloads/picard.jar FastqToSam F1=SQ6981_S1_L00X_MASTER_R2_001.fastq.gz O=SQ6981_S1_L00X_MASTER_R2_001.unmapped.bam SM=SQ6981_R2

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
java -jar /home/alexander/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 SQ6981_S1_L00X_MASTER_R1_001.fastq.gz SQ6981_S1_L00X_MASTER_R2_001.fastq.gz SQ6981_S1_L00X_MASTER_R1_001_paired.fq.gz SQ6981_S1_L00X_MASTER_R1_001_unpaired.fq.gz SQ6981_S1_L00X_MASTER_R2_001_paired.fq.gz SQ6981_S1_L00X_MASTER_R2_001_unpaired.fq.gz ILLUMINACLIP:/home/alexander/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

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
bwa mem -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_SQ6981.sam
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
#Picard - https://github.com/broadinstitute/picard/releases/tag/2.20.0
#Version 2.20.0
#https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
java -jar /home/alexander/Downloads/picard.jar MarkDuplicates I=aligned_SQ6981.bam O=aligned_SQ6981_picard_markedduplicates.bam M=marked_dup_metrics.txt
      
java -jar picard.jar

##########################################################
# BaseRecalibration
##########################################################

java -jar GenomeAnalysisTK.jar \
   -T PrintReads \
   -R reference.fasta \
   -I input.bam \
   -BQSR recalibration_report.grp \
   -o output.bam

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
# add coverage depth with -C 30; freebayes -f ref.fa -C 5 aln.bam >var.vcf
freebayes -C 30 -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Analysis/aligned_markdup_SQ6981.bam > Freebayes_30X_aligned_markdup_SQ6981.bam.vcf

#bcftools call - https://samtools.github.io/bcftools/bcftools.html
#Version: 1.7 (using htslib 1.7-2)
#Input: bcftools mpileup
#bcftools call -v -m SQ8992_L001_R1_R2_hg38_sorted.mpileup > SQ8992_L001_R1_R2_hg38_sorted.mpileup_variants.vcf
#bcftools call -v -m samtools_aligned_markdup_SQ6981.mpileup > samtools_aligned_markdup_SQ6981.mpileup_variants.vcf
#bcftools call -v -m bcftools_aligned_markdup_SQ6981.mpileup > bcftools_aligned_markdup_SQ6981.mpileup_variants.vcf
bcftools mpileup --threads 6 -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | bcftools call --threads 6 -Ov -mv > bcftools_aligned_markdup_SQ6981.mpileup.vcf
#30X
bcftools mpileup --threads 6 -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | bcftools call --threads 6 -Ov -mv | bcftools filter -i 'QUAL>20 && DP>30' > bcftools_30X_aligned_markdup_SQ6981.mpileup.vcf

#VarScan - http://varscan.sourceforge.net/using-varscan.html
#v2.4.3 seems to not like samtools mpileup, using bcftools
#Input: bcftools mpileup

#varscan mpileup2snp samtools_aligned_markdup_SQ6981.mpileup --min-coverage 30 --output-vcf 1 > varscan_samtools_aligned_markdup_SQ6981.mpileup.vcf
#varscan mpileup2snp bcftools_aligned_markdup_SQ6981.mpileup --min-coverage 30 --output-vcf 1 > varscan_bcftools_aligned_markdup_SQ6981.mpileup.vcf
#bcftools mpileup --threads 6 -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2indel > varscan_indel_bcftools_aligned_markdup_SQ6981.mpileup.vcf 
#bcftools mpileup -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2snp > varscan_snp_bcftools_aligned_markdup_SQ6981.mpileup.vcf
samtools mpileup -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2snp --output-vcf 1 > varscan_snp_samtools_aligned_markdup_SQ6981.mpileup.vcf
#30X
samtools mpileup -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2snp --min-coverage 30 --output-vcf 1 > varscan_snp_30X_samtools_aligned_markdup_SQ6981.mpileup.vcf

#varscan mpileup2indel samtools_aligned_markdup_SQ6981.mpileup --min-coverage 30
#varscan mpileup2indel bcftools_aligned_markdup_SQ6981.mpileup --min-coverage 30
#bcftools mpileup -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2indel > varscan_indel_bcftools_aligned_markdup_SQ6981.mpileup.vcf
samtools mpileup -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2indel --output-vcf 1 > varscan_indel_samtools_aligned_markdup_SQ6981.mpileup.vcf
#30X
samtools mpileup -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam | varscan mpileup2indel --min-coverage 30 --output-vcf 1 > varscan_indel_30X_samtools_aligned_markdup_SQ6981.mpileup.vcf
#Naive Variant Caller - ignore

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
# VarScan - Combine the varscan vcf's
##########################################################
bcftools concat -o 

##########################################################
# VCF Statistics
##########################################################
bcftools stats SQ6981_S1.vcf | head -n 30

##########################################################
# Filtering
##########################################################
#Version: 1.7 (using htslib 1.7-2)
#bcftools filter -s LowQual -e '%QUAL<20 || DP>100' > var.flt.vcf
bcftools filter -i 'QUAL>20 && DP>30' GATK_output.vcf.gz > var.flt.vcf

##########################################################
# VCF Intersects
##########################################################
#https://samtools.github.io/bcftools/bcftools.html#isec
bgzip the vcfs
bcftools index the vcfs
bcftools isec -p dir -n=3 bcftools_aligned_markdup_SQ6981.mpileup.vcf.gz Freebayes_aligned_markdup_SQ6981.bam.vcf.gz GATK_output.vcf.gz
bcftools isec -p 30X -n=4 bcftools_aligned_markdup_SQ6981.mpileup.vcf.gz Freebayes_30X_aligned_markdup_SQ6981.bam.vcf.gz GATK_output.vcf.gz varscan_snpandindels_30X_samtools_aligned_markdup_SQ6981.mpileup.vcf.gz
#-n=4 means that the variants reported are in all 4 files. This is what I will pass to VEP.


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


