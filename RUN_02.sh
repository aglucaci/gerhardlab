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
fastqc SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
fastqc SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

##########################################################
# Mapping
##########################################################

#Mapping
#Reference downloaded Illumina iGenomes hg19
#Reference Directory: /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
#/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

#Index our hg19 genome fasta (take a while.)
bwa index genome.fa

#Reads Directory: /media/alexander/Elements/RQ534361-KA/Data
#/media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
#/media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

#BWA-Mem mapping, use -M flag to make this alignment Picard friendly (Useful for GATK)
#bwa mem -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_SQ6981.sam
bwa mem -M -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_SQ6981.sam

#SAM to BAM
samtools view -bS aligned_SQ6981.sam > aligned_SQ6981.bam

##########################################################
# MarkDuplicates & Sorting
##########################################################

#Sort our BAM file
samtools sort aligned_SQ6981.bam > aligned_sorted_SQ6981.bam

#Index our sorted BAM
samtools index aligned_sorted_SQ6981.bam

#BEDTools
bedtools genomecov -ibam aligned_sorted_SQ6981.bam -bg > aligned_sorted_SQ6981.bg

#Alignment flagstats
samtools flagstat aligned_sorted_SQ6981.bam > flagstat_aligned_sorted_SQ6981.txt

# -- MarkDuplicates - samtools
#http://www.htslib.org/doc/samtools.html
# The first sort can be omitted if the file is already name ordered
samtools sort -n -o aligned_namesorted_SQ6981.bam aligned_SQ6981.bam
# Add ms and MC tags for markdup to use later
samtools fixmate -m aligned_namesorted_SQ6981.bam aligned_fixmate_namesorted_SQ6981.bam
# Markdup needs position order
samtools sort -o aligned_positionsort_SQ6981.bam aligned_fixmate_namesorted_SQ6981.bam
# Finally mark duplicates
samtools markdup aligned_positionsort_SQ6981.bam aligned_markdup_SQ6981.bam

# -- MarkDuplicates - Picard

##########################################################
# mpileup
##########################################################

#samtools mpileup
#http://samtools.sourceforge.net/mpileup.shtml
samtools mpileup -E -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_sorted_SQ6981.bam > samtools_aligned_sorted_SQ6981.mpileup
#MarkDuplicates Version  
samtools mpileup -E -uf /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam > samtools_aligned_markdup_SQ6981.mpileup

#bcftools mpileup
#https://samtools.github.io/bcftools/bcftools.html#mpileup
bcftools mpileup -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_sorted_SQ6981.bam
bcftools mpileup -Ou -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_markdup_SQ6981.bam

##########################################################
# Variant Calling
##########################################################

#FreeBayes
#https://github.com/ekg/freebayes
#freebayes -f ref.fa aln.bam > var.vcf
freebayes -f /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa aligned_sorted_SQ6981.bam > Freebayes_aligned_sorted_SQ6981.vcf

#bcftools call
bcftools call -v -m SQ8992_L001_R1_R2_hg38_sorted.mpileup > SQ8992_L001_R1_R2_hg38_sorted.mpileup_variants.vcf
bcftools call -v -m samtools_aligned_markdup_SQ6981.mpileup > samtools_aligned_markdup_SQ6981.mpileup_variants.vcf

#varscan

#naive variant caller

#GATK

##########################################################
# Filtering
##########################################################

##########################################################
# Annotation
##########################################################

##########################################################
# VCF Intersects
##########################################################


##########################################################
# End of Script.
##########################################################


