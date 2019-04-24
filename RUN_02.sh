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

# Start script 
clear

# Add title and other info to display here.
echo . 

#COMBINE RUNS INTO LARGER FASTQ
cat SQ6981_S1_L001_R1_001.fastq.gz SQ6981_S1_L002_R1_001.fastq.gz SQ6981_S1_L003_R1_001.fastq.gz SQ6981_S1_L004_R1_001.fastq.gz > SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
cat SQ6981_S1_L001_R2_001.fastq.gz SQ6981_S1_L002_R2_001.fastq.gz SQ6981_S1_L003_R2_001.fastq.gz SQ6981_S1_L004_R2_001.fastq.gz > SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

#RUN FASTQC
fastqc SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
fastqc SQ6981_S1_L00X_MASTER_R2_001.fastq.gz

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
bwa mem -M -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_SQ6981.sam

#SAM to BAM
samtools view -bS aligned_SQ6981.sam > aligned_SQ6981.bam

#Sort our BAM file
samtools sort aligned_SQ6981.bam > aligned_sorted_SQ6981.bam

#Index our BAM
samtools index aligned_sorted_SQ6981.bam

#Alignment flagstats
samtools flagstat aligned_sorted_SQ6981.bam


