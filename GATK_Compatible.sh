#Whole exome to variant calling
#Following the GATK4 Best Practices
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

#PATH VARIABLE

#FILES

##########################################################
# Mapping
##########################################################

#Reference downloaded Illumina iGenomes hg19
#Reference Directory: /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta
#/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

#Index our hg19 genome fasta (take a while.)
#Version: 0.7.17-r1188
bwa index genome.fa

#bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta
#bwa mem -M -t 7 /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz /media/alexander/Elements/RQ534361-KA/Data/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > aligned_M_SQ6981.sam

cd /media/alexander/Elements/RQ534361-KA/Analysis/
bwa mem -K 100000000 -p -v 3 -t 6 -Y /media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa SQ6981_S1_L00X_MASTER_R1_001.fastq.gz SQ6981_S1_L00X_MASTER_R2_001.fastq.gz > SQ6981.aligned.unsorted.sam
