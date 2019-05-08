##SOFTWARE USED
#BWA, SAMTOOLS is installed locally.

PICARD="/home/alexander/Downloads/picard.jar"
GATK="/home/alexander/Downloads/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar"
REFERENCE="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
GNOMAD_VCF="/media/alexander/Elements/gnomad.exomes.r2.1.1.sites.vcf.bgz"

READ1="SQ6981_S1_L00X_MASTER_R1_001.fastq.gz"
READ2="SQ6981_S1_L00X_MASTER_R2_001.fastq.gz"

WD = ""
#######################################################################
# Initialize this script.
#######################################################################
clear

echo .
echo "## Initializing script - Whole Exome to Variant Calling (c) 2019"
echo "## by Alexander G. Lucaci"
echo "## Dataset: RQ534361-KA"
echo ""

#######################################################################
# Data preprocessing
#######################################################################
echo "Changing directory to Data"
cd ../Data

echo "Combining FASTQ files.."
if [[ ! -e "../Analysis/Run_03/"$READ1 ]]
then
    cat SQ6981_S1_L001_R1_001.fastq.gz SQ6981_S1_L002_R1_001.fastq.gz SQ6981_S1_L003_R1_001.fastq.gz SQ6981_S1_L004_R1_001.fastq.gz > ../Analysis/Run_03/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
    cat SQ6981_S1_L001_R2_001.fastq.gz SQ6981_S1_L002_R2_001.fastq.gz SQ6981_S1_L003_R2_001.fastq.gz SQ6981_S1_L004_R2_001.fastq.gz > ../Analysis/Run_03/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz
fi

#if [ ! -f "/etc/bebebe" ]; then echo "The file does not exist"; fi

echo "Changing directory to Data to ../Analysis/Run_03"
cd ../Analysis/Run_03

##########################################################
# FASTQ to uBAM
##########################################################
#echo "Generating unmapped BAMs from FASTQ"
#java -jar $PICARD FastqToSam F1=$READ1 O=SQ6981_S1_L00X_MASTER_R1_001.unmapped.bam SM=SQ6981_R1
#java -jar $PICARD FastqToSam F1=$READ2 O=SQ6981_S1_L00X_MASTER_R2_001.unmapped.bam SM=SQ6981_R2

##########################################################
# Run FASTQC
##########################################################
#Version 0.11.5
#echo "Generating FASTQC files"
#pwd
#fastqc $READ1
#fastqc $READ2
#Giving me a Java 9/8 error. See how to resolve this.

##########################################################
# Mapping 
##########################################################
#echo "Indexing reference genome"
#bwa index $REFERENCE

echo "Mapping reads to reference"
if [[ ! -e "aligned_SQ6981.sam" ]]
then
    bwa mem -M -t 7 -R "@RG\tID:SQ6981\tSM:SQ6981\tPL:Illumina" $REFERENCE $READ1 $READ2 > aligned_SQ6981.sam 
fi

##########################################################
# Prepping for analysis 
##########################################################
echo "Converting SAM to BAM"
if [[ ! -e "aligned_SQ6981.bam" ]]
then
    samtools view -bS aligned_SQ6981.sam > aligned_SQ6981.bam
fi

echo "Sorting BAM"
if [[ ! -e "sorted_aligned_SQ6981.bam" ]]
then
    java -jar $PICARD SortSam I=aligned_SQ6981.bam O=sorted_aligned_SQ6981.bam SORT_ORDER=coordinate
fi

echo "Marking PCR Duplicates"
if [[ ! -e "marked_duplicates_sorted_aligned_SQ6981.bam" ]]
then
    java -jar $PICARD MarkDuplicates I=sorted_aligned_SQ6981.bam O=marked_duplicates_sorted_aligned_SQ6981.bam M=marked_dup_metrics.txt
fi

echo "Indexing the gnomAD VCF"
java -jar $GATK IndexFeatureFile -F $GNOMAD_VCF

#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
echo "Generated recalibration table based on gnomAD"
java -jar $GATK BaseRecalibrator -I marked_duplicates_sorted_aligned_SQ6981.bam -R $REFERENCE --known-sites $GNOMAD_VCF -O recal_data.table

echo "Applying base recalibration"
java -jar $GATK ApplyBQSR -R $REFERENCE -I marked_duplicates_sorted_aligned_SQ6981.bam --bqsr-recal-file recal_data.table -O BQSR_marked_duplicates_sorted_aligned_SQ6981.bam

##########################################################
# Variant Calling
##########################################################
QUERY="BQSR_marked_duplicates_sorted_aligned_SQ6981.bam"

freebayes -C 30 -f $REFERENCE $QUERY > Freebayes_30X_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.vcf

bcftools mpileup --threads 6 -Ou -f $REFERENCE $QUERY | bcftools call --threads 6 -Ov -mv | bcftools filter -i 'QUAL>20 && DP>30' > bcftools_30X_aligned_markdup_SQ6981.mpileup.vcf

samtools mpileup -f $REFERENCE $QUERY | varscan mpileup2snp --min-coverage 30 --output-vcf 1 > varscan_snp_30X_samtools_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.mpileup.vcf

samtools mpileup -f $REFERENCE $QUERY | varscan mpileup2indel --min-coverage 30 --output-vcf 1 > varscan_indel_30X_samtools_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.mpileup.vcf

java -jar $GATK --java-options "-Xmx4g" HaplotypeCaller -R $REFERENCE-I $QUERY -O GATK_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.vcf.gz -bamout GATK_bamout.bam

##########################################################
# Consolidate VCFs
##########################################################




##########################################################
# Intersect VCFs
##########################################################



##########################################################
# Annotate VCFs - VEP
##########################################################


##########################################################
# Filter VCF
##########################################################




