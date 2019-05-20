## REPOSITORY (RQ553841)
## Start date 5/18/2019
## PROJECT NAME - Ca2+

#LINK TO REFERENCE - Illumina iGenomes
#ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz

#LINK TO GNOMAD VCF - https://gnomad.broadinstitute.org/downloads
#https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

## Software Versions
PICARD="/home/alexander/Downloads/picard.jar"
GATK="/home/alexander/Downloads/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar"


#REFERENCE="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#REFERENCE="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#REFERENCE_INDEX="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
#REFERENCE_DICT="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"

## Reference files ##
REFERENCE="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
REFERENCE_BWA="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
REFERENCE_BWA_INDEX="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.bwt"
REFERENCE_SAMTOOLS_INDEX="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
REFERENCE_DICT="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"

## gnomAD Annotation ##
GNOMAD_VCF="/media/alexander/Elements/gnomad.exomes.r2.1.1.sites.vcf.bgz"
GNOMAD_VCF_INDEX="/media/alexander/Elements/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi"

#CONCATENATED FORWARD AND REVERSE READS
READ1="SQ6981_S1_L00X_MASTER_R1_001.fastq.gz"
READ2="SQ6981_S1_L00X_MASTER_R2_001.fastq.gz"

#Set WD
WD="/media/alexander/Elements/RQ553841_Ca2+/Scripts"

#Set Project name
PROJECT_NAME="RQ553841-Ca2+"
PROJECT_SQ="SQ7711"

#######################################################################
# Initialize this script.
#######################################################################
clear

cd $WD
echo $WD
echo .
echo "## Initializing script - Whole Exome to Variant Calling (c) 2019"
echo "## by Alexander G. Lucaci"
echo "## Dataset: "$PROJECT_NAME
echo ""

#######################################################################
# Data preprocessing
#######################################################################
echo "# STARTING WES -> VCF PIPELINE"
echo "() Combining FASTQ files - "$PROJECT_NAME


echo "() Combining forward reads"
if [[ ! -e "../Analysis/"$READ1 ]]
then
    cat ../Data/SQ7711_S1_L001_R1_001.fastq.gz ../Data/SQ7711_S1_L002_R1_001.fastq.gz > ../Analysis/$READ1
fi

echo "() Combining reverse reads"
if [[ ! -e "../Analysis/"$READ2 ]]
then
    cat ../Data/SQ7711_S1_L001_R2_001.fastq.gz ../Data/SQ7711_S1_L002_R2_001.fastq.gz > ../Analysis/$READ2
fi

##########################################################
# Mapping 
##########################################################
#Had to replace chr1 to 1 so that the gnomad VCF agrees with the reference, gnomAD uses hg19.
#cat $REFERENCE | perl -pe 's/chr//g' > $REFERENCE

#bwa index $REFERENCE
#head $REFERENCE

cd ../Analysis

echo "() Mapping reads to reference"
if [[ ! -e "aligned_"$PROJECT_SQ".sam" ]]
then
    bwa mem -M -t 7 -R "@RG\tID:"$PROJECT_SQ"\tSM:"$PROJECT_SQ"\tPL:Illumina" $REFERENCE ../Analysis/$READ1 ../Analysis/$READ2 > ../Analysis/"aligned_"$PROJECT_SQ".sam" 
fi

echo "() Converting SAM to BAM"
if [[ ! -e "aligned_"$PROJECT_SQ".bam" ]]
then
    samtools view -bS "aligned_"$PROJECT_SQ".sam" > "aligned_"$PROJECT_SQ".bam"
fi

##########################################################
# Analyze the Mapping
##########################################################
echo "() Generating alignment flagstats"
if [[ ! -e "flagstat_aligned_"$PROJECT_SQ".bam.flagstats" ]]
then
    samtools flagstat "aligned_"$PROJECT_SQ".bam" > "flagstat_aligned_"$PROJECT_SQ".bam.flagstats"
fi

echo "() BEDTools - Generate a bedgraph of coverage"
if [[ ! -e "aligned_"$PROJECT_SQ".bam.bedgraph" ]]
then
    bedtools genomecov -ibam "sorted_aligned_"$PROJECT_SQ".bam" -bg > "aligned_"$PROJECT_SQ".bam.bedgraph"
fi

##########################################################
# Prepping for analysis 
##########################################################
echo "() Sorting BAM"

if [[ ! -e "sorted_aligned_"$PROJECT_SQ".bam" ]]
then
    java -jar $PICARD SortSam I="aligned_"$PROJECT_SQ".bam" O="sorted_aligned_"$PROJECT_SQ".bam" SORT_ORDER=coordinate
fi

echo "() Marking PCR Duplicates"
if [[ ! -e "marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam" ]]
then
    java -jar $PICARD MarkDuplicates I="sorted_aligned_"$PROJECT_SQ".bam" O="marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam" M=marked_dup_metrics.txt
fi

echo "() Generating recalibration table based on gnomAD"
if [[ ! -e "recal_data.table" ]]
then
    java -jar $GATK BaseRecalibrator -I "marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam" -R $REFERENCE --known-sites $GNOMAD_VCF -O recal_data.table
fi

echo "() Applying base recalibration"
if [[ ! -e "BQSR_marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam" ]]
then
    java -jar $GATK ApplyBQSR -R $REFERENCE -I "marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam" --bqsr-recal-file recal_data.table -O "BQSR_marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam"
fi

##########################################################
# Variant Calling
##########################################################
QUERY="BQSR_marked_duplicates_sorted_aligned_"$PROJECT_SQ".bam"
GATK_VCF="GATK_BQSR_sorted_marked_duplicates_aligned_"$PROJECT_SQ".bam.vcf.gz"

echo "() Variant calling with GATK - HaplotypeCaller"
if [[ ! -e $GATK_VCF ]]
then
    #java -jar $GATK --java-options "-Xmx4g" HaplotypeCaller -R $REFERENCE-I $QUERY -O $GATK_VCF -bamout GATK_bamout.bam
    java -jar $GATK HaplotypeCaller -R $REFERENCE -I $QUERY -O $GATK_VCF -bamout GATK_bamout.bam
fi

##########################################################
# VCF Statistics
##########################################################
bcftools stats $GATK_VCF > GATK_VCF.stats


##########################################################
# Annotate VCFs - VEP
##########################################################
#vep


echo "() PIPELINE COMPLETE"


##########################################################
# End of File.
##########################################################
