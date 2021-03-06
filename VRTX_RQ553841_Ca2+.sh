
## REPOSITORY (RQ553841)
## PROJECT NAME - Ca2+

#LINK TO REFERENCE - Illumina iGenomes
#https://support.illumina.com/sequencing/sequencing_software/igenome.html
#ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz

#LINK TO GNOMAD VCF - https://gnomad.broadinstitute.org/downloads
#https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

## Software Versions
#FASTQC
#PICARD
#GATK 4.1.2.0
#samtools Version: 0.1.19-44428cd

## -- SOFTWARE -- ##
FASTQC="/home/alexander/Downloads/FastQC/fastqc"
PICARD="/home/alexander/Downloads/picard.jar"
GATK="/home/alexander/Downloads/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"

## Reference files ##
REFERENCE="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
REFERENCE_INDEX="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
REFERENCE_DICT="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"
REFERENCE_BWA="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"

## gnomAD Annotation ##
GNOMAD_VCF="/run/media/alexander/4TB-VD1/Projects/GNOMAD/gnomad.exomes.r2.1.1.sites.vcf.bgz"
GNOMAD_INDEX="/run/media/alexander/4TB-VD1/Projects/GNOMAD/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi"

##Ammended gnomAD Annotation ##
#Change 1 to chr1.
GNOMAD_AMMENDED_VCF="/run/media/alexander/4TB-VD1/Projects/GNOMAD/AMMENDED_gnomad.exomes.r2.1.1.sites.vcf.bgz"
GNOMAD_AMMENDED_INDEX="/run/media/alexander/4TB-VD1/Projects/GNOMAD/AMMENDED_gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi"

#Set Working Directory
WD="/run/media/alexander/4TB-VD1/Projects/RQ553841_Ca2+/Scripts"
cd $WD

#CAT the fastq files together into forward and reverse
READ1="SQ7711_S1_L00X_COMBINE_R1_001.fastq.gz"
READ2="SQ7711_S1_L00X_COMBINE_R2_001.fastq.gz"

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
# Data QC
#######################################################################

## -- Files
#SQ7711_S1_L001_R1_001.fastq.gz
#SQ7711_S1_L001_R2_001.fastq.gz
#SQ7711_S1_L002_R1_001.fastq.gz
#SQ7711_S1_L002_R2_001.fastq.gz

# -- SHA CHECK
#sha256sum SQ7711_S1_L001_R1_001.fastq.gz
#sha256sum SQ7711_S1_L001_R2_001.fastq.gz
#sha256sum SQ7711_S1_L002_R1_001.fastq.gz
#sha256sum SQ7711_S1_L002_R2_001.fastq.gz

# -- FASTQC
#$FASTQC SQ7711_S1_L001_R1_001.fastq.gz
#$FASTQC SQ7711_S1_L001_R2_001.fastq.gz
#$FASTQC SQ7711_S1_L002_R1_001.fastq.gz
#$FASTQC SQ7711_S1_L002_R2_001.fastq.gz

#######################################################################
# Data Preprocessing
#######################################################################
echo "() Changing directory to ../Data"
cd ../Data

echo "() Combining RQ553841 FASTQ files.."
echo "() Combining forward reads"
if [[ ! -e "../Analysis/"$READ1 ]]
then
    cat SQ7711_S1_L001_R1_001.fastq.gz SQ7711_S1_L002_R1_001.fastq.gz > ../Analysis/$READ1
fi

echo "() Combining reverse reads"
if [[ ! -e "../Analysis/"$READ2 ]]
then
    cat SQ7711_S1_L001_R2_001.fastq.gz SQ7711_S1_L002_R2_001.fastq.gz > ../Analysis/$READ2
fi

echo "() Changing directory to Data to ../Analysis"
cd ../Analysis

##########################################################
# Mapping 
##########################################################
echo "() Mapping reads to reference"

if [[ ! -e "aligned_SQ7711.sam" ]]
then
    bwa mem -M -t 25 -R "@RG\tID:SQ7711\tSM:SQ7711\tPL:Illumina" $REFERENCE_BWA $READ1 $READ2 > aligned_SQ7711.sam 
fi

echo "() Converting SAM to BAM"
if [[ ! -e "aligned_SQ7711.bam" ]]
then
    samtools view -bS aligned_SQ7711.sam > aligned_SQ7711.bam
fi

##########################################################
# Analyze the Mapping
##########################################################
echo "() Generating alignment flagstats"
if [[ ! -e "flagstat_aligned_SQ7711.bam.txt" ]]
then
    samtools flagstat aligned_SQ7711.bam > flagstat_aligned_SQ7711.bam.txt
fi

##########################################################
# Prepping for analysis 
##########################################################
echo "() Sorting BAM"

if [[ ! -e "sorted_aligned_SQ7711.bam" ]]
then
    java -jar $PICARD SortSam I=aligned_SQ7711.bam O=sorted_aligned_SQ7711.bam SORT_ORDER=coordinate
    #samtools sort aligned_SQ7711.bam > sorted_aligned_SQ7711.bam
fi

#echo "BEDTools - Generate a bedgraph of coverage"
#if [[ ! -e "aligned_SQ6981.bam.bedgraph" ]]
#then
#    bedtools genomecov -ibam sorted_aligned_SQ6981.bam -bg > aligned_SQ6981.bam.bedgraph
#fi

echo "Marking PCR Duplicates"
if [[ ! -e "marked_duplicates_sorted_aligned_SQ7711.bam" ]]
then
    java -jar $PICARD MarkDuplicates I=sorted_aligned_SQ7711.bam O=marked_duplicates_sorted_aligned_SQ7711.bam M=marked_dup_metrics.txt
fi

#Ammend gnomAD VCF
echo "Ammending gnomAD VCF - for agreement"
if [[ ! -e $GNOMAD_AMMENDED_VCF ]]
then
    python ../Scripts/ammend_gnomad.py > $GNOMAD_AMMENDED_VCF
fi

#Index the ammended gnomAD VCF
echo "Indexing the Ammended gnomAD VCF"
if [[ ! -e $GNOMAD_AMMENDED_INDEX ]]
then
    java -jar $GATK IndexFeatureFile -F $GNOMAD_AMMENDED_INDEX
    continue
fi

echo "Generating recalibration table based on gnomAD"
if [[ ! -e "recal_data.table" ]]
then
    java -jar $GATK BaseRecalibrator -I marked_duplicates_sorted_aligned_SQ7711.bam -R $REFERENCE --known-sites $GNOMAD_AMMENDED_VCF -O recal_data.table
fi

echo "Applying base recalibration"
if [[ ! -e "BQSR_marked_duplicates_sorted_aligned_SQ6981.bam" ]]
then
    java -jar $GATK ApplyBQSR -R $REFERENCE -I marked_duplicates_sorted_aligned_SQ7711.bam --bqsr-recal-file recal_data.table -O BQSR_marked_duplicates_sorted_aligned_SQ7711.bam
fi

##########################################################
# Variant Calling
##########################################################
QUERY="BQSR_marked_duplicates_sorted_aligned_SQ7711.bam"

#######################################################################
# End of file
#######################################################################



