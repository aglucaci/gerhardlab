## REPOSITORY (RQ534361)
## PROJECT NAME - KA

## --- SOFTWARE USED --- ##
#BWA - Version: 0.7.17-r1188
#SAMTOOLS - Version: 1.7 (using htslib 1.7-2)
#PICARD - Version: 2.20.0
#GATK - Version: 4.1.1.0
#BCFTOOLS - Version:   v2.26.0
#BEDTOOLS - Version:   v2.26.0
#FASTQC - Version: 0.11.5

#Computer: this script is configured to run on the ubuntu-deploy server.
#java -version
#openjdk version "1.8.0_191"
#OpenJDK Runtime Environment (build 1.8.0_191-8u191-b12-2ubuntu0.18.04.1-b12)
#OpenJDK 64-Bit Server VM (build 25.191-b12, mixed mode)

#LINK TO REFERENCE - Illumina iGenomes
#https://support.illumina.com/sequencing/sequencing_software/igenome.html
#ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz

#LINK TO GNOMAD VCF - https://gnomad.broadinstitute.org/downloads
#https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

PICARD="/home/alexander/Downloads/picard.jar"
GATK="/home/alexander/Downloads/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar"
#cat genome_test.fa | perl -pe 's/chr//g' > genome_chr_replaced.fa
#Had to replace chr1 to 1 so that the gnomad VCF agrees with the reference, gnomAD uses hg19.

REFERENCE="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
REFERENCE_BWA

REFERENCE_BWA_INDEX="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.bwt"
REFERENCE_SAMTOOLS_INDEX="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
REFERENCE_DICT="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"
GNOMAD_VCF="/media/alexander/Elements/gnomad.exomes.r2.1.1.sites.vcf.bgz"
GNOMAD_VCF_INDEX="/media/alexander/Elements/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi"

READ1="SQ6981_S1_L00X_MASTER_R1_001.fastq.gz"
READ2="SQ6981_S1_L00X_MASTER_R2_001.fastq.gz"

WD="/media/alexander/Elements/RQ534361-KA/Scripts"

#######################################################################
# Initialize this script.
#######################################################################
clear
echo "Set WD:"
cd $WD
echo $WD
echo .
echo "## Initializing script - Whole Exome to Variant Calling (c) 2019"
echo "## by Alexander G. Lucaci"
echo "## Dataset: RQ534361-KA"
echo ""

#echo "Number of Steps:"

#######################################################################
# Data preprocessing
#######################################################################
echo "() Changing directory to ../Data"
cd ../Data

echo "() Combining FASTQ files.."
if [[ ! -e "../Analysis/Run_03/"$READ1 ]]
then
    cat SQ6981_S1_L001_R1_001.fastq.gz SQ6981_S1_L002_R1_001.fastq.gz SQ6981_S1_L003_R1_001.fastq.gz SQ6981_S1_L004_R1_001.fastq.gz > ../Analysis/Run_03/SQ6981_S1_L00X_MASTER_R1_001.fastq.gz
fi

if [[ ! -e "../Analysis/Run_03/"$READ2 ]]
then
    cat SQ6981_S1_L001_R2_001.fastq.gz SQ6981_S1_L002_R2_001.fastq.gz SQ6981_S1_L003_R2_001.fastq.gz SQ6981_S1_L004_R2_001.fastq.gz > ../Analysis/Run_03/SQ6981_S1_L00X_MASTER_R2_001.fastq.gz
fi

echo "() Changing directory to Data to ../Analysis/Run_03"
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
echo "() Indexing reference genome"
if [[ ! -e $REFERENCE_BWA_INDEX ]]
then
    bwa index $REFERENCE
    #continue
fi

echo "() Mapping reads to reference"
if [[ ! -e "aligned_SQ6981.sam" ]]
then
    bwa mem -M -t 7 -R "@RG\tID:SQ6981\tSM:SQ6981\tPL:Illumina" $REFERENCE $READ1 $READ2 > aligned_SQ6981.sam 
    #Don't need to store the SAM file: 
    #on Ubuntu this gave a memory error.
    #bwa mem -M -t 7 -R "@RG\tID:SQ6981\tSM:SQ6981\tPL:Illumina" $REFERENCE $READ1 $READ2 | samtools view -bS > aligned_SQ6981.bam
fi

echo "() Converting SAM to BAM"
if [[ ! -e "aligned_SQ6981.bam" ]]
then
    samtools view -bS aligned_SQ6981.sam > aligned_SQ6981.bam
fi

##########################################################
# Analyze the Mapping
##########################################################
echo "() Generating alignment flagstats"
if [[ ! -e "flagstat_aligned_SQ6981.bam.txt" ]]
then
    samtools flagstat aligned_SQ6981.bam > flagstat_aligned_SQ6981.bam.txt
fi

##########################################################
# Prepping for analysis 
##########################################################

echo "() Sorting BAM"
#could also do samtools sort
if [[ ! -e "sorted_aligned_SQ6981.bam" ]]
then
    java -jar $PICARD SortSam I=aligned_SQ6981.bam O=sorted_aligned_SQ6981.bam SORT_ORDER=coordinate
fi

echo "() BEDTools - Generate a bedgraph of coverage"
if [[ ! -e "aligned_SQ6981.bam.bedgraph" ]]
then
    bedtools genomecov -ibam sorted_aligned_SQ6981.bam -bg > aligned_SQ6981.bam.bedgraph
fi

echo "() Marking PCR Duplicates"
if [[ ! -e "marked_duplicates_sorted_aligned_SQ6981.bam" ]]
then
    java -jar $PICARD MarkDuplicates I=sorted_aligned_SQ6981.bam O=marked_duplicates_sorted_aligned_SQ6981.bam M=marked_dup_metrics.txt
fi

#echo "Cleaning up Reference for GATK"
#cat genome_test.fa | perl -pe 's/chr//g' > genome_chr_replaced.fa
#Had to replace chr1 to 1 so that the VCF agrees with the reference. gnomAD uses hg19.
#REFERENCE_CHR_REPLACED="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_chr_replaced.fa"
#echo "samtools faidx on the Reference chr replaced"
#samtools faidx $REFERENCE_CHR_REPLACED
#echo "Creating new dict"
#java -jar $PICARD CreateSequenceDictionary R=$REFERENCE_CHR_REPLACED O="/media/alexander/Elements/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome_chr_replaced.dict"

echo "() Indexing the gnomAD VCF"
if [[ ! -e $GNOMAD_VCF_INDEX ]]
then
    java -jar $GATK IndexFeatureFile -F $GNOMAD_VCF
    continue
fi

echo "() Samtools indexing the Reference fasta"
if [[ ! -e $REFERENCE_SAMTOOLS_INDEX ]]
then
    samtools faidx $REFERENCE
fi

echo "() Creating new dict"
if [[ ! -e $REFERENCE_DICT ]]
then
    java -jar $PICARD CreateSequenceDictionary R=$REFERENCE O=$REFERENCE_DICT
    #java -XX:ParallelGCThreads=<num of  threads> -jar <picard-package-name>.jar
fi

echo "() Generating recalibration table based on gnomAD"
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165
#https://www.biostars.org/p/312190/
if [[ ! -e "recal_data.table" ]]
then
    java -jar $GATK BaseRecalibrator -I marked_duplicates_sorted_aligned_SQ6981.bam -R $REFERENCE --known-sites $GNOMAD_VCF -O recal_data.table
fi

echo "() Applying base recalibration"
if [[ ! -e "BQSR_marked_duplicates_sorted_aligned_SQ6981.bam" ]]
then
    java -jar $GATK ApplyBQSR -R $REFERENCE -I marked_duplicates_sorted_aligned_SQ6981.bam --bqsr-recal-file recal_data.table -O BQSR_marked_duplicates_sorted_aligned_SQ6981.bam
fi

##########################################################
# Variant Calling
##########################################################
QUERY="BQSR_marked_duplicates_sorted_aligned_SQ6981.bam"

echo "() Variant calling with FreeBayes.."
freebayes -C 30 -f $REFERENCE $QUERY > Freebayes_30X_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.vcf

echo "() Variant calling with bcftools.."
#can switch to samtools mpileup, any difference?
bcftools mpileup --threads 6 -Ou -f $REFERENCE $QUERY | bcftools call --threads 7 -Ov -mv | bcftools filter -i 'QUAL>20 && DP>30' > bcftools_30X_aligned_markdup_SQ6981.mpileup.vcf

echo "() Variant calling with VarScan - SNPs.."
samtools mpileup -f $REFERENCE $QUERY | varscan mpileup2snp --min-coverage 30 --output-vcf 1 > varscan_snp_30X_samtools_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.mpileup.vcf

echo "() Variant calling with VarScan - InDels.."
samtools mpileup -f $REFERENCE $QUERY | varscan mpileup2indel --min-coverage 30 --output-vcf 1 > varscan_indel_30X_samtools_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.mpileup.vcf

echo "() Variant calling with GATK - HaplotypeCaller"
java -jar $GATK --java-options "-Xmx4g" HaplotypeCaller -R $REFERENCE-I $QUERY -O GATK_BQSR_sorted_marked_duplicates_aligned_SQ6981.bam.vcf.gz -bamout GATK_bamout.bam

##########################################################
# Consolidate VCFs
##########################################################
echo "() Concatenating VarScan output into one VCF (SNPs + Indels"
#bcftools concat -o 



##########################################################
# Intersect VCFs
##########################################################
#Where do they all agree?

#What is unique to each method?


##########################################################
# Annotate VCFs - VEP
##########################################################


##########################################################
# Filter VCF
##########################################################




echo "() PIPELINE COMPLETE"


##########################################################
# End of File.
##########################################################
