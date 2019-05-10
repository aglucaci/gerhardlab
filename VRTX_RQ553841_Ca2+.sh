
## REPOSITORY (RQ553841)
## PROJECT NAME - Ca2+

#LINK TO REFERENCE - Illumina iGenomes
#https://support.illumina.com/sequencing/sequencing_software/igenome.html
#ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz

#LINK TO GNOMAD VCF - https://gnomad.broadinstitute.org/downloads
#https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

FASTQC="/home/alexander/Downloads/FastQC/fastqc"
REFERENCE="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
REFERENCE_INDEX="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai"
REFERENCE_DICT="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.dict"
REFERENCE_BWA="/run/media/alexander/4TB-VD1/Projects/REFERENCE/Homo_sapiens_UCSC_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
WD="/run/media/alexander/4TB-VD1/Projects/RQ553841_Ca2+/Scripts"

READ1="SQ7711_S1_L00X_COMBINE_R1_001.fastq.gz"
READ2="SQ7711_S1_L00X_COMBINE_R2_001.fastq.gz"

cd $WD

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
echo "() Read 1"
if [[ ! -e "../Analysis/"$READ1 ]]
then
    cat SQ7711_S1_L001_R1_001.fastq.gz SQ7711_S1_L002_R1_001.fastq.gz > ../Analysis/$READ1
fi

echo "() Read 2"
if [[ ! -e "../Analysis/"$READ2 ]]
then
    cat SQ7711_S1_L001_R2_001.fastq.gz SQ7711_S1_L002_R2_001.fastq.gz > ../Analysis/$READ2
fi

echo "Changing directory to Data to ../Analysis"
cd ../Analysis

##########################################################
# Mapping 
##########################################################
echo "Mapping reads to reference"

if [[ ! -e "aligned_SQ7711.sam" ]]
then
    bwa mem -M -t 25 -R "@RG\tID:SQ7711\tSM:SQ7711\tPL:Illumina" $REFERENCE_BWA $READ1 $READ2 > aligned_SQ7711.sam 
    #Don't need to store the SAM file: 
    #on Ubuntu this gave a memory error.
    #VRTX COMMAND
    #bwa mem -M -t 25 -R "@RG\tID:SQ7711\tSM:SQ7711\tPL:Illumina" $REFERENCE_BWA $READ1 $READ2 | samtools view -bS > aligned_SQ7711.bam
fi

#exit 1

#######################################################################
# End of file
#######################################################################



