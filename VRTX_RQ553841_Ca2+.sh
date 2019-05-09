
FASTQC="/home/alexander/Downloads/FastQC/fastqc"
WD="/run/media/alexander/4TB-VD1/Projects/RQ553841_Ca2+/Scripts"
cd $WD

cd ../Data

## Files
#SQ7711_S1_L001_R1_001.fastq
#SQ7711_S1_L001_R2_001.fastq
#SQ7711_S1_L002_R1_001.fastq
#SQ7711_S1_L002_R2_001.fastq


#sha256sum SQ7711_S1_L001_R1_001.fastq.gz
#sha256sum SQ7711_S1_L001_R2_001.fastq.gz
#sha256sum SQ7711_S1_L002_R1_001.fastq.gz
#sha256sum SQ7711_S1_L002_R2_001.fastq.gz


$FASTQC SQ7711_S1_L001_R1_001.fastq.gz
$FASTQC SQ7711_S1_L001_R2_001.fastq.gz
$FASTQC SQ7711_S1_L002_R1_001.fastq.gz
$FASTQC SQ7711_S1_L002_R2_001.fastq.gz
