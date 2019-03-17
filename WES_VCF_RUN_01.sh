##Alexander Lucaci
##Temple University
##Gerhard Lab, Department of Medical Genetics and Molecular Biochemistry

#Bulk of calculations done on bioinformatics server, bioinform.cst.temple.edu

#Initial WES to VCF runthrough. (Experiment #01)
#Protocol as follows

#GATHER INITIAL INFO ON FASTQ files using FASTQC
fastqc ../Data/SQ8992_L001_R1_001.fastq  (all the files)

#Reference genome from GRC
#Build the index
bowtie2-build -f hg38.fa hg38_index

#Map our reads, with SAM file output
bowtie2 -x hg38_index -1 SQ8992_L001_R1_001.fastq -2 SQ8992_L001_R2_001.fastq > SQ8992_hg38.sam

#Modifiy the name, this will be excluded in the next run.
Cp SQ8992_hg38.sam SQ8992_L001_R1_R2_hg38.sam

#SAM to BAM
samtools view -bS SQ8992_L001_R1_R2_hg38.sam >  SQ8992_L001_R1_R2_hg38.bam

#Sort our BAM file
samtools sort SQ8992_L001_R1_R2_hg38.bam > SQ8992_L001_R1_R2_hg38_sorted.bam

#Index our BAM
samtools index SQ8992_L001_R1_R2_hg38_sorted.bam

#Pileup our reads, will probably be switched to bcftools mpileup in follow on experiment.
samtools mpileup -E -uf hg38.fa SQ8992_L001_R1_R2_hg38_sorted.bam > SQ8992_L001_R1_R2_hg38_sorted.mpileup

#Make variant calls.
bcftools call -v -m SQ8992_L001_R1_R2_hg38_sorted.mpileup > SQ8992_L001_R1_R2_hg38_sorted.mpileup_variants.vcf

#Output is VCF, annotation done with VEP, check results with gnomAD.
#https://useast.ensembl.org/Homo_sapiens/Tools/VEP/
#https://gnomad.broadinstitute.org/

