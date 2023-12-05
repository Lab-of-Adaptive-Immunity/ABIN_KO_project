# Author: Juraj Michalik, Lab of Adaptive Immunity, Institute of Molecular Genetics of the Czech Academy of Sciences
# LICENCE: MIT (see LICENCE.md)

mkdir STAR1
mkdir STAR1/GTKO_1ku1_B26-6
mkdir STAR1/GTKO_1ku1_B26-7
mkdir STAR1/GTKO_1ku1_B26-8
mkdir STAR1/GTKO_4ku1_B26-6
mkdir STAR1/GTKO_4ku1_B26-7
mkdir STAR1/GTKO_4ku1_B26-8
mkdir STAR1/GTKO_non-act_B26-7
mkdir STAR1/GTKO_non-act_B26-8
mkdir STAR1/GTKO_non-act_B26-9
mkdir STAR1/WT_1ku1_B26-6
mkdir STAR1/WT_1ku1_B26-7
mkdir STAR1/WT_1ku1_B26-8
mkdir STAR1/WT_4ku1_B26-6
mkdir STAR1/WT_4ku1_B26-7
mkdir STAR1/WT_4ku1_B26-8
mkdir STAR1/WT_non-act_B26-7
mkdir STAR1/WT_non-act_B26-8
mkdir STAR1/WT_non-act_B26-9

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_1ku1_B26-6.fastq.gz --outFileNamePrefix STAR1/GTKO_1ku1_B26-6/GTKO_1ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_1ku1_B26-7.fastq.gz --outFileNamePrefix STAR1/GTKO_1ku1_B26-7/GTKO_1ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_1ku1_B26-8.fastq.gz --outFileNamePrefix STAR1/GTKO_1ku1_B26-8/GTKO_1ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_4ku1_B26-6.fastq.gz --outFileNamePrefix STAR1/GTKO_4ku1_B26-6/GTKO_4ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_4ku1_B26-7.fastq.gz --outFileNamePrefix STAR1/GTKO_4ku1_B26-7/GTKO_4ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_4ku1_B26-8.fastq.gz --outFileNamePrefix STAR1/GTKO_4ku1_B26-8/GTKO_4ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_non-act_B26-7.fastq.gz --outFileNamePrefix STAR1/GTKO_non-act_B26-7/GTKO_non-act_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_non-act_B26-8.fastq.gz --outFileNamePrefix STAR1/GTKO_non-act_B26-8/GTKO_non-act_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/GTKO_non-act_B26-9.fastq.gz --outFileNamePrefix STAR1/GTKO_non-act_B26-9/GTKO_non-act_B26-9 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_1ku1_B26-6.fastq.gz --outFileNamePrefix STAR1/WT_1ku1_B26-6/WT_1ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_1ku1_B26-7.fastq.gz --outFileNamePrefix STAR1/WT_1ku1_B26-7/WT_1ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_1ku1_B26-8.fastq.gz --outFileNamePrefix STAR1/WT_1ku1_B26-8/WT_1ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_4ku1_B26-6.fastq.gz --outFileNamePrefix STAR1/WT_4ku1_B26-6/WT_4ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_4ku1_B26-7.fastq.gz --outFileNamePrefix STAR1/WT_4ku1_B26-7/WT_4ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_4ku1_B26-8.fastq.gz --outFileNamePrefix STAR1/WT_4ku1_B26-8/WT_4ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_non-act_B26-7.fastq.gz --outFileNamePrefix STAR1/WT_non-act_B26-7/WT_non-act_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_non-act_B26-8.fastq.gz --outFileNamePrefix STAR1/WT_non-act_B26-8/WT_non-act_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged1/WT_non-act_B26-9.fastq.gz --outFileNamePrefix STAR1/WT_non-act_B26-9/WT_non-act_B26-9 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

# copy bams into new dir + index them
mkdir Bamfiles1
cp STAR1/*/*sortedByCoord.out.bam Bamfiles1
cd Bamfiles1
ls *.bam | xargs -n1 -P5 samtools index
