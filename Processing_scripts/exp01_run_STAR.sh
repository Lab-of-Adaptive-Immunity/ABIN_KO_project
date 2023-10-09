mkdir STAR
mkdir STAR/GTKO_1ku1_B26_6
mkdir STAR/GTKO_1ku1_B26_7
mkdir STAR/GTKO_1ku1_B26_8
mkdir STAR/GTKO_4ku1_B26_6
mkdir STAR/GTKO_4ku1_B26_7
mkdir STAR/GTKO_4ku1_B26_8
mkdir STAR/GTKO_non-act_B26_7
mkdir STAR/GTKO_non-act_B26_8
mkdir STAR/GTKO_non-act_B26_9
mkdir STAR/WT_1ku1_B26_6
mkdir STAR/WT_1ku1_B26_7
mkdir STAR/WT_1ku1_B26_8
mkdir STAR/WT_4ku1_B26_6
mkdir STAR/WT_4ku1_B26_7
mkdir STAR/WT_4ku1_B26_8
mkdir STAR/WT_non-act_B26_7
mkdir STAR/WT_non-act_B26_8
mkdir STAR/WT_non-act_B26_9

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_1ku1_B26-6.fastq.gz --outFileNamePrefix STAR/GTKO_1ku1_B26_6/GTKO_1ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_1ku1_B26-7.fastq.gz --outFileNamePrefix STAR/GTKO_1ku1_B26_7/GTKO_1ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_1ku1_B26-8.fastq.gz --outFileNamePrefix STAR/GTKO_1ku1_B26_8/GTKO_1ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_4ku1_B26-6.fastq.gz --outFileNamePrefix STAR/GTKO_4ku1_B26_6/GTKO_4ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_4ku1_B26-7.fastq.gz --outFileNamePrefix STAR/GTKO_4ku1_B26_7/GTKO_4ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_4ku1_B26-8.fastq.gz --outFileNamePrefix STAR/GTKO_4ku1_B26_8/GTKO_4ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_non-act_B26-7.fastq.gz --outFileNamePrefix STAR/GTKO_non-act_B26-7/GTKO_non-act_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_non-act_B26-8.fastq.gz --outFileNamePrefix STAR/GTKO_non-act_B26-8/GTKO_non-act_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/GTKO_non-act_B26-9.fastq.gz --outFileNamePrefix STAR/GTKO_non-act_B26-9/GTKO_non-act_B26-9 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_1ku1_B26-6.fastq.gz --outFileNamePrefix STAR/WT_1ku1_B26_6/WT_1ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_1ku1_B26-7.fastq.gz --outFileNamePrefix STAR/WT_1ku1_B26_7/WT_1ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_1ku1_B26-8.fastq.gz --outFileNamePrefix STAR/WT_1ku1_B26_8/WT_1ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_4ku1_B26-6.fastq.gz --outFileNamePrefix STAR/WT_4ku1_B26_6/WT_4ku1_B26-6 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_4ku1_B26-7.fastq.gz --outFileNamePrefix STAR/WT_4ku1_B26_7/WT_4ku1_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_4ku1_B26-8.fastq.gz --outFileNamePrefix STAR/WT_4ku1_B26_8/WT_4ku1_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_non-act_B26-7.fastq.gz --outFileNamePrefix STAR/WT_non-act_B26_6/WT_non-act_B26-7 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_non-act_B26-8.fastq.gz --outFileNamePrefix STAR/WT_non-act_B26_7/WT_non-act_B26-8 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/WT_non-act_B26-9.fastq.gz --outFileNamePrefix STAR/WT_non-act_B26_8/WT_non-act_B26-9 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

# copy bams into new dir + index them
mkdir Bamfiles
cp STAR/*/*sortedByCoord.out.bam Bamfiles
cd Bamfiles
ls *.bam | xargs -n1 -P5 samtools index
