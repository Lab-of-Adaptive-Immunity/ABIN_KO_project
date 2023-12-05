# Author: Juraj Michalik, Lab of Adaptive Immunity, Institute of Molecular Genetics of the Czech Academy of Sciences
# LICENCE: MIT (see LICENCE.md)

mkdir STAR2 
mkdir STAR2/A1WT1ku1_activated_B26_17
mkdir STAR2/A2WT1ku1_activated_B26_18
mkdir STAR2/A3WT1ku1_activated_B26_19
mkdir STAR2/B1WT1ku1_activatedplusp38i_B26_17
mkdir STAR2/B2WT1ku1_activatedplusp38i_B26_18
mkdir STAR2/B3WT1ku1_activatedplusp38i_B26_19
mkdir STAR2/C1GTKO1ku1_activated_B26_17
mkdir STAR2/C2GTKO1ku1_activated_B26_18
mkdir STAR2/C3GTKO1ku1_activated_B26_19
mkdir STAR2/D1GTKO1ku1_activatedplusp38i_B26_17
mkdir STAR2/D2GTKO1ku1_activatedplusp38i_B26_18
mkdir STAR2/D3GTKO1ku1_activatedplusp38i_B26_19

STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/A1WT1ku1_activated_B26_17.fastq.gz --outFileNamePrefix STAR2/A1WT1ku1_activated_B26_17/A1WT1ku1_activated_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/A2WT1ku1_activated_B26_18.fastq.gz --outFileNamePrefix STAR2/A2WT1ku1_activated_B26_18/A2WT1ku1_activated_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/A3WT1ku1_activated_B26_19.fastq.gz --outFileNamePrefix STAR2/A3WT1ku1_activated_B26_19/A3WT1ku1_activated_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/B1WT1ku1_activatedplusp38i_B26_17.fastq.gz --outFileNamePrefix STAR2/B1WT1ku1_activatedplusp38i_B26_17/B1WT1ku1_activatedplusp38i_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/B2WT1ku1_activatedplusp38i_B26_18.fastq.gz --outFileNamePrefix STAR2/B2WT1ku1_activatedplusp38i_B26_18/B2WT1ku1_activatedplusp38i_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/B3WT1ku1_activatedplusp38i_B26_19.fastq.gz --outFileNamePrefix STAR2/B3WT1ku1_activatedplusp38i_B26_19/B3WT1ku1_activatedplusp38i_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/C1GTKO1ku1_activated_B26_17.fastq.gz --outFileNamePrefix STAR2/C1GTKO1ku1_activated_B26_17/C1GTKO1ku1_activated_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/C2GTKO1ku1_activated_B26_18.fastq.gz --outFileNamePrefix STAR2/C2GTKO1ku1_activated_B26_18/C2GTKO1ku1_activated_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/C3GTKO1ku1_activated_B26_19.fastq.gz --outFileNamePrefix STAR2/C3GTKO1ku1_activated_B26_19/C3GTKO1ku1_activated_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/D1GTKO1ku1_activatedplusp38i_B26_17.fastq.gz --outFileNamePrefix STAR2/D1GTKO1ku1_activatedplusp38i_B26_17/D1GTKO1ku1_activatedplusp38i_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/D2GTKO1ku1_activatedplusp38i_B26_18.fastq.gz --outFileNamePrefix STAR2/D2GTKO1ku1_activatedplusp38i_B26_18/D2GTKO1ku1_activatedplusp38i_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir GRCm39_106/ --readFilesIn Fastqs_merged2/D3GTKO1ku1_activatedplusp38i_B26_19.fastq.gz --outFileNamePrefix STAR2/D3GTKO1ku1_activatedplusp38i_B26_19/D3GTKO1ku1_activatedplusp38i_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

# copy bams into new dir + index them
mkdir Bamfiles2
cp STAR2/*/*sortedByCoord.out.bam Bamfiles2
cd Bamfiles2
ls *.bam | xargs -n1 -P5 samtools index
