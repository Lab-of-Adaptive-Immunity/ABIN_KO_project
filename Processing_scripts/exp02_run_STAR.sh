mkdir STAR 
mkdir STAR/A1WT1ku1_activated_B26_17
mkdir STAR/A2WT1ku1_activated_B26_18
mkdir STAR/A3WT1ku1_activated_B26_19
mkdir STAR/B1WT1ku1_activatedplusp38i_B26_17
mkdir STAR/B2WT1ku1_activatedplusp38i_B26_18
mkdir STAR/B3WT1ku1_activatedplusp38i_B26_19
mkdir STAR/C1GTKO1ku1_activated_B26_17
mkdir STAR/C2GTKO1ku1_activated_B26_18
mkdir STAR/C3GTKO1ku1_activated_B26_19
mkdir STAR/D1GTKO1ku1_activatedplusp38i_B26_17
mkdir STAR/D2GTKO1ku1_activatedplusp38i_B26_18
mkdir STAR/D3GTKO1ku1_activatedplusp38i_B26_19

STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/A1WT1ku1_activated_B26_17.fastq.gz --outFileNamePrefix STAR/A1WT1ku1_activated_B26_17/A1WT1ku1_activated_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/A2WT1ku1_activated_B26_18.fastq.gz --outFileNamePrefix STAR/A2WT1ku1_activated_B26_18/A2WT1ku1_activated_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/A3WT1ku1_activated_B26_19.fastq.gz --outFileNamePrefix STAR/A3WT1ku1_activated_B26_19/A3WT1ku1_activated_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/B1WT1ku1_activatedplusp38i_B26_17.fastq.gz --outFileNamePrefix STAR/B1WT1ku1_activatedplusp38i_B26_17/B1WT1ku1_activatedplusp38i_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/B2WT1ku1_activatedplusp38i_B26_18.fastq.gz --outFileNamePrefix STAR/B2WT1ku1_activatedplusp38i_B26_18/B2WT1ku1_activatedplusp38i_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/B3WT1ku1_activatedplusp38i_B26_19.fastq.gz --outFileNamePrefix STAR/B3WT1ku1_activatedplusp38i_B26_19/B3WT1ku1_activatedplusp38i_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/C1GTKO1ku1_activated_B26_17.fastq.gz --outFileNamePrefix STAR/C1GTKO1ku1_activated_B26_17/C1GTKO1ku1_activated_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/C2GTKO1ku1_activated_B26_18.fastq.gz --outFileNamePrefix STAR/C2GTKO1ku1_activated_B26_18/C2GTKO1ku1_activated_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/C3GTKO1ku1_activated_B26_19.fastq.gz --outFileNamePrefix STAR/C3GTKO1ku1_activated_B26_19/C3GTKO1ku1_activated_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/D1GTKO1ku1_activatedplusp38i_B26_17.fastq.gz --outFileNamePrefix STAR/D1GTKO1ku1_activatedplusp38i_B26_17/D1GTKO1ku1_activatedplusp38i_B26_17 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/D2GTKO1ku1_activatedplusp38i_B26_18.fastq.gz --outFileNamePrefix STAR/D2GTKO1ku1_activatedplusp38i_B26_18/D2GTKO1ku1_activatedplusp38i_B26_18 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 10 --genomeDir /mnt/scratch/beluga/michalik/STAR_references_Conda/GRCm39_106/ --readFilesIn Fastqs_merged/D3GTKO1ku1_activatedplusp38i_B26_19.fastq.gz --outFileNamePrefix STAR/D3GTKO1ku1_activatedplusp38i_B26_19/D3GTKO1ku1_activatedplusp38i_B26_19 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

# copy bams into new dir + index them
mkdir Bamfiles_X
cp STAR/*/*sortedByCoord.out.bam Bamfiles_X
cd Bamfiles_X
ls *.bam | xargs -n1 -P5 samtools index
