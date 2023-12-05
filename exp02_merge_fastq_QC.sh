# Author: Juraj Michalik, Lab of Adaptive Immunity, Institute of Molecular Genetics of the Czech Academy of Sciences
# LICENCE: MIT (see LICENCE.md)

mkdir Fastqs_merged2
cat Fastqs/A1WT1ku1_activated_B26_17* > Fastqs_merged2/A1WT1ku1_activated_B26_17.fastq.gz
cat Fastqs/A2WT1ku1_activated_B26_18* > Fastqs_merged2/A2WT1ku1_activated_B26_18.fastq.gz
cat Fastqs/A3WT1ku1_activated_B26_19* > Fastqs_merged2/A3WT1ku1_activated_B26_19.fastq.gz
cat Fastqs/B1WT1ku1_activatedplusp38i_B26_17* > Fastqs_merged2/B1WT1ku1_activatedplusp38i_B26_17.fastq.gz
cat Fastqs/B2WT1ku1_activatedplusp38i_B26_18* > Fastqs_merged2/B2WT1ku1_activatedplusp38i_B26_18.fastq.gz
cat Fastqs/B3WT1ku1_activatedplusp38i_B26_19* > Fastqs_merged2/B3WT1ku1_activatedplusp38i_B26_19.fastq.gz
cat Fastqs/C1GTKO1ku1_activated_B26_17* > Fastqs_merged2/C1GTKO1ku1_activated_B26_17.fastq.gz
cat Fastqs/C2GTKO1ku1_activated_B26_18* > Fastqs_merged2/C2GTKO1ku1_activated_B26_18.fastq.gz
cat Fastqs/C3GTKO1ku1_activated_B26_19* > Fastqs_merged2/C3GTKO1ku1_activated_B26_19.fastq.gz
cat Fastqs/D1GTKO_1ku1_activatedplusp38i_B26_17* >  Fastqs_merged2/D1GTKO1ku1_activatedplusp38i_B26_17.fastq.gz
cat Fastqs/D2GTKO1ku1_activatedplusp38i_B26_18* > Fastqs_merged2/D2GTKO1ku1_activatedplusp38i_B26_18.fastq.gz
cat Fastqs/D3GTKO1ku1_activatedplusp38i_B26_19* > Fastqs_merged2/D3GTKO1ku1_activatedplusp38i_B26_19.fastq.gz

# QC analysis of merged data
mkdir FastQC2
fastqc -o FastQC2 Fastqs_merged2/*

# generating multiQC from above
mkdir MultiQC2
cd MultiQC2
multiqc ../FastQC2/
