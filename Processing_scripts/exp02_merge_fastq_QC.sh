mkdir Fastqs_merged
cat Fastqs/A1WT1ku1_activated_B26_17_L00*/* > Fastqs_merged/A1WT1ku1_activated_B26_17.fastq.gz
cat Fastqs/A2WT1ku1_activated_B26_18_L00*/* > Fastqs_merged/A2WT1ku1_activated_B26_18.fastq.gz
cat Fastqs/A3WT1ku1_activated_B26_19_L00*/* > Fastqs_merged/A3WT1ku1_activated_B26_19.fastq.gz
cat Fastqs/B1WT1ku1_activatedplusp38i_B26_17*/* > Fastqs_merged/B1WT1ku1_activatedplusp38i_B26_17.fastq.gz
cat Fastqs/B2WT1ku1_activatedplusp38i_B26-18*/* > Fastqs_merged/B2WT1ku1_activatedplusp38i_B26_18.fastq.gz
cat Fastqs/B3WT1ku1_activatedplusp38i_B26_19*/* > Fastqs_merged/B3WT1ku1_activatedplusp38i_B26_19.fastq.gz
cat Fastqs/C1GTKO1ku1_activated_B26_17_L00*/* > Fastqs_merged/C1GTKO1ku1_activated_B26_17.fastq.gz
cat Fastqs/C2GTKO1ku1_activated_B26_18_L00*/* > Fastqs_merged/C2GTKO1ku1_activated_B26_18.fastq.gz
cat Fastqs/C3GTKO1ku1_activated_B26_19_L00*/* > Fastqs_merged/C3GTKO1ku1_activated_B26_19.fastq.gz
cat Fastqs/D1GTKO_1ku1_activatedplusp38i_B26_17*/* >  Fastqs_merged/D1GTKO1ku1_activatedplusp38i_B26_17.fastq.gz
cat Fastqs/D2GTKO1ku1_activatedplusp38i_B26_18*/* > Fastqs_merged/D2GTKO1ku1_activatedplusp38i_B26_18.fastq.gz
cat Fastqs/D3GTKO1ku1_activatedplusp38i_B26_19*/* > Fastqs_merged/D3GTKO1ku1_activatedplusp38i_B26_19.fastq.gz

# QC analysis of merged data
mkdir FastQC
fastqc -o FastQC Fastqs_merged/*

# generating multiQC from above
mkdir MultiQC
cd MultiQC
multiqc ../FastQC/