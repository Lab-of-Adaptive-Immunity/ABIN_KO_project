mkdir Fastqs_merged1
cat Fastqs/GTKO_1ku1_B26_6* > Fastqs_merged1/GTKO_1ku1_B26-6.fastq.gz
cat Fastqs/GTKO_1ku1_B26_7* > Fastqs_merged1/GTKO_1ku1_B26-7.fastq.gz
cat Fastqs/GTKO_1ku1_B26_8* > Fastqs_merged1/GTKO_1ku1_B26-8.fastq.gz
cat Fastqs/GTKO_4ku1_B26_6* > Fastqs_merged1/GTKO_4ku1_B26-6.fastq.gz
cat Fastqs/GTKO_4ku1_B26_7* > Fastqs_merged1/GTKO_4ku1_B26-7.fastq.gz
cat Fastqs/GTKO_4ku1_B26_8* > Fastqs_merged1/GTKO_4ku1_B26-8.fastq.gz
cat Fastqs/GTKO_non_act_B26_7* > Fastqs_merged1/GTKO_non-act_B26-7.fastq.gz
cat Fastqs/GTKO_non_act_B26_8* > Fastqs_merged1/GTKO_non-act_B26-8.fastq.gz
cat Fastqs/GTKO_non_act_B26_9* > Fastqs_merged1/GTKO_non-act_B26-9.fastq.gz
cat Fastqs/WT_1ku1_B26_6* > Fastqs_merged1/WT_1ku1_B26-6.fastq.gz
cat Fastqs/WT_1ku1_B26_7* > Fastqs_merged1/WT_1ku1_B26-7.fastq.gz
cat Fastqs/WT_1ku1_B26_8* > Fastqs_merged1/WT_1ku1_B26-8.fastq.gz
cat Fastqs/WT_4ku1_B26_6* > Fastqs_merged1/WT_4ku1_B26-6.fastq.gz
cat Fastqs/WT_4ku1_B26_7* > Fastqs_merged1/WT_4ku1_B26-7.fastq.gz
cat Fastqs/WT_4ku1_B26_8* > Fastqs_merged1/WT_4ku1_B26-8.fastq.gz
cat Fastqs/WT_non_act_B26_7* > Fastqs_merged1/WT_non-act_B26-7.fastq.gz
cat Fastqs/WT_non_act_B26_8* > Fastqs_merged1/WT_non-act_B26-8.fastq.gz
cat Fastqs/WT_non_act_B26_9* > Fastqs_merged1/WT_non-act_B26-9.fastq.gz

# QC analysis of merged data
mkdir FastQC1
fastqc -o FastQC1 Fastqs_merged1/*

# generating multiQC from above
mkdir MultiQC1
cd MultiQC1
multiqc ../FastQC1/
