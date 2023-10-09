mkdir Fastqs_merged
cat Fastqs/GTKO_1ku1_B26-6_L00*/* > Fastqs_merged/GTKO_1ku1_B26-6.fastq.gz
cat Fastqs/GTKO_1ku1_B26-7_L00*/* > Fastqs_merged/GTKO_1ku1_B26-7.fastq.gz
cat Fastqs/GTKO_1ku1_B26-8_L00*/* > Fastqs_merged/GTKO_1ku1_B26-8.fastq.gz
cat Fastqs/GTKO_4ku1_B26-6_L00*/* > Fastqs_merged/GTKO_4ku1_B26-6.fastq.gz
cat Fastqs/GTKO_4ku1_B26-7_L00*/* > Fastqs_merged/GTKO_4ku1_B26-7.fastq.gz
cat Fastqs/GTKO_4ku1_B26-8_L00*/* > Fastqs_merged/GTKO_4ku1_B26-8.fastq.gz
cat Fastqs/GTKO_non-act_B26-7_L00*/* > Fastqs_merged/GTKO_non-act_B26-7.fastq.gz
cat Fastqs/GTKO_non-act_B26-8_L00*/* > Fastqs_merged/GTKO_non-act_B26-8.fastq.gz
cat Fastqs/GTKO_non-act_B26-9_L00*/* > Fastqs_merged/GTKO_non-act_B26-9.fastq.gz
cat Fastqs/WT_1ku1_B26-6_L00*/* > Fastqs_merged/WT_1ku1_B26-6.fastq.gz
cat Fastqs/WT_1ku1_B26-7_L00*/* > Fastqs_merged/WT_1ku1_B26-7.fastq.gz
cat Fastqs/WT_1ku1_B26-8_L00*/* > Fastqs_merged/WT_1ku1_B26-8.fastq.gz
cat Fastqs/WT_4ku1_B26-6_L00*/* > Fastqs_merged/WT_4ku1_B26-6.fastq.gz
cat Fastqs/WT_4ku1_B26-7_L00*/* > Fastqs_merged/WT_4ku1_B26-7.fastq.gz
cat Fastqs/WT_4ku1_B26-8_L00*/* > Fastqs_merged/WT_4ku1_B26-8.fastq.gz
cat Fastqs/WT_non-act_B26-7_L00*/* > Fastqs_merged/WT_non-act_B26-7.fastq.gz
cat Fastqs/WT_non-act_B26-8_L00*/* > Fastqs_merged/WT_non-act_B26-8.fastq.gz
cat Fastqs/WT_non-act_B26-9_L00*/* > Fastqs_merged/WT_non-act_B26-9.fastq.gz

# QC analysis of merged data
mkdir FastQC
fastqc -o FastQC Fastqs_merged/*

# generating multiQC from above
mkdir MultiQC
cd MultiQC
multiqc ../FastQC/