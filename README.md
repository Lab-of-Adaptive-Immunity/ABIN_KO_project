## About

This depository contians scripts to re-do analyses and figures for Janusova et al. paper: 

ABIN1 is a negative regulator of effector functions in cytotoxic T cells

Link: Coming SOON!

Data Accession: [GSE245397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE245397)

## Running these scripts

1. Clone this project. You'll end up with ABIN_KO_project directory.

2. Download file Mus_musculus.GRCm39.106.gtf from Ensembl in ABIN_KO_project directory. The file can be obtained [here](https://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz). 

3. Create Fastqs directory in ABIN_KO_project, go in it (cd ABIN_KO_project) and download all .fastq files from GSE245397 accession number to it.

4. In ABIN_KO_project, create indices for STAR 2.7.10a aligner form GRCm39 v106 mouse assembly (Ensembl). Files needed can be downloaded [here](https://ftp.ensembl.org/pub/release-106/) (navigate to gtf and fasta directories and download files mendtioned in command). The command to generate indices stored in directory GRCm39 is: STAR --runThreadN 20 --runMode genomeGenerate --genomeDir GRCm39 v106 --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm39.106.gtf


5. Return to ABIN_KO_project directory (cd ..) and successively run following files. This will merge fastq files from the same samples and then align them :

* exp01_merge_fastq_QC.sh
* exp01_run_STAR.sh
* exp02_merge_fastq_QC.sh
* exp02_run_STAR.sh

6. The above should create bamfiles (directories Bamfiles1 and Bamfiles2) needed for analyses. The analyses are contained in files exp1_RNAseq.Rmd, exp2_RNAseq.Rmd and exp1_exp2_combined_analysis.Rmd.

Now you should have files with complete analyses. Notably, tehre should be directories Figures, Rds_data, Tables and GSEA_tables that respectively contain figures, rds files with DESeq (result) objects, count tables and list of genes used for GSEA.

