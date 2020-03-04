# Whole-Exome Sequencing Nextflow Pipeline
Pipeline takes fastq file from Tumor and Normal samples for indels and SNPs 
predictions. 

The pipeline uses 5 different somatic variant callers:
* MuTect2
* MuTect1
* VarScan2
* Strekla2/Manta

It outputs a vcf files with the annotated and filtered SNPs and Indels, which
where called with each of the callers and a high confidence vcf file (hc) in
which only variants that were called by a minimum of 2 of the callers are listed.
All vcf files are annotatd with VEP. In addition the germline variants are called
using HaploTypeCaller and a phased vcf for pVACseq is generated as well.

![Beschreibung](img/flowchart.png)

## 1. Installation

## 1.1 Singularity Image

To run the pipeline using singularity, the image has to be build:
```
singularity build coolNameHere.sif docker://icbi/wes
```
## 1.2 Software
To run the pipeline local, the required software has to be installed:
* FASTQC        (Version >= 0.11.8)
* FLEXBAR        (Version 3.5)
* JAVA7 			 (Version 1.7)
* JAVA8 			 (Version 1.8)
* BWA 			 (Version 0.7.17)
* SAMTOOLS 		 (Version 1.9)
* PICARD 			 (Version 2.20.0)
* GATK3 			 (Version 3.8-0)
* GATK4 			 (Version 4.1.4.1)
* VARSCAN 		 (Version 2.4.3)
* MUTECT1 		 (Version 1.1.7)
* BAMREADCOUNT 		 (Version 0.8.0)
* VEP 			 (Version 2.0)
* BGZIP
* TABIX
* BCFTOOLS
* MANTA
* STRELKA


## 1.2 References
the Pipeline requires different references and databases:

References:
* RefFasta - reference.fa - Reference Genome; FASTA file
* RefIdx - reference.fai - Referene Genome Index, FAI file
* RefDict - reference.dict - Reference Genome Dictionary, DICT file
* BwaRef - reference.{amb,sa,ann,pac,bwt} - Reference Genome prepared for BWA mem
* VepFasta - reference.toplevel.fa - Reference genome; ENSEMBL toplevel

Databases:
* MillsGold/Idx - Mills_and_1000G_gold_standard.indels.vcf/idx -  Gold standard Indels database, VCF file, IDX File
* 1000G high confidence SNPs - VCF + IDX
* HapMap VCF + IDF 
* Cosmic/Idx - CosmicCodingMuts.vcf - Cosmic conding mutations, VCF file, IDX file
* DBSNP/Idx - Homo_sapiens_assembly.dbsnp.vcf/idx - SNPS, microsatellites, and small-scale insertions and deletions, VCF file, IDX file
* GnomAD/Idx - small_exac_common_3.vcf/idx - exonix sites only for contamination estimation from GATK, VCF file, IDX file
* GnomADfull VCF + IDX
* KnownIdenls/Idx - Homo_sapiens_assembly.known_indels.vcf/idx - Known Indels from GATK resource Bundle, VCF file, IDX file

see also:
<https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle>
<https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/>
<ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/>

## 2. Usage
Before running the pipeline, the config files has to be edited. In the
params.config parameters like references, databases and samples are defined. The sge.config 
is a template for the configuration to run the pipeline on cluster.

Every parameter can be edited in the params file or with the command lind by using --NameOfTheParameter given in the params.config.
References, Databases and Software should be edited in the params.config.

```
nextflow run wes.nf "--readsTumor <tumorFastq> --readsNormal <nomralFastq>"|--batchFile <batchFile.csv>" "--BaitsBed" "--RegionsBed" [--single_end]
```
#### Singularity
The singularity mode has to be anabled in the params.config file and the path to the image has to be edited.

#### Single-end reads:
**--single_end:** sets parameter to TRUE (default false)

#### Mandatory arguments:
**--readsTumor:** 		 reads_{1,2}.fastq or reads_1.fastq; 		 paired-end or single-end reads; FASTA 
**--readsNormal:** 		 reads_{1,2}.fastq or reads_1.fastq; 		 paired-end or single-end reads; FASTA files (can be zipped)

or

**--batchFile:**
* CSV-file, paired-end T/N reads:

 tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,group
 sample1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,normal1,Normal1_reads_1.fastq,Normal1_reads_2.fastq,group1
 sample2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,normal2,Normal2_reads_1.fastq,Normal2_reads_2.fastq,group1
 ...
 sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,normalN,NormalN_reads_1.fastq,NormalN_reads_2.fastq,groupX

* CSV-file, single-end T/N reads:

 tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,group
 sample1,Tumor1_reads_1.fastq,None,normal1,Normal1_reads_1.fastq,None,group1
 sample2,Tumor2_reads_1.fastq,None,normal2,Normal2_reads_1.fastq,None,group1
 ...
 sampleN,TumorN_reads_1.fastq,None,normalN,NormalN_reads_1.fastq,None,groupX

**--BaitsBed:** 	 baits.bed; 		 baits.bed file for Exon baits

**--RegionsBed:** 		 regions.bed; 			 regions.bed file for Exon targeting

#### Optional argument:
**--tumorSampleName**          tumor sample name. If not specified samples will be named according to the fastq filenames.  

**--normalSampleName**          normal sample name. If not specified samples will be named according to the fastq filenames.  

**--trim_adapters**          If true Illumina universal adpter (AGATCGGAAGAG) will be trimmed from reads unless
                             --adapterSeq (string of atapter sequence) or --adapterSeqFile (fasta file with adapter sequences) is provided.

**--adapterSeq**             String of atapter sequence (see --trim_adapers)
**--adapterSeqFile**         Fasta file with atapter sequence(s) (see --trim_adapers)

## 3. Output
The Pipeline creates an ouput directory with the following structure:
```
RESULTS
├── 00_prepare_Intervals
│   ├── S07604514_Covered.list
│   ├── S07604514_Regions.list
│   ├── S07604514_Regions_merged_padded.bed
│   ├── S07604514_Regions_merged_padded.bed.gz
│   ├── S07604514_Regions_merged_padded.bed.gz.tbi
│   ├── S07604514_Regions_merged_padded.list
│   └── SplitIntervals
├── CRC26
│   ├── 01_preprocessing
│   ├── 02_QC
│   ├── 03_manta_somatic
│   ├── 03_mutect1
│   ├── 03_mutect2
│   ├── 03_strelka_somatic
│   ├── 03_varscan
│   ├── 04_haplotypeCaller
│   ├── 05_hcVCF
│   ├── 06_vep
│   └── 07_PhasedVCF
├── Documentation
│   ├── pipeline_report.html
│   └── pipeline_report.txt
└── pipeline_info
    └── icbi
```
