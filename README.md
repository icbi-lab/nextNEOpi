# Whole-Exome Sequencing Nextflow Pipeline
Pipeline takes fastq file from tumor (and control) sampels for indels and SNPs 
predictions. 

Therefore 3 different callers are used, depending on input:
* MuTect2
* MuTect1
* VarScan2 (only when matched control sample is given)

It outputs a txt file with the annotated and filtered SNPs and Indels, which
where called with a minimum of 2 of the callers.

![Beschreibung](img/flowchart.png)

## 1. Installation

## 1.1 Singularity Image

To run the pipeline using singularity, the image has to be build:
```
singularity build coolNameHere.sif docker://icbi/wes
```
## 1.2 Software
To run the pipeline local, the required software has to be installed:
* JAVA7 			 (Version 1.7)
* JAVA8 			 (Version 1.8)
* BWA 			 (Version 0.7.17)
* SAMTOOLS 		 (Version 1.9)
* PICARD 			 (Version 2.20.0)
* GATK3 			 (Version 3.8-0)
* GATK4 			 (Version 4.1.3.0)
* VARSCAN 		 (Version 2.4.3)
* MUTECT1 		 (Version 1.1.7)
* BAMREADCOUNT 		 (Version 0.8.0)
* VEP 			 (Version 2.0)


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
* Cosmic/Idx - CosmicCodingMuts.vcf - Cosmic conding mutations, VCF file, IDX file
* DBSNP/Idx - Homo_sapiens_assembly.dbsnp.vcf/idx - SNPS, microsatellites, and small-scale insertions and deletions, VCF file, IDX file
* GnomAD/Idx - small_exac_common_3.vcf/idx - exonix sites only for contamination estimation from GATK, VCF file, IDX file
* KnownIdenls/Idx - Homo_sapiens_assembly.known_indels.vcf/idx - Known Indels from GATK resource Bundle, VCF file, IDX file

## 2. Usage
Before running the pipeline, the config files has to be edited. In the
params.config parameters like references, databases and samples are defined. The sge.config 
is a template for the configuration to run the pipeline on cluster.

Every parameter can be edited in the params file or with the command lind by using --NameOfTheParameter given in the params.config.
References, Databases and Software should be edited in the params.config.

```
nextflow run wes.nf "--readsTumor|--batchFile" "[--readsControl]" "--IntervalsList" "--IntervalsBed" [--single_end]
```
#### Singularity
The singularity mode has to be anabled in the params.config file and the path to the image has to be edited.

#### Single-end reads:
**--single_end:** sets parameter to TRUE (default false)

#### Mandatory arguments:
**--readsTumor:** 		 reads_{1,2}.fastq or reads_1.fastq; 		 paired-end or single-end reads; FASTA file (can be zipped)

or

**--batchFile:**
* CSV-file, paired-end T/N reads:

 sampleId,readsTumorFWD,readsTumorREV,readsControlFWD,readsControlREV,group
 sample1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,Control1_reads_1.fastq,Control1_reads_2.fastq,group1
 sample2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,Control2_reads_1.fastq,Control2_reads_2.fastq,group1
 ...
 sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,ControlN_reads_1.fastq,ControlN_reads_2.fastq,groupX

* CSV-file, single-end T only reads:

 sampleId,readsTumorFWD,readsTumorREV,readsControlFWD,readsControlREV,group
 sample1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,NO_FILE,NO_FILE,group1
 sample2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,NO_FILE,NO_FILE,group1
 ...
 sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,NO_FILE,,NO_FILE,groupX

**--IntervalsList:** 	 intervals.list; 		 interval.list file for targeting

**--IntervalsBed:** 		 intervals.bed; 			 interval.bed file for targeting

#### Optional argument:
**--readsControl:** 		 reads_{1,2}.fastq or reads_1.fastq; 		 paired-end or single-end reads; FASTA file (can be zipped)
**--sampleName**          sample name. If not specified samples will be named according to the fastq filenames.  

**--sampleName**          sample name. If not specified samples will be named according to the fastq filenames.  

## 3. Output
The Pipeline creates an ouput directory with the following structure:
```
RESULTS
├──TumorSample
│   ├── 1_preprocessing
│   ├── 2_QC
│   ├── 3_mutect2
│   │   ├── 1_processing
│   ├── 4_mutect1
│   ├── 5_varscan
│   │   ├── 1_processing
│   └── 6_vep
├──ControlSample
│   ├── 1_preprocessing
│   ├── 3_mutect2
│   │   ├── 1_processing
│   └── 6_vep
│       ├── 1_processing
└──SplitIntervals

```
