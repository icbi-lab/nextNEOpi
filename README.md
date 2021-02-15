# NeoEpitope predictions Nextflow Pipeline
Pipeline takes fastq file from Tumor and Normal samples and optionally RNAseq from tumor
to predict neoepitopes.

The pipeline uses the following tools:
* MuTect2
* MuTect1
* VarScan2
* Strekla2/Manta
* Sequenza
* ASCAT
* CNVkit
* OptiType
* HLAHD
* pVACseq (netMHCpan, netMHCIIpan, mhcflurry)
* NeoFuse
* mixMHC2pred
* mixcr

It outputs a vcf files with the annotated and filtered SNPs and Indels, which
where called with each of the callers and a high confidence vcf file (hc) in
which only variants that were called by a minimum of 2 of the callers are listed.
All vcf files are annotatd with VEP. In addition the germline variants are called
using HaploTypeCaller and a phased vcf for pVACseq is generated as well.
Copy number variations are analyzed using CNVkit, ASCAT, and sequenza. Tumor purity
is estimated by ASCAT and Sequenza and is used to derive the clonality measure for
the predicted neoantigens. Tumor mutational burden (TMB) is calculated for all
variants over the entire read covered genome and for coding variants on read covered
exons.
HLA class I and class II alleles are predicted with OptiType and HLAHD.
Class I and Class II neoepitopes are predicted with pVACseq using netMHCpan,
netMHCIIpan and mhcflurry. In addition mixMHC2pred is used as complement Class II
neoepitope predictor. Fusion neoantiges are predicted with NeoFuse.
CSiN immunogenicity score is reported for Class I, Class II and combined neoepitopes.
A GBM model is be used to predict immunogenicity scores for MHC class I single nucleotide
variant (SNV) neoantigens 8-11 amino acid residues in length. https://github.com/vincentlaboratories/.
Finally mixcr is run to predict TCRs.


![Beschreibung](img/flowchart.png)

## 1. Installation

## 1.1 Nextflow

Please see the installation instructions at:
https://www.nextflow.io/index.html#GetStarted


```
curl -s https://get.nextflow.io | bash

```

## 1.2 Software
The pipeline will install almost all required tools via conda environments or Singularity images.
The only software that needs to be available is Java (minimum version 8), Nextflow (see above), conda,
Singularity.

Due to license concerns you also need to download and install HLA-HD by your on.

[Not recommended]:
If you can not run either conda or Singularity you need to install the required software tools
locally.

* FASTQC        (Version >= 0.11.8)
* FASTP         (Version >= v0.20.1)
* JAVA7 		(Version 1.7)
* JAVA8 		(Version 1.8)
* BWA 			(Version 0.7.17)
* SAMTOOLS 		(Version 1.9)
* GATK3 		(Version 3.8-0)
* GATK4 		(Version >= 4.1.7.0)
* VARSCAN 		(Version 2.4.3)
* MUTECT1 		(Version 1.1.7) ---- optional
* BAMREADCOUNT 	(Version 0.8.0)
* VEP 			(Version v102)
* BGZIP
* TABIX
* BCFTOOLS
* MANTA
* STRELKA
* SAMBAMBA
* OPTITYPE
* PYTHON
* PERL
* CONDA
* YARA
* HLAHD
* ALLELECOUNT
* RSCRIPT
* SEQUENZA
* CNVkit
[End not recommended]


## 1.2 References
the Pipeline requires different references and databases:

please see ```resources.config``

We prepared a bundle with all needed references, indexes and databases which can be obtained from:

https://apps-01.i-med.ac.at/resources/nextNEOpi/nextNEOpi_resources.tar.gz

download and extract the contents of the archive into the directory you specified for ```resourcesBaseDir```

The structure should look as shown blow:

```
├── {resourcesBaseDir}
    ├── databases
    ├── ExomeCaptureKits
    └── references
```


Ref:
<https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle>
<https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/>
<ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/>


## 2. Usage
Before running the pipeline, the config files may need to be edited. In the
params.config parameters default settings are defined. The process.config
is a template for the configuration of the single processes, you may check
the number of cpus assigned for each process.

Every parameter can be edited in the params file or with the command lind by using --NameOfTheParameter given in the params.config.
References, Databases and Software should be edited in the resources.config.

```
nextflow run nextNEOpi.nf "--readsTumor <tumorFastq> --readsNormal <nomralFastq>"|--batchFile <batchFile.csv>" [--single_end] -profile singularity|conda,[cluster] [-resume]
```

#### Profiles: conda or singularity
We highly recommend to use either the ```singularity``` or ```conda``` profile. You can specify one of the two profiles using the option ```-profile singularity``` or ```-profile conda```

#### Profiles: cluster
We recommend to run the pipeline on a HPC cluster. You can enable runs in cluster mode by the option ```-profile singularity,cluster``` or ```-profile conda,cluster```
Please see profiles. to adjust the cluster profile to your schedulin system.


#### Single-end reads:
**--single_end:** sets parameter to TRUE (default false)

#### Mandatory arguments:
**--readsTumor:** 		 reads_{1,2}.fastq or reads_1.fastq; 		 paired-end or single-end reads; FASTA
**--readsNormal:** 		 reads_{1,2}.fastq or reads_1.fastq; 		 paired-end or single-end reads; FASTA files (can be zipped)

or

**--batchFile:**
* CSV-file, paired-end T/N reads, paired-end RNAseq reads:

 tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,readsRNAseqFWD,readsRNAseqREV,HLAfile,gender,group
 sample1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,normal1,Normal1_reads_1.fastq,Normal1_reads_2.fastq,Tumor1_RNAseq_reads_1.fastq,Tumor1_RNAseq_reads_2.fastq,None,XX,group1
 sample2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,normal2,Normal2_reads_1.fastq,Normal2_reads_2.fastq,Tumor2_RNAseq_reads_1.fastq,Tumor2_RNAseq_reads_2.fastq,None,XY,group1
 ...
 sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,normalN,NormalN_reads_1.fastq,NormalN_reads_2.fastq,TumorN_RNAseq_reads_1.fastq,TumorN_RNAseq_reads_2.fastq,XX,groupX

* CSV-file, single-end T/N reads, single-end RNAseq reads:

 tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,readsRNAseqFWD,readsRNAseqREV,HLAfile,gender,group
 sample1,Tumor1_reads_1.fastq,None,normal1,Normal1_reads_1.fastq,None,Tumor1_RNAseq_reads_1.fastq,None,None,XX,group1
 sample2,Tumor2_reads_1.fastq,None,normal2,Normal2_reads_1.fastq,None,Tumor1_RNAseq_reads_1.fastq,None,None,XY,group1
 ...
 sampleN,TumorN_reads_1.fastq,None,normalN,NormalN_reads_1.fastq,None,Tumor1_RNAseq_reads_1.fastq,None,None,None,groupX

* CSV-file, single-end T/N reads, NO RNAseq reads:

 tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,readsRNAseqFWD,readsRNAseqREV,HLAfile,gender,group
 sample1,Tumor1_reads_1.fastq,None,normal1,Normal1_reads_1.fastq,None,None,None,None,XX,group1
 sample2,Tumor2_reads_1.fastq,None,normal2,Normal2_reads_1.fastq,None,None,None,None,XY,group1
 ...
 sampleN,TumorN_reads_1.fastq,None,normalN,NormalN_reads_1.fastq,None,None,None,None,XX,groupX


Note: You must not mix samples with single-end and paired-end reads in a batch file. Though, it is possible to have for e.g. all
DNA reads paired-end and all RNAseq reads single-end or vice-versa.

Note: in the HLAfile coulumn a user suppiled HLA types file may be specified for a given sample, see also --customHLA option below

Note: gender can be XX or Female, XY or Male. If not specified or "None" Male is assumed

#### Optional argument:
**--tumorSampleName**       tumor sample name. If not specified samples will be named according to the fastq filenames.

**--normalSampleName**      normal sample name. If not specified samples will be named according to the fastq filenames.

**--trim_adapters**         If true adpter sequences are automatically determined and will be trimmed from reads. If
                            --adapterSeq (string of atapter sequence) or --adapterSeqFile (fasta file with adapter sequences) is provided
                            then adapters will be used as specified (no automatic detection).
                            Default: false

**--trim_adapters_RNAseq**  If true adpter sequences are automatically determined and will be trimmed from RNAseq reads. If
                            --adapterSeqRNAseq (string of atapter sequence) or --adapterSeqFileRNAseq (fasta file with adapter
                            sequences) is provided then adapters will be used as specified (no automatic detection).
                            Default: false

**--adapterSeq**            String of atapter sequence (see --trim_adapers)
**--adapterSeqFile**        Fasta file with atapter sequence(s) (see --trim_adapers)

**--adapterSeqRNAseq**      String of atapter sequence (see --trim_adapers_RNAseq)
**--adapterSeqFileRNAseq**  Fasta file with atapter sequence(s) (see --trim_adapers_RNAseq)

**--mutect2ponFile**        Panel of Normals file for Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360040510131-CreateSomaticPanelOfNormals-BETA-)
                            Default: false

**--priorityCaller**        Set the variant caller used as base for the hc variants. Only variants that are confirmed by any of the two confirming
                            callers (e..g. mutect1, varscan) will be retained. m2 = mutect2, m1 = mutect1, vs = varscan, st = strelka
                            Default: m2

**--minAD**                 Minimum allelic depth (reads covering a variant)
                            Default: 5

**-use_NetChop**            Use NetChop to generate peptides
                            Default: false

**--TCR**                   Run mixcr for TCR prediction
                            Default: true
**--customHLA**             Provide a custom HLA types file. One type per line in 4 digit format (e.g. HLA-A*01:01:01)

**--gender**                Provide the gender of the sample (XX or Female, XY or Male)

**Further options:**        There are many more options that can be set in the params.conf file or specified on the commandline
                            (see params.conf)

## 3. Output  (this is not uptodate, will be changed when we have the final structure)
The Pipeline creates an ouput directory with the following structure:
```
RESULTS
├── 00_prepare_Intervals
│   └── SplitIntervals
├── CRC01
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
│   ├── 07_PhasedVCF
│   ├── 08_OptiType
│   ├── 09_HLA_HD
│   ├── 10_NeoFuse
│   ├── 11_pVACseq
│   ├── 12_mixMHC2pred
│   └── 13_TCRs
├── Documentation
│   ├── pipeline_report.html
│   └── pipeline_report.txt
└── pipeline_info
    └── icbi
```
