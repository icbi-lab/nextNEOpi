![nextNEOpi overview](img/nextNEOpi_noBG.png)
# NeoEpitope predictions Nextflow Pipeline
Pipeline takes fastq files from Tumor and Normal samples (WES or WGS) and optionally RNAseq from tumor
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
* HLA-HD
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
HLA class I and class II alleles are predicted with OptiType and HLA-HD.
Class I and Class II neoepitopes are predicted with pVACseq using netMHCpan,
netMHCIIpan and mhcflurry. In addition mixMHC2pred is used as complement Class II
neoepitope predictor. Fusion neoantiges are predicted with NeoFuse.
CSiN immunogenicity score is reported for Class I, Class II and combined neoepitopes.
A GBM model [1] is be used to predict immunogenicity scores for MHC class I single nucleotide
variant (SNV) neoantigens 8-11 amino acid residues in length. Finally mixcr is run to predict the
TCR and BCR repertoire.

[1] https://github.com/vincentlaboratories/neoag/.

![nextNEOpi overview](img/nextNEOpi_small.png)



## 1. Installation

### 1.1 Nextflow

The command below may be used to install Nextflow. Please see also the installation instructions at:
https://www.nextflow.io/index.html#GetStarted


```
curl -s https://get.nextflow.io | bash

```

### 1.2 Analysis tools and software packages

The pipeline will install almost all required tools via conda environments or Singularity images.

The software that needs to be present on the system is **Java** (minimum version 8), **Nextflow** (see above), **Conda**,
**Singularity**.

**Optional but recommended:**
Due to license restrictions you may also need to download and install **HLA-HD** by your own, and set the installation path in ```conf/params.config```. _If HLA-HD is not available Class II neoepitopes will NOT be predicted_

_**[Manual installaton: Not recommended]:**_

If you prefer local installation of the analysis tools please install the following software:

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
* VEP 			(Version v103)
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
* HLA-HD
* ALLELECOUNT
* RSCRIPT (R > 3.6.1)
* SEQUENZA (3.0)
* CNVkit

all these tools need be available via the $PATH environment variable. However, you still need Java, Nextflow, Conda and Singularity installed on your system.

_**[End manual installation: not recommended]**_


### 1.2 References
The pipeline requires different reference files, indexes and databases:

please see ```conf/resources.config```

We prepared a bundle with all needed references, indexes and databases which can be obtained from:

https://apps-01.i-med.ac.at/resources/nextneopi/nextNEOpi_resources.tar.gz

download and extract the contents of the archive into the directory you specified for ```resourcesBaseDir``` in the ```conf/params.config``` file.

The structure should look as shown blow:

```
├── {resourcesBaseDir}
    ├── databases
    ├── ExomeCaptureKits
    └── references
```


**Notes**
1. You may also provide your own versions of these files. To do so, please change the ```conf/resources.config``` accordingly.
2. Due to license restriction, we do not provide a copy of the optional COSMIC database. If you also want to include COSMIC data, you may get a copy at https://cancer.sanger.ac.uk/cosmic
3. We provide the region and bait files for two different Exome capturing kits from Agilent:
   - SureSelect Human All Exon V6 exome
   - SureSelect Human All Exon V7 exome
   - Twist Human comprehensive exome

You may add your own region and bait files by defining an entry in ```conf/resources.config```



Refs:
- <https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle>
- <https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/>
- <ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/>
- <https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files>
- <https://www.gencodegenes.org/human/>


## 2. Usage
Before running the pipeline, the config files in the ```conf/``` directory may need to be edited. In the
```params.config``` parameters default settings are defined. The ```process.config```
is a template for the configuration of the single processes, you may check
the number of CPUs assigned for each process and adjust according to your systems capabilities.

Most pipeline parameters can be edited in the ```params.config``` file or changed on run time with command line options by using ```--NameOfTheParameter``` given in the ```params.config```.
References, databases should be edited in the ```resources.config``` file.

```
nextflow run nextNEOpi.nf --readsTumor <tumorFastq> --readsNormal <normalFastq> [--reads_RNAseq] | --batchFile <batchFile_FASTQ.csv | batchFile_BAM.csv> -profile singularity|conda,[cluster] [-resume] -config conf/params.config
```

**Profiles:** conda or singularity

We highly recommend to use either the ```singularity``` or ```conda``` profile. You can specify one of the two profiles using the option ```-profile singularity``` or ```-profile conda```. This way you do not have to care about installing all the required software including all
its dependencies.

**Profiles:** cluster

We strongly recommend to run the pipeline on a HPC cluster. You can enable runs in cluster mode by using a profile named e.g. **cluster** and the option ```-profile singularity,cluster``` or ```-profile conda,cluster```

For an example SGE cluster profile, please see ```profiles``` in ```conf/profiles.config```. You may uncomment and adjust the cluster profile to your scheduling system.


**Sequencing data input:**

Besides raw reads in **FASTQ** fromated files, input data may also be provided in **BAM** format.

**RNA reads from tag seq library i.e. 3-prime end sequencing protocol**

```--RNA_tag_seq``` turns off the "--trna-vaf" and "--trna-cov" filter from pVACseq epitope filtering. It also turns of HLA typing from RNAseq data. 3-prime end sequencing does not cover the entire transcript.

**Mandatory arguments:**

```--batchFile``` _[recommended]_

Make sure that your batchFile CSV includes the column names as shown in the examples below as header line. See also `example_batchFile_FASTQ.csv` or `example_batchFile_BAM.csv`

**FASTQ raw reads**
* e.g.: CSV-file with Tumor/Normal WES/WGS, and RNAseq reads, all paired end reads:

 | tumorSampleName | readsTumorFWD | readsTumorREV | normalSampleName | readsNormalFWD | readsNormalREV | readsRNAseqFWD | readsRNAseqREV | HLAfile | sex |
 | --------------- | ------------- | ------------- | ---------------- | -------------- | -------------- | -------------- | -------------- | ------- | ------ |
 | sample1 | Tumor1_reads_1.fastq | Tumor1_reads_2.fastq | normal1 | Normal1_reads_1.fastq | Normal1_reads_2.fastq | Tumor1_RNAseq_reads_1.fastq | Tumor1_RNAseq_reads_2.fastq | None | XX
 sample2 | Tumor2_reads_1.fastq | Tumor2_reads_2.fastq | normal2 | Normal2_reads_1.fastq | Normal2_reads_2.fastq | Tumor2_RNAseq_reads_1.fastq | Tumor2_RNAseq_reads_2.fastq | None | XY
 |... |
 sampleN | TumorN_reads_1.fastq | TumorN_reads_2.fastq | normalN | NormalN_reads_1.fastq | NormalN_reads_2.fastq | TumorN_RNAseq_reads_1.fastq | TumorN_RNAseq_reads_2.fastq | custom_HLAs.txt | XX


* e.g.:CSV-file with Tumor/Normal WES/WGS, and RNAseq reads, e.g. all single end reads:

 | tumorSampleName | readsTumorFWD | readsTumorREV | normalSampleName | readsNormalFWD | readsNormalREV | readsRNAseqFWD | readsRNAseqREV | HLAfile | sex |
 | --------------- | ------------- | ------------- | ---------------- | -------------- | -------------- | -------------- | -------------- | ------- | ------ |
 | sample1 | Tumor1_reads_1.fastq | None | normal1 | Normal1_reads_1.fastq | None | Tumor1_RNAseq_reads_1.fastq | None | None | XX
 sample2 | Tumor2_reads_1.fastq | None | normal2 | Normal2_reads_1.fastq | None | Tumor2_RNAseq_reads_1.fastq | None | None | XY
 |... |
 sampleN | TumorN_reads_1.fastq | None | normalN | NormalN_reads_1.fastq | None | TumorN_RNAseq_reads_1.fastq | None | custom_HLAs.txt | XX



* e.g.:CSV-file with Tumor/Normal WES/WGS, NO RNAseq reads, e.g. all single end reads:

 | tumorSampleName | readsTumorFWD | readsTumorREV | normalSampleName | readsNormalFWD | readsNormalREV | readsRNAseqFWD | readsRNAseqREV | HLAfile | sex |
 | --------------- | ------------- | ------------- | ---------------- | -------------- | -------------- | -------------- | -------------- | ------- | ------ |
 | sample1 | Tumor1_reads_1.fastq | None | normal1 | Normal1_reads_1.fastq | None | None | None | None | XX
 sample2 | Tumor2_reads_1.fastq | None | normal2 | Normal2_reads_1.fastq | None | None | None | None | XY
 |... |
 sampleN | TumorN_reads_1.fastq | None | normalN | NormalN_reads_1.fastq | None | None | None | custom_HLAs.txt | XX



or


```--readsTumor``` reads_{1,2}.fastq or reads_1.fastq; paired-end or single-end reads; FASTQ (may be gziped)

```--readsNormal``` reads_{1,2}.fastq or reads_1.fastq; paired-end or single-end reads; FASTA files (may be gziped)

(_optional but recommended_:
```--readsRNAseq``` reads_{1,2}.fastq or reads_1.fastq; paired-end or single-end reads; FASTQ (may be gziped))


**BAM files**

**Note:** If BAM files are used it is very much recommended that they also include also the unmapped and multimapping reads. These reads can be
helpful for HLA-typing.

* e.g.: CSV-file with Tumor/Normal WES/WGS, and RNAseq data:

 | tumorSampleName | bamTumor | normalSampleName | bamNormal | bamRNAseq | HLAfile | sex |
 | --------------- | ---------| ---------------- | --------- | --------- | ------- | --- |
 | sample1 | Tumor1.bam | normal1 | Normal1.bam | Tumor1_RNAseq.bam | None | XX
 | sample2 | Tumor2.bam | normal2 | Normal2.bam | Tumor2_RNAseq.bam | None | XY
 |... |
 | sampleN | TumorN.bam | normalN | NormalN.bam | TumorN_RNAseq.bam | None | XX


* e.g.:CSV-file with Tumor/Normal WES/WGS, NO RNAseq data:

 | tumorSampleName | bamTumor | normalSampleName | bamNormal | bamRNAseq | HLAfile | sex |
 | --------------- | ---------| ---------------- | --------- | --------- | ------- | --- |
 | sample1 | Tumor1.bam | normal1 | Normal1.bam | None | None | XX
 | sample2 | Tumor2.bam | normal2 | Normal2.bam | None | None | XY
 |... |
 | sampleN | TumorN.bam | normalN | NormalN.bam | None | None | XX



or


```--tumorBam``` tumor_1.bam

```--normalBam``` normal_1.bam

(_optional but recommended_:
```--rnaBam``` tumor_1_RNAseq.bam


**Notes**
- _You must not mix samples with single-end and paired-end reads in a batch file. Though, it is possible to have for e.g. all
DNA reads paired-end and all RNAseq reads single-end or vice-versa._

- in the ```HLAfile``` coulumn a user suppiled HLA types file may be specified for a given sample, see also ```--customHLA``` option below

- the ```sex``` column can be "XX", "female" or "Female", "XY", "male" or "Male". If not specified or "None" Male is assumed

- when providing paired end fastq files via commandline options (```--readsTumor, --readsNormal, --readsRNA```), please make sure you put the filename pattern into qoutes: e.g. ```--readsTumor "reads_{1,2}.fastq.gz"```


**Example run command with batchfile:**
```
nextflow run nextNEOpi.nf \
    --batchFile batchfile.csv \
    -config conf/params.config \
    --outputDir /data/results/nextNEOpi/myResults \
    --trim_adapters true \
    --trim_adapters_RNAseq true \
    --use_NetChop false \
    -profile singularity,cluster \
    -resume
```

**Optional argument:**

```--tumorSampleName```       tumor sample name. If not specified samples will be named according to the fastq filenames.

```--normalSampleName```      normal sample name. If not specified samples will be named according to the fastq filenames.

```--trim_adapters```         If true adpter sequences are automatically determined and will be trimmed from reads. If
                            ```--adapterSeq``` (string of atapter sequence) or ```--adapterSeqFile``` (fasta file with adapter sequences) is provided then adapters will be used as specified (no automatic detection).
                            Default: false

```--trim_adapters_RNAseq```  If true adpter sequences are automatically determined and will be trimmed from RNAseq reads. If
                            ```--adapterSeqRNAseq``` (string of atapter sequence) or ```--adapterSeqFileRNAseq``` (fasta file with adapter
                            sequences) is provided then adapters will be used as specified (no automatic detection).
                            Default: false

```--adapterSeq```            String of atapter sequence (see ```--trim_adapers```)
```--adapterSeqFile```        Fasta file with atapter sequence(s) (see ```--trim_adapers```)

```--adapterSeqRNAseq```      String of atapter sequence (see ```--trim_adapers_RNAseq```)
```--adapterSeqFileRNAseq```  Fasta file with atapter sequence(s) (see ```--trim_adapers_RNAseq```)

```--mutect2ponFile```        Panel of Normals file for Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-)
                            Default: false

```--priorityCaller```        Set the variant caller used as base for the hc variants. Only variants that are confirmed by any of the two confirming
                            callers (e..g. mutect1, varscan) will be retained. M2 = mutect2, M1 = mutect1, VS = varscan, ST = strelka
                            Default: M2

```--minAD```                 Minimum allelic depth (reads covering a variant)
                            Default: 5

```--use_NetChop```            Use NetChop to generate peptides
                            Default: false

```--TCR```                   Run mixcr for TCR prediction
                            Default: true

```--customHLA```             Provide a custom HLA types file. The HLA types in this file will be used in addition to those derived from the sequencing data in the WES/WGS/RNAseq fastq files. One type per line in 4 digit format (e.g. HLA-A*01:01:01)

```--HLAHD_DIR``` Specify the path to your HLA-HD installation. Needed if Class II neoantigens should be predicted.

```--HLA_force_RNA``` Use only RNAseq for HLA typing. Default: false

```--HLA_force_DNA``` Use only WES/WGS for HLA typing. Default: false

```--run_HLAHD_RNA``` Run HLA-HD also on RNAseq. Highly accurate but can be very slow
                      on larger fastq files. Default: false

```--disable_OptiType``` Disable OptiType for HLA typing. If set, **HLA-HD** or a user supplied **custom HLA file** must be available (see ```--HLAHD_DIR``` and/or ```--customHLA```)

```--sex```                Provide the sex of the sample (XX or Female, XY or Male, None)

```--pVACseq_filter_set```   Can be one of [standard, relaxed, custom]. The ```standard``` filter set is using the pVACseq default filters. The ```relaxed``` filter set is filtering only for ic50 < 500 & rank < 2 & expn-val > 2. With filter set ```custom``` users can define a custom set of filters by providing the desired filters (space separated) using the ```--pVACseq_custom_filters``` option. E.g. ```--pVACseq_filter_set custom --pVACseq_custom_filters "--binding-threshold 250 --percentile-threshold 1"```. For filter options please see also the pVACseq manual. Default: standard

```--pVACseq_custom_filters``` See ```--pVACseq_filter_set```

**Further options:**        There are many more options that can be set in the params.conf file or specified on the commandline
                            (see ```conf/params.config```)

## 3. Output
The Pipeline stores its ouput in the following structure:
```
RESULTS
├── analyses
│   ├── Subject_01
│   │   ├── 01_preprocessing
│   │   ├── 02_alignments
│   │   ├── 03_baserecalibration
│   │   ├── 03_realignment
│   │   ├── 04_expression
│   │   ├── 04_variations
│   │   ├── 05_vep
│   │   ├── 06_proteinseq
│   │   ├── 07_MutationalBurden
│   │   ├── 08_CNVs
│   │   ├── 09_CCF
│   │   ├── 10_HLA_typing
│   │   ├── 11_Fusions
│   │   ├── 12_pVACseq
│   │   ├── 13_mixMHC2pred
│   │   ├── 14_CSiN
│   │   ├── 14_IGS
│   │   ├── 15_BCR_TCR
│   │   └── QC
│   ├── Subject_02
│   │   ├── [...]
│   ├── [...]
│   │   ├── [...]
│   ├── Subject_n
│   │   ├── [...]
├── Documentation
├── neoantigens
│   ├── Subject_ID
│   │   ├── Class_I
│   │   ├── Class_II
│   │   └── Final_HLAcalls
│   ├── Subject_02
│   │   ├── [...]
│   ├── [...]
│   │   ├── [...]
│   ├── Subject_n
│   │   ├── [...]
├── pipeline_info
│   └── icbi
└── supplemental
    ├── 00_prepare_Intervals
    └── 01_prepare_CNVkit
```
