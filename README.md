Pipeline for RNA-seq FASTQ and VCF files processing

The pipeline takes as input:

1. RNAseq fastq files (single- or paired-end)
2. VCF file *

*Make sure that the VCF file is VEP annotated

## 1. Usage

Before running the pipeline, the config files has to be edited. In the
params.config parameters like references, databases and samples are defined. The sge.config
is a template for the configuration to run the pipeline on cluster.
Every parameter can be edited in the params file or with the command lind by using --NameOfTheParameter given in the params.config.
References, Databases and Software should be edited in the params.config.

`nextflow run rna.nf --reads1 --reads2 --tsv_file --vcf_file --single_end --phased_vcf --phased_vcf_tbi`

### Singularity
The singularity mode has to be enabled in the params.config file and the path to the image has to be edited.

### Single-end reads:
--single_end: sets parameter to TRUE (default false)

## 2. Output

```
/output/dir/
├── NeoFuse
│   ├── Sample_name
│   │   ├── Arriba
│   │   ├── LOGS
│   │   ├── NeoFuse
│   │   ├── OptiType
│   │   └── TPM
│   └── STAR_idx
└── pVACseq
    └── Sample_name
        └── MHC_Class_I
```

