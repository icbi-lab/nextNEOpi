/*
 * Define default parameters
 */

 params.genome     = "data/genome.fa"
 params.variants   = "data/known_variants.vcf.gz"
 params.blacklist  = "data/blacklist.bed"
 params.reads      = "data/reads/ENCSR000COQ1_{1,2}.fastq.gz"
 params.results    = "results"
 params.gatk       = "/opt/broad/GenomeAnalysisTK.jar" 

/*
 * Parse input parameters
 */

genome_file = file(params.genome)
varaints_file = file(params.variants)
intervals_file = file(params.intervals)
reads_ch = Channel.fromFilePairs(params.reads)
GATK = params.gatk

/*
 * Process 1A: Create a Fasta genome indes with samtools
 */

process '1A_prepare_genome_samtools' {

  input:
    file genome from genome_file

  output:
    file "${genome}.fai" into genome_index_ch

  script:
  """
  samtools faidx ${genome}
  """
}

/*
 * Process 1B: Create a Fasta genome sequqnce dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {

  input:
    file genome from genome_file

  output:
    file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  PICARD = `which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}

/*
 * Process 1C: Create the genome index file for STAR
 */

process '1C_prepare_star_genome_index' {

  input:
    file genome from genome_file

  output:
    file genome_dir into genome_dir_ch

  script:
  """
  mkdir genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${tasl.cpus}
"""
}
