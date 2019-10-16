#!/usr/bin/env nextflow

/*
________________________________________________________________________________

                           C O N F I G U R A T I O N
________________________________________________________________________________
*/


if (params.readsTumor != "NO_FILE") {
  tumor_ch = Channel.fromFilePairs(params.readsTumor, flat:true)
} else  exit 1, "No tumor sample defined"

control_ch = Channel.fromFilePairs(params.readsControl, flat:true)

scatter_count = Channel.from(params.scatter_count)

reference = defineReference()
database = defineDatabases()

BWA           = file(params.BWA)
VARSCAN       = file(params.VARSCAN)
GATK4         = file(params.GATK4)
GATK3         = file(params.GATK3)
MUTECT1       = file(params.MUTECT1)
SAMTOOLS      = file(params.SAMTOOLS)
VEP           = file(params.VEP)
PICARD        = file(params.PICARD)
BAMREADCOUNT  = file(params.BAMREADCOUNT)
JAVA8         = file(params.JAVA8)

/*
________________________________________________________________________________

                              P R O C E S S E S
________________________________________________________________________________
*/

/*
*********************************************
**       P R E P R O C E S S I N G         **
*********************************************
*/
process 'SplitIntervals' {

  tag "SplitIntervals"
  publishDir "$params.outputDir/SplitIntervals/", mode: params.publishDirMode

  input:
  set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict,
    reference.IntervalsList])
  val x from scatter_count

  output:
  file "${IntervalName}/*-scattered.interval_list" into interval_ch

  script:
  IntervalName = IntervalsList.baseName
  """
  $GATK4 SplitIntervals \
  -R ${RefFasta}  \
  -scatter ${x} \
  -L ${IntervalsList} \
  -O ${IntervalName}
  """
}

process 'BwaTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/preprocessing/",
   mode: params.publishDirMode

  input:
  set TumorReplicateId, file(readsFWD), file(readsREV) from tumor_ch
  set file(RefFasta), file(RefIdx), file(RefDict), file(BwaRef) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])

  output:
  set TumorReplicateId, file("${TumorReplicateId}_aligned_sort.bam"),
   file("${TumorReplicateId}_aligned_sort.bai") into BwaSortTumor

  script:
  """
  $BWA mem \
  -R "@RG\\tID:${TumorReplicateId}\\tLB:${TumorReplicateId}\\tSM:${TumorReplicateId}\\tPL:ILLUMINA" \
  -M ${RefFasta} \
  -t ${task.cpus} \
  ${readsFWD} \
  ${readsREV} |java -jar $PICARD SortSam \
  INPUT=/dev/stdin \
  OUTPUT=${TumorReplicateId}_aligned_sort.bam \
  SORT_ORDER=coordinate \
  VALIDATION_STRINGENCY=LENIENT \
  CREATE_INDEX=true
  """
}

process 'MarkDuplicatesTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/preprocessing/",
   mode: params.publishDirMode

  input:
  set TumorReplicateId, file(bam), file(bai) from BwaSortTumor

  output:
  set TumorReplicateId, file("${TumorReplicateId}_aligned_sort_mkdp.bam"),
   file("${TumorReplicateId}_aligned_sort_mkdp.bai"),
   file("${TumorReplicateId}_aligned_sort_mkdp.txt") into MarkDuplicatesTumor1,
    MarkDuplicatesTumor2

  script:
  """
  java -jar $PICARD MarkDuplicates \
  INPUT=${bam} \
  OUTPUT=${TumorReplicateId}_aligned_sort_mkdp.bam \
  METRICS_FILE=${TumorReplicateId}_aligned_sort_mkdp.txt \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT
  """
}

if (params.readsControl != "NO_FILE") {
  process 'BwaControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/preprocessing/",
     mode: params.publishDirMode

    input:
    set ControlReplicateId, file(readsFWD), file(readsREV) from control_ch
    set file(RefFasta), file(RefIdx), file(RefDict), file(BwaRef) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])

    output:
    set ControlReplicateId, file("${ControlReplicateId}_aligned_sort.bam"),
     file("${ControlReplicateId}_aligned_sort.bai") into BwaSortControl

    script:
    """
    $BWA mem \
    -R "@RG\\tID:${ControlReplicateId}\\tLB:${ControlReplicateId}\\tSM:${ControlReplicateId}\\tPL:ILLUMINA" \
    -M ${RefFasta} \
    -t ${task.cpus} \
    ${readsFWD} \
    ${readsREV} |java -jar $PICARD SortSam \
    INPUT=/dev/stdin \
    OUTPUT=${ControlReplicateId}_aligned_sort.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true
    """
  }

  process 'MarkDuplicatesControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/preprocessing/",
     mode: params.publishDirMode

    input:
    set ControlReplicateId, file(bam), file(bai) from BwaSortControl

    output:
    set ControlReplicateId, file("${ControlReplicateId}_aligned_sort_mkdp.bam"),
     file("${ControlReplicateId}_aligned_sort_mkdp.bai"),
     file("${ControlReplicateId}_aligned_sort_mkdp.txt") into MarkDuplicatesControl1,
      MarkDuplicates2

    script:
    """
    java -jar $PICARD MarkDuplicates \
    INPUT=${bam} \
    OUTPUT=${ControlReplicateId}_aligned_sort_mkdp.bam \
    METRICS_FILE=${ControlReplicateId}_aligned_sort_mkdp.txt \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT
    """
  }
}

/*
*********************************************
**             M U T E C T 2               **
*********************************************
*/

process 'BaseRecalApplyTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/mutect2/processing/",
   mode: params.publishDirMode

  input:
  set TumorReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesTumor1
  set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
  set file(MilesGold), file(MilesGoldIdx), file(DBSNP), file(DBSNPIdx),
    file(KnownIndels), file(KnownIndelsIdx) from Channel.value(
        [database.MilesGold, database.MilesGoldIdx,
        database.DBSNP, database.DBSNPIdx,
        database.KnownIndels, database.KnownIndelsIdx])

  output:
  set TumorReplicateId, file("${TumorReplicateId}_bqsr.table") into BaseRecalibratorTumor
  set TumorReplicateId, file("${TumorReplicateId}_recal4.bam"),
   file("${TumorReplicateId}_recal4.bai") into ApplyTumor1, ApplyTumor2,
    ApplyTumor3, ApplyTumor4


  script:
  """
  $GATK4 BaseRecalibrator \
   -I ${bam} \
   -R ${RefFasta} \
   -L ${IntervalsList} \
   -O ${TumorReplicateId}_bqsr.table \
   --known-sites ${DBSNP} \
   --known-sites ${KnownIndels} \
   --known-sites ${MilesGold} && \
   $GATK4 ApplyBQSR \
   -I ${bam} \
   -R ${RefFasta} \
   -L ${IntervalsList} \
   -O ${TumorReplicateId}_recal4.bam \
   --bqsr-recal-file ${TumorReplicateId}_bqsr.table
  """
}

process 'GetPileupTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/mutect2/processing/",
   mode: params.publishDirMode

  input:
  set file(GnomAD), file(GnomADIdx) from Channel.value(
    [database.GnomAD, database.GnomADIdx])
  file(IntervalsList) from Channel.value([reference.IntervalsList])
  set TumorReplicateId, file(bam), file(bai) from ApplyTumor1

  output:
  set TumorReplicateId, file("${TumorReplicateId}_pileup.table") into PileupTumor

  script:
  """
  $GATK4 GetPileupSummaries \
  -I ${bam} -O ${TumorReplicateId}_pileup.table \
  -L ${IntervalsList} --variant ${GnomAD}
  """
}

process 'AnalyzeCovariates' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/QC/", mode: params.publishDirMode

  input:
  set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
  set file(DBSNP), file(DBSNPIdx), file(KnownIndels), file(KnownIndelsIdx),
    file(MilesGold), file(MilesGoldIdx) from Channel.value(
    [database.DBSNP, database.DBSNPIdx, database.KnownIndels,
    database.KnownIndelsIdx, database.MilesGold, database.MilesGoldIdx])
  set TumorReplicateId, file(recalTable) from BaseRecalibratorTumor
  set TumorReplicateId, file(bam), file(bai) from ApplyTumor2

  output:
  set TumorReplicateId, file("${TumorReplicateId}_postbqsr.table") into AnalyzeCovariates

  script:
  """
  $GATK4 BaseRecalibrator \
  -I ${bam} -R ${RefFasta} \
  -L ${IntervalsList} -O ${TumorReplicateId}_postbqsr.table \
  --known-sites ${DBSNP} \
  --known-sites ${KnownIndels} \
  --known-sites ${MilesGold} && \
  $GATK4 AnalyzeCovariates \
  -before ${recalTable} \
  -after ${TumorReplicateId}_postbqsr.table \
  -csv ${TumorReplicateId}_BQSR.csv \
  -plots ${TumorReplicateId}_BQSR.pdf
  """
}

process 'CollectSequencingArtifactMetrics' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/QC/", mode: params.publishDirMode

  input:
  set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict])
  set TumorReplicateId, file(bam), file(bai) from ApplyTumor3

  output:
  file("${TumorReplicateId}.pre_adapter_detail_metrics") into CollectSequencingArtifactMetrics

  script:
  """
  java -jar $PICARD CollectSequencingArtifactMetrics \
  I=${bam} R=${RefFasta} O=${TumorReplicateId}
  """
}

process 'Mutect2Tumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/mutect2/processing/",
   mode: params.publishDirMode

  input:
  set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict])
  set TumorReplicateId, file(Tumorbam), file(Tumorbai),
    file(IntervalsList) from ApplyTumor4.combine(interval_ch.flatten())

  output:
  file("${IntervalsList}.vcf.gz") into Mutect2Vcf
  file("${IntervalsList}.vcf.gz.stats") into Mutect2Stats
  set TumorReplicateId, file("${IntervalsList}.vcf.gz.tbi") into Mutect2Idx

  script:
  """
  $GATK4 Mutect2 \
  -R ${RefFasta} \
  -I ${Tumorbam} -tumor ${TumorReplicateId} \
  -L ${IntervalsList} --native-pair-hmm-threads ${task.cpus} \
  -O ${IntervalsList}.vcf.gz
  """
}

process 'MergeTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/mutect2/", mode: params.publishDirMode

  input:
  file(vcf) from Mutect2Vcf.collect()
  file (stats) from Mutect2Stats.collect()
  set TumorReplicateId, file(idx) from Mutect2Idx.groupTuple()

  output:
  set TumorReplicateId, file("${TumorReplicateId}_mutect2_raw.vcf.gz"),
  file("${TumorReplicateId}_mutect2_raw.vcf.gz.tbi"),
  file("${TumorReplicateId}_mutect2_raw.vcf.gz.stats") into mutect2

  script:
  """
  $GATK4 MergeVcfs \
  -I ${vcf.join(" -I ")} \
  -O ${TumorReplicateId}_mutect2_raw.vcf.gz

  $GATK4 MergeMutectStats \
  --stats ${stats.join(" --stats ")} \
  -O ${TumorReplicateId}_mutect2_raw.vcf.gz.stats
  """
}

process 'FilterMutec2Tumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/mutect2/", mode: params.publishDirMode

  input:
  set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict])
  set TumorReplicateId, file(pileup) from PileupTumor
  set TumorReplicateId, file(vcf), file(vcfIdx), file(vcfStats) from mutect2
  file(preAdapterDetail) from CollectSequencingArtifactMetrics

  output:
  set file("${TumorReplicateId}_mutect2_final.vcf.gz"),
    file("${TumorReplicateId}_mutect2_final.vcf.gz.tbi") into FilterMutect2

  script:
  """
  $GATK4 CalculateContamination \
  -I ${pileup} -O ${TumorReplicateId}_cont.table && \
  $GATK4 FilterMutectCalls \
  -R ${RefFasta} -V ${vcf} \
  --contamination-table ${TumorReplicateId}_cont.table \
  -O ${TumorReplicateId}_oncefiltered.vcf.gz && \
  $GATK4 FilterByOrientationBias \
  -V ${TumorReplicateId}_oncefiltered.vcf.gz \
  -P ${preAdapterDetail} \
  -O ${TumorReplicateId}_twicefitlered.vcf.gz && \
  $GATK4 SelectVariants \
  --variant ${TumorReplicateId}_twicefitlered.vcf.gz \
  -R ${RefFasta} --exclude-filtered true \
  --output ${TumorReplicateId}_mutect2_pass.vcf && \
  $GATK4 VariantFiltration \
  --variant ${TumorReplicateId}_mutect2_pass.vcf \
  -R ${RefFasta} \
  --genotype-filter-expression 'g.getAD().1 < 2' \
  --genotype-filter-name "AD.1_2" \
  --output ${TumorReplicateId}_mutect2_final.vcf.gz
  """
}

if (params.readsControl != 'NO_FILE') {
  process 'BaseRecalApplyControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/mutect2/processing/",
     mode: params.publishDirMode

    input:
    set ControlReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesControl1
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(MilesGold), file(MilesGoldIdx), file(DBSNP), file(DBSNPIdx),
     file(KnownIndels), file(KnownIndelsIdx) from Channel.value(
          [database.MilesGold, database.MilesGoldIdx,
          database.DBSNP, database.DBSNPIdx,
          database.KnownIndels, database.KnownIndelsIdx])

    output:
    set ControlReplicateId,
     file("${ControlReplicateId}_bqsr.table") into BaseRecalibratorControl
    set ControlReplicateId, file("${ControlReplicateId}_recal4.bam"),
     file("${ControlReplicateId}_recal4.bai") into ApplyControl

    script:
    """
    $GATK4 BaseRecalibrator \
     -I ${bam} \
     -R ${RefFasta} \
     -L ${IntervalsList} \
     -O ${ControlReplicateId}_bqsr.table \
     --known-sites ${DBSNP} \
     --known-sites ${KnownIndels} \
     --known-sites ${MilesGold} && \
     $GATK4 ApplyBQSR \
     -I ${bam} \
     -R ${RefFasta} \
     -L ${IntervalsList} \
     -O ${ControlReplicateId}_recal4.bam \
     --bqsr-recal-file ${ControlReplicateId}_bqsr.table
    """
  }

  process 'GetPileupControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/mutect2/processing/",
     mode: params.publishDirMode

    input:
    set file(GnomAD), file(GnomADIdx) from Channel.value(
      [database.GnomAD, database.GnomADIdx])
    file(intervals) from Channel.value([reference.IntervalsList])
    file ControlReplicateId, file(bam), file(bai) from ApplyControl

    output:
    set ControlReplicateId, file("${ControlReplicateId}_pileup.table") into PileupControl

    script:
    """
    $GATK4 GetPileupSummaries \
    -I ${bam} -O ${ControlReplicateId}_pileup.table \
    -L ${intervals} --variant ${GnomAD}
    """
  }
  process 'Mutect2' {
    tag "$TumorreplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/mutect2/processing/",
     mode: params.publishDirMode

    input:
    set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict])
    set TumorReplicateId, file(Tumorbam), file(Tumorbai),
      file(intervals) from ApplyTumor4.combine(interval_ch.flatten())
    set ControlReplicateId, file(Controlbam), file(Controlbai) from ApplyControl

    output:
    file("${intervals}.vcf.gz") into Mutect2Vcf
    file("${intervals}.vcf.gz.stats") into Mutect2Stats
    set TumorReplicateId, ControlReplicateId,
     file("${intervals}.vcf.gz.tbi") into Mutect2Idx

    script:
    """
    $GATK4 Mutect2 \
    -R ${RefFasta} \
    -I ${Tumorbam} -tumor ${TumorReplicateId} \
    -I ${Controlbam} -normal ${ControlReplicateId} \
    -L ${intervals} --native-pair-hmm-threads ${task.cpus} \
    -O ${intervals}.vcf.gz
    """
  }

  process 'Merge' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/mutect2/",
     mode: params.publishDirMode

    input:
    file(vcf) from Mutect2Vcf.collect()
    file (stats) from Mutect2Stats.collect()
    set TumorReplicateId, ControlReplicateId,
     file(idx) from Mutect2Idx.groupTuple()

    output:
    set TumorReplicateId, ControlReplicateId,
     file("${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz"),
  	 file("${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz.tbi"),
  	 file("${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz.stats") into
      mutect2

    script:
    """
    $GATK4 MergeVcfs \
    -I ${vcf.join(" -I ")} \
    -O ${TumorReplicateId}_mutect2_raw.vcf.gz

    $GATK4 MergeMutectStats \
    --stats ${stats.join(" --stats ")} \
    -O ${TumorReplicateId}_mutect2_raw.vcf.gz.stats
    """
  }

  process 'FilterMutec2' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/mutect2/",
     mode: params.publishDirMode

    input:
    set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference,RefDict])
    set TumorReplicateId, file(pileupTumor) from PileupTumor
    set ControlReplicateId, file(pileupControl) from PileupControl
    set TumorReplicateId, ControlReplicateId, file(vcf), file(vcfIdx),
     file(vcfStats) from mutect2
    set TumorReplicateId, file(preAdapterDetail) from CollectSequencingArtifactMetrics

    output:
    set TumorReplicateId, ControlReplicateId,
      file("${TumorReplicateId}_${ControlReplicateId}_mutect2_final.vcf.gz"),
      file("${TumorReplicateId}_${ControlReplicateId}_mutect2_final.vcf.gz.tbi") into
       FilterMutect2

    script:
    """
    $GATK4 CalculateContamination \
    -I ${pileupTumor} --matched-normal ${pileupControl} \
    -O ${TumorReplicateId}_${ControlReplicateId}_cont.table && \
    $GATK4 FilterMutectCalls \
    -R ${RefFasta} -V ${vcf} \
    --contamination-table ${TumorReplicateId}_${ControlReplicateId}_cont.table \
    -O ${TumorReplicateId}_${ControlReplicateId}_oncefiltered.vcf.gz && \
    $GATK4 FilterByOrientationBias \
    -V ${TumorReplicateId}_${ControlReplicateId}_oncefiltered.vcf.gz \
    -P ${preAdapterDetail} \
    -O ${TumorReplicateId}_${ControlReplicateId}_twicefitlered.vcf.gz && \
    $GATK4 SelectVariants \
    --variant ${TumorReplicateId}_${ControlReplicateId}_twicefitlered.vcf.gz \
    -R ${RefFasta} --exclude-filtered true \
    --output ${TumorReplicateId}_${ControlReplicateId}_mutect2_pass.vcf && \
    $GATK4 VariantFiltration \
    --variant ${TumorReplicateId}_${ControlReplicateId}_mutect2_pass.vcf \
    -R ${RefFasta} \
    --genotype-filter-expression 'g.getAD().1 < 2' \
    --genotype-filter-name "AD.1_2" \
    --output ${TumorReplicateId}_${ControlReplicateId}_mutect2_final.vcf.gz
    """
  }
}

/*
*********************************************
**             V A R S C A N               **
*********************************************
*/

process 'IndelRealignerTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/varscan/processing/",
   mode: params.publishDirMode

  input:
  set TumorReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesTumor2
  set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
  set file(KnownIndels), file(KnownIndelsIdx), file(MilesGold),
   file(MilesGoldIdx) from Channel.value([database.KnownIndels, database.KnownIndelsIdx,
      database.MilesGold, database.MilesGoldIdx])

  output:
  set TumorReplicateId, file("${TumorReplicateId}_mkdp_realign.bam"),
   file("${TumorReplicateId}_mkdp_realign.bai") into IndelRealignerTumor

  script:
  """
  $JAVA8 -jar $GATK3 \
  -T RealignerTargetCreator \
  --known ${MilesGold} \
  --known ${KnownIndels} \
  -R ${RefFasta} \
  -L ${IntervalsList} \
  -I ${bam} \
  -o target.list \
  -nt ${task.cpus} && \
  $JAVA8 -jar $GATK3 \
  -T IndelRealigner \
  -R ${RefFasta} \
  -L ${IntervalsList} \
  -I ${bam} \
  -targetIntervals target.list \
  -known ${KnownIndels} \
  -known ${MilesGold} \
  -nWayOut _realign.bam && \
  rm target.list
  """
}

process 'FixMateBaseRecalTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/varscan/processing/",
   mode: params.publishDirMode

  input:
  set TumorReplicateId, file(bam), file(bai) from IndelRealignerTumor
  set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
  set file(DBSNP), file(DBSNPIdx), file(KnownIndels), file(KnownIndelsIdx),
    file(MilesGold), file(MilesGoldIdx) from Channel.value(
    [database.DBSNP, database.DBSNPIdx, database.KnownIndels,
    database.KnownIndelsIdx, database.MilesGold, database.MilesGoldIdx])

  output:
  set TumorReplicateId, file("${TumorReplicateId}_fixmate.bam"),
   file("${TumorReplicateId}_fixmate_baserecal3.bai"),
   file("${TumorReplicateId}_bqsr3.table") into FixMateBaseRecalTumor

  script:
  """
  java -XX:ParallelGCThreads=8 -jar $PICARD FixMateInformation \
    I=${bam} \
    O=${TumorReplicateId}_fixmate.bam \
    CREATE_INDEX=true && \
    $JAVA8 -jar $GATK3\
    -T BaseRecalibrator \
    -R ${RefFasta} \
    -L ${IntervalsList} \
    -I ${TumorReplicateId}_fixmate.bam \
    -knownSites ${DBSNP} \
    -knownSites ${KnownIndels} \
    -knownSites ${MilesGold} \
    -nct ${task.cpus} \
    -o ${TumorReplicateId}_bqsr3.table
  """
}

process 'PrintReadsTumor' {
  tag "$TumorReplicateId"
  publishDir "$params.outputDir/$TumorReplicateId/varscan/processing",
   mode: params.publishDirMode

  input:
  set TumorReplicateId, file(bam), file(bai), file(bqsr3) from FixMateBaseRecalTumor
  set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.RefIdx,
     reference.IntervalsList, reference.IntervalsList])

  output:
  set TumorReplicateId, file("${TumorReplicateId}_printreads.bam"),
   file("${TumorReplicateId}_printreads.bai") into PrintReadsTumor

  script:
  """
  $JAVA8 -jar $GATK3  \
  -T PrintReads \
  -R ${RefFasta} \
  -L ${IntervalsList} \
  -I ${bam} \
  -BQSR ${bqsr3} \
  -nct ${task.cpus} \
  -o ${TumorReplicateId}_printreads.bam
  """
}

if (params.readsControl != "NO_FILE") {

  process 'IndelRealignerControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/varscan/processing/",
     mode: params.publishDirMode

    input:
    set ControlReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesControl2
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(KnownIndels), file(KnownIndelsIdx), file(MilesGold),
      file(MilesGoldIdx) from Channel.value(
      [database.KnownIndels, database.KnownIndelsIdx,
      database.MilesGold, database.MilesGoldIdx])

    output:
    set ControlReplicateId, file("${ControlReplicateId}_mkdp_realign.bam"),
     file("${ControlReplicateId}_mkdp_realign.bai") into IndelRealignerControl

    script:
    """
    $JAVA8 -jar $GATK3 \
    -T RealignerTargetCreator \
    --known ${MilesGold} \
    --known ${KnownIndels} \
    -R ${RefFasta} \
    -L ${IntervalsList} \
    -I ${bam} \
    -o target.list \
    -nt ${task.cpus} && \
    $JAVA8 -jar $GATK3 \
    -T IndelRealigner \
    -R ${RefFasta} \
    -L ${IntervalsList} \
    -I ${bam} \
    -targetIntervals target.list \
    -known ${KnownIndels} \
    -known ${MilesGold} \
    -nWayOut ${ControlReplicateId}_mkdp_realign.bam && \
    rm target.list
    """
  }

  process 'FixMateBaseRecalControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/varscan/processing/",
     mode: params.publishDirMode

    input:
    set ControlReplicateId, file(bam), file(bai) from IndelRealignerControl
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(DBSNP), file(DBSNPIdx), file(KnownIndels), file(KnownIndelsIdx),
      file(MilesGold), file(MilesGoldIdx) from Channel.value(
      [database.DBSNP, database.DBSNPIdx, database.KnownIndels, database.KnownIndelsIdx,
      database.MilesGold, database.MilesGoldIdx])

    output:
    set ControlReplicateId, file("${ControlReplicateId}_fixmate.bam"),
      file("${ControlReplicateId}_fixmate.bai"),
       file("${ControlReplicateId}_bqsr3.table") into FixMateBaseRecalControl

    script:
    """
    java -XX:ParallelGCThreads=8 -jar $PICARD FixMateInformation \
      I=${bam} \
      O=${ControlReplicateId}_fixmate.bam \
      CREATE_INDEX=true && \
      $JAVA8 -jar $GATK3\
      -T BaseRecalibrator \
      -R ${RefFasta} \
      -L ${IntervalsList} \
      -I ${ControlReplicateId}_fixmate.bam \
      -knownSites ${DBSNP} \
      -knownSites ${KnownIndels} \
      -knownSites ${MilesGold} \
      -nct ${task.cpus} \
      -o ${ControlReplicateId}_bqsr3.table
    """
  }

  process 'PrintReadsControl' {
    tag "$ControlReplicateId"
    publishDir "$params.outputDir/$ControlReplicateId/varscan/processing/",
     mode: params.publishDirMode

    input:
    set ControlReplicateId, file(bam), file(bai), file(bqsr3) from FixMateBaseRecalControl
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])

    output:
    set ControlReplicateId, file("${ControlReplicateId}_printreads.bam"),
     file("${ControlReplicateId}_printreads.bai") into PrintReadsControl

    script:
    """
    $JAVA8 -jar $GATK3  \
    -T PrintReads \
    -R ${RefFasta} \
    -L ${IntervalsList} \
    -I ${bam} \
    -BQSR ${bqsr3} \
    -nct ${task.cpus} \
    -o ${ControlReplicateId}_printreads.bam
    """
  }

  process 'ReorderSam' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/varscan/processing/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bamTumor), file(baiTumor) from PrintReadsTumor
    set ControlReplicateId, file(bamControl), file(baiControl) from PrintReadsControl
    set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict])

    output:
    set TumorReplicateId, file("${TumorReplicateId}_reorder_tumor.bam"),
     file("${TumorReplicateId}_reorder_tumor.bai") into ReorderSamTumor
    set ControlReplicateId, file("${ControlReplicateId}_reorder_control.bam"),
     file("${ControlReplicateId}_reorder_control.bai") into ReorderSamControl

    script:
    """
    java -jar $PICARD ReorderSam \
    REFERENCE=${RefFasta} \
    INPUT=${bamTumor} \
    CREATE_INDEX=true \
    OUTPUT=${TumorReplicateId}_reorder_tumor.bam && \
    java -jar $PICARD ReorderSam \
    REFERENCE=${RefFasta} \
    INPUT=${bamControl} \
    CREATE_INDEX=true \
    OUTPUT=${ControlReplicateId}_reorder_tumor.bam
    """
  }

  process 'VarscanSomatic' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/varscan/", mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bamTumor), file(baiTumor) from ReorderSamTumor
    set ControlReplicateId, file(bamControl), file(baiControl) from ReorderSamControl
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsBed) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsBed])

    output:
    set TumorReplicateId, ControlReplicateId,
     file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.vcf"),
     file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.vcf") into VarscanSomatic

    script:
    """
    mkfifo ${TumorReplicateId}_${ControlReplicateId}_mpileup.fifo
    $SAMTOOLS mpileup \
    -q ${params.q} \
    -f ${RefFasta} \
    -l ${IntervalsBed} \
    ${bamControl} ${bamTumor} > ${TumorReplicateId}_${ControlReplicateId}_mpileup.fifo &
    java -jar $VARSCAN somatic \
    ${TumorReplicateId}_${ControlReplicateId}_mpileup.fifo \
    ${TumorReplicateId}_${ControlReplicateId}_varscan \
    --output-vcf 1 \
    --mpileup 1 \
    --min-coverage ${params.min_cov} \
    --min-coverage-normal ${params.min_cov_normal} \
    --min-coverage-tumor ${params.min_cov_tumor} \
    --min-freq-for-hom ${params.min_freq_for_hom} \
    --tumor-purity ${params.tumor_purity} \
    --p-value ${params.Somatic_pvalue} \
    --somatic-p-value ${params.Somatic_somaticpvalue} \
    --strand-filter ${params.strand_filter} && \
    rm -f ${TumorReplicateId}_${ControlReplicateId}_mpileup.fifo
    """
  }

  process 'ProcessSomatic' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/varscan/processing",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, ControlReplicateId, file(snp), file(indel) from VarscanSomatic

    output:
    set TumorReplicateId, ControlReplicateId,
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.hc.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.LOH.hc.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscna.snp.LOH.hc.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Germline.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Germline.hc.vcf") into ProcessSomaticSNP
    set TumorReplicateId, ControlReplicateId,
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Somatic.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Somatic.hc.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.LOH.hc.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscna.indel.LOH.hc.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Germline.vcf"),
      file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Germline.hc.vcf") into ProcessSomaticIndel

    script:
    """
    java -jar $VARSCAN processSomatic \
    ${snp} \
    --min-tumor-freq ${params.min_tumor_freq} \
    --max-normal-freq ${params.max_normal_freq} \
    --p-value ${paramsProcessSomatic_pvalue} && \
    java -jar $VARSCAN processSomatic \
    ${indel} \
    --min-tumor-freq ${params.min_tumor_freq} \
    --max-normal-freq ${params.max_normal_freq} \
    --p-value ${paramsProcessSomatic_pvalue}
    """
  }

  process 'FilterVarscan' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/varscan/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam), file(bai) from ReorderSamTumor
    set TumorReplicateId, ControlReplicateId, file(snpSomatic), file(snpSomaticHc),
      file(snpLOH), file(snpLOHhc), file(snpGerm), file(snpGemHc) from ProcessSomaticSNP
    set TumorReplicateId, ControlReplicateId, file(indelSomatic), file(indelSomaticHc),
      file(indelLOH), file(indelLOHhc), file(indelGerm),
      file(indelGemHc) from ProcessSomaticIndel
    set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict])

    output:
    set TumorReplicateId, ControlReplicateId,
     file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.hc.filtered.vcf")
    set TumorReplicateId, ControlReplicateId,
     file("${TumorReplicateId}_${ControlReplicateId_varscan.indel.Somatic.hc.filtered.vcf}")

    script:
    """
    cat ${snpSomaticHc} | awk '{if (!/^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' |$BAMREADCOUNT \
    -q${params.min_map_cov} -b${params.min_base_q} \
    -w1 -l /dev/stdin \
    -f ${RefFasta} ${bam} |java -jar $VARSCAN fpfilter \
    ${snpSomaticHc} \
    /dev/stdin \
    --output-file ${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.hc.filtered.vcf && \
    cat ${indelSomaticHc} | awk '{if (! /^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' |$BAMREADCOUNT \
    -q${params.min_map_cov} -b${params.min_base_q} \
    -w1 -l /dev/stdin \
    -f ${RefFasta} ${bam} |java -jar $VARSCAN fpfilter \
    ${indelSomaticHc} \
    /dev/stdin \
    --output-file ${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Somatic.hc.filtered.vcf
    """
  }
}

/*
*********************************************
**             M U T E C T 1               **
*********************************************
*/

if (params.readsControl != "NO_FILE") {
  process 'MutectTumorNormal' {
    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/mutect1/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bamTumor), file(baiTumor) from PrintReadsTumor
    set ControlReplicateId, file(bamControl), file(baiControl) from PrintReadsControl
    set file(RefFasta), file(RefIdx), file(RefDict) file(IntervalsList) from Channel.value(
      [refernce.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(DBSNP), file(DBSNPIdx), file(Cosmic), file(CosmicIdx) from Channel.value(
      [database.DBSNP, database.DBSNPIdx, database.Cosmic, database.CosmicIdx])

    output:
    set TumorReplicateId, ControlReplicateId,
      file("${TumorReplicateId}_${ControlReplicateId}_mutect1_raw.vcf.gz"),
      file("${TumorReplicateId}_${ControlReplicateId}_mutect1_raw.vcf.gz.idx")
      file("${TumorReplicateId}_${ControlReplicateId}_mutect1_raw.stats.txt") into Mutect1raw
    set TumorReplicateId, ControlReplicateId,
      file("${TumorReplicateId}_${ControlReplicateId}_mutect1.final.vcf.gz"),
      file("${TumorReplicateId}_${ControlReplicateId}_mutect1.final.vcf.gz.idx")
      file("${TumorReplicateId}_${ControlReplicateId}_mutect1_final.stats.txt") into Mutect1filter

    script:
    """
    $JAVA8 -jar $MUTECT1 \
    --analysis_type MuTect \
    --reference_sequence ${RefFasta} \
    --cosmic ${Cosmic} \
    --dbsnp ${DBSNP} \
    -L ${IntervalsList} \
    --input_file:normal ${bamControl} \
    --input_file:tumor ${bamTumor} \
    --out ${TumorReplicateId}_${ControlReplicateId}_mutect1.stats.txt \
    --vcf ${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.vcf.gz && \
    $GATK4 SelectVariants \
    --variant ${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.vcf.gz \
    -R ${RefFasta} \
    --exclude-filtered true \
    --output ${TumorReplicateId}_${ControlReplicateId}_mutect1_pass.vcf && \
    $GATK4 VariantFiltration \
    --variant ${TumorReplicateId}_${ControlReplicateId}_mutect1_pass.vcf \
    -R ${RefFasta} \
    --genotype-filter-expression "g.getAD().1 < 2" \
    --genotype-filter-name "AD.1_2" \
    --output ${TumorReplicateId}_${ControlReplicateId}_mutect1_final.vcf.gz
    """
  }
}else {
  process 'MutectTumor' {
    tag "$TumorReplicateId"
    publishDir "$params.publishDir/$TumorReplicateId/mutect1/", mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bamTumor), file(baiTumor) from PrintReadsTumor
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
      [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(DBSNP), file(DBSNPIdx), file(Cosmic), file(CosmicIdx) from Channel.value(
      [database.DBSNP, database.DBSNPIdx, database.Cosmic, database.CosmicIdx])

    output:
    set TumorReplicateId,
      file("${TumorReplicateId}mutect1_raw.vcf.gz"),
      file("${TumorReplicateId}_mutect1_raw.vcf.gz.idx")
      file("${TumorReplicateId}_mutect1_raw.stats.txt") into Mutect1raw
    set TumorReplicateId,
      file("${TumorReplicateId}_mutect1.final.vcf.gz"),
      file("${TumorReplicateId}_mutect1.final.vcf.gz.idx")
      file("${TumorReplicateId}_mutect1_final.stats.txt") into Mutect1filter

    script:
    """
    $JAVA8 -jar $MUTECT1 \
    --analysis_type MuTect \
    --reference_sequence ${RefFasta} \
    --cosmic ${Cosmic} \
    --dbsnp ${DBSNP} \
    -L ${IntervalsList} \
    --input_file:tumor ${bamTumor} \
    --out ${TumorReplicateId}_mutect1.stats.txt \
    --vcf ${TumorReplicateId}_mutect1.raw.vcf.gz && \
    $GATK4 SelectVariants \
    --variant ${TumorReplicateId}_mutect1.raw.vcf.gz \
    -R ${RefFasta} \
    --exclude-filtered true \
    --output ${TumorReplicateId}_mutect1_pass.vcf && \
    $GATK4 VariantFiltration \
    --variant ${TumorReplicateId}_mutect1_pass.vcf \
    -R ${RefFasta} \
    --genotype-filter-expression "g.getAD().1 < 2" \
    --genotype-filter-name "AD.1_2" \
    --output ${TumorReplicateId}_mutect1_final.vcf.gz
    """
  }
}

/*
________________________________________________________________________________

                            F U N C T I O N S
________________________________________________________________________________

*/

def checkParamReturnFileReferences(item) {
  params."${item}" = params.references."${item}"
  return file(params."${item}")
}

def checkParamReturnFileDatabases(item) {
  params."${item}" = params.databases."${item}"
  return file(params."${item}")
}

def defineReference() {
  if (params.references.size() != 7) exit 1, """
  ERROR: Not all References needed found in configuration
  Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
  """
  return [
    'RefFasta'   : checkParamReturnFileReferences("RefFasta"),
    'RefIdx'       : checkParamReturnFileReferences("RefIdx"),
    'RefDict'    : checkParamReturnFileReferences("RefDict"),
    'BwaRef'     : checkParamReturnFileReferences("BwaRef"),
    'IntervalsList': checkParamReturnFileReferences("IntervalsList"),
    'VepFasta'     : checkParamReturnFileReferences("VepFasta"),
    'IntervalsBed' : checkParamReturnFileReferences("IntervalsBed")
  ]
}

def defineDatabases() {
  if (params.databases.size() !=10) exit 1, """
  ERROR: Not all Databases needed found in configuration
  Please check if Mills_and_1000G_gold_standard, CosmicCodingMuts, DBSNP, GnomAD, and knownIndels are given.
  """
  return [
  'MilesGold'   : checkParamReturnFileDatabases("MilesGold"),
  'MilesGoldIdx'   : checkParamReturnFileDatabases("MilesGoldIdx"),
  'Cosmic'      : checkParamReturnFileDatabases("Cosmic"),
  'CosmicIdx'      : checkParamReturnFileDatabases("CosmicIdx"),
  'DBSNP'       : checkParamReturnFileDatabases("DBSNP"),
  'DBSNPIdx'       : checkParamReturnFileDatabases("DBSNPIdx"),
  'GnomAD'      : checkParamReturnFileDatabases("GnomAD"),
  'GnomADIdx'      : checkParamReturnFileDatabases("GnomADIdx"),
  'KnownIndels' : checkParamReturnFileDatabases("KnownIndels"),
  'KnownIndelsIdx' : checkParamReturnFileDatabases("KnownIndelsIdx")
  ]
}
