#!/usr/bin/env nextflow

log.info ""
log.info " NEXTFLOW ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
log.info "-------------------------------------------------------------------------"
log.info "         W H O L E - E X O M E   P I P E L I N E             "
log.info "-------------------------------------------------------------------------"
log.info ""
log.info " Finding SNPs and Indels in tumor only or tumor + matched control samples"
log.info " using three different caller: \n \t * MuTect1 \n \t * MuTect2 \n \t * VarScan2 (only with matched control samples)"
log.info ""
log.info "-------------------------------------------------------------------------"
log.info "C O N F I G U R A T I O N"
log.info ""
log.info "Command Line: \t\t " + workflow.commandLine
log.info "Working Directory: \t " + params.workDir
log.info "Output Directory: \t " + params.outputDir
log.info ""
log.info "I N P U T"
log.info ""
if  (params.readsControl != "NO_FILE") log.info " Reads Tumor: \t\t " + params.readsTumor
if  (params.readsControl != "NO_FILE") log.info "Reads Control: \t\t " + params.readsControl
log.info ""
log.info "Please check --help for further instruction"
log.info "-------------------------------------------------------------------------"

/*
________________________________________________________________________________

                            C O N F I G U R A T I O N
________________________________________________________________________________
*/
if (params.help) exit 0, helpMessage()


// default is not to process a batchfile
params.batchFile = false

// set single_end variable to supplied param
single_end = params.single_end

// set control reads variable to supplied param
readsControl = params.readsControl

// did not get a CSV batch file: just run a single sample
if (! params.batchFile) {
        // create channel with sampleName/reads file set
        if (params.readsTumor != "NO_FILE") {
            // Sample name and reads file is passed via cmd line options
            // Sample name to use, if not given uses fastq file simpleName
            params.sampleName = "undefined"
            if(!params.single_end) {
                sampleName = params.sampleName != "undefined" ? params.sampleName : file(params.readsTumor)[0].simpleName
            } else {
                sampleName = params.sampleName != "undefined" ? params.sampleName : file(params.readsTumor).simpleName
            }

            Channel
                   .fromFilePairs(params.readsTumor)
                   .map { reads -> tuple(sampleName, reads[1][0], reads[1][1], "None") }
                   .set { tumor_ch }

            if (params.readsControl != "NO_FILE") {
                Channel
                       .fromFilePairs(params.readsControl)
                       .map { reads -> tuple(sampleName, reads[1][0], reads[1][1], "None") }
                       .set { control_ch }
            }
        } else  {
            exit 1, "No tumor sample defined"
        }

} else {
        // batchfile ()= csv with sampleId and T reads [N reads] [and group]) was provided
        // create channel with all sampleId/reads file sets from the batch file
        
        // check if reverse reads are specified, if not set up single end processing
        // check if Control reads are specified, if not set up no control processing
        // attention: if one of the samples is no-control or single end, all others
        // will be handled as such. You might want to process mixed sample data as
        // separate batches.
        batchCSV  = file(params.batchFile).splitCsv(header:true)
        for ( row in batchCSV ) {
            if (row.readsTumorREV == "None") {
                single_end = true
            }
            if (row.readsControlFWD != "None") {
                readsControl = "FILE"
            }                
        }

        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.sampleId, 
                                    file(row.readsTumorFWD),
                                    file(row.readsTumorREV),
                                    row.group) }
                .set { tumor_ch }


        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.sampleId, 
                                    file(row.readsControlFWD),
                                    file(row.readsControlREV),
                                    row.group) }
                .set { control_ch }

}


reference = defineReference()
database = defineDatabases()
mkTmpDir()

scatter_count = Channel.from(params.scatter_count)

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
JAVA7         = file(params.JAVA7)
PERL          = file(params.PERL)

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
// Splitting interval file in 20(default) files for scattering Mutect2

    tag "SplitIntervals"
    publishDir "$params.outputDir/SplitIntervals/", mode: params.publishDirMode

    input:
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict,
        reference.IntervalsList])
    val x from scatter_count

    output:
    file "${IntervalName}/*-scattered.interval_list" into interval_ch1, interval_ch2

    script:
    IntervalName = IntervalsList.baseName
    """
    mkdir -p ${params.tmpDir}

    $GATK4 SplitIntervals \
        --tmp-dir ${params.tmpDir} \
        -R ${RefFasta}  \
        -scatter ${x} \
        -L ${IntervalsList} \
        -O ${IntervalName}
    """
}

if (single_end) {
    process 'BwaTumorSingle' {
    // Aligning tumor reads to reference, sort and index; create BAMs

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
         mode: params.publishDirMode

        input:
        set TumorReplicateId, file(readsFWD) from tumor_ch
        set file(RefFasta), file(RefIdx), file(RefDict), file(BwaRef) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])

        output:
        set TumorReplicateId, file("${TumorReplicateId}_aligned.bam") into BwaTumor_ch

        script:
        """
        $BWA mem \
        -R "@RG\\tID:${TumorReplicateId}\\tLB:${TumorReplicateId}\\tSM:${TumorReplicateId}\\tPL:ILLUMINA" \
        -M ${RefFasta} \
        -t ${task.cpus} \
        ${readsFWD}  | \
        $SAMTOOLS view -Shb -o ${TumorReplicateId}_aligned.bam -
        """
    }
} else {
    process 'BwaTumor' {
    // Aligning tumor reads to reference, sort and index; create BAMs

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
         mode: params.publishDirMode

        input:
        set TumorReplicateId, file(readsFWD), file(readsREV) from tumor_ch
        set file(RefFasta), file(RefIdx), file(RefDict), file(BwaRef) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])

        output:
        set TumorReplicateId, file("${TumorReplicateId}_aligned.bam") into BwaTumor_ch

        script:
        """
        $BWA mem \
        -R "@RG\\tID:${TumorReplicateId}\\tLB:${TumorReplicateId}\\tSM:${TumorReplicateId}\\tPL:ILLUMINA" \
        -M ${RefFasta} \
        -t ${task.cpus} \
        ${readsFWD} \
        ${readsREV}  | \
        $SAMTOOLS view -Shb -o ${TumorReplicateId}_aligned.bam -
        """
    }
}

process 'MarkDuplicatesTumor' {
// Mark duplicates with Picard

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam) from BwaTumor_ch

    output:
    set TumorReplicateId, file("${TumorReplicateId}_aligned_sort_mkdp.bam"),
     file("${TumorReplicateId}_aligned_sort_mkdp.bam.bai") into MarkDuplicatesTumor0
    set TumorReplicateId, file("${TumorReplicateId}_aligned_sort_mkdp.bam"),
     file("${TumorReplicateId}_aligned_sort_mkdp.bam.bai"),
     file("${TumorReplicateId}_aligned_sort_mkdp.txt") into MarkDuplicatesTumor1,
        MarkDuplicatesTumor2

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 MarkDuplicatesSpark \
        --tmp-dir ${params.tmpDir} \
        -I ${bam} \
        -O ${TumorReplicateId}_aligned_sort_mkdp.bam \
        -M ${TumorReplicateId}_aligned_sort_mkdp.txt \
        --create-output-bam-index true \
        --read-validation-stringency LENIENT \
        --conf 'spark.executor.cores=${task.cpus}' 2> /dev/stdout
    """
}

process 'alignmentMetrics' {
// Generate HS metrics using picard

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/02_QC/",
        mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam), file(bai) from MarkDuplicatesTumor0
    set file(RefFasta), file(RefIdx), file(IntervalsList), file(BaitIntervalsList) from Channel.value(
    [reference.RefFasta, reference.RefIdx, reference.IntervalsList, reference.BaitIntervalsList])


    output:
    file("${TumorReplicateId}.HS.metrics.txt")
    file("${TumorReplicateId}.perTarget.coverage.txt")
    file("${TumorReplicateId}.AS.metrics.txt")
    file("${TumorReplicateId}.flagstat.txt")

    script:
    """
    mkdir -p ${params.tmpDir}
    java -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectHsMetrics \
        TMP_DIR=${params.tmpDir} \
        INPUT=${bam} \
        OUTPUT=${TumorReplicateId}.HS.metrics.txt \
        R=${RefFasta} \
        BAIT_INTERVALS=${BaitIntervalsList} \
        TARGET_INTERVALS=${IntervalsList} \
        PER_TARGET_COVERAGE=${TumorReplicateId}.perTarget.coverage.txt && \
    java -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectAlignmentSummaryMetrics \
        TMP_DIR=${params.tmpDir} \
        INPUT=${bam} \
        OUTPUT=${TumorReplicateId}.AS.metrics.txt \
        R=${RefFasta} &&
    $SAMTOOLS flagstat -@${task.cpus} ${bam} > ${TumorReplicateId}.flagstat.txt
    """
}


if (readsControl != "NO_FILE" && single_end) {
    process 'BwaControlSinlge' {
    // Aligning control reads to reference, sort and index; create BAMs

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/1_preprocessing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(readsFWD) from control_ch
        set file(RefFasta), file(RefIdx), file(RefDict), file(BwaRef) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])

        output:
        set ControlReplicateId, file("${ControlReplicateId}_Control_aligned.bam") into BwaControl_ch

        script:
        """
        $BWA mem \
        -R "@RG\\tID:${ControlReplicateId}\\tLB:${ControlReplicateId}\\tSM:${ControlReplicateId}\\tPL:ILLUMINA" \
        -M ${RefFasta} \
        -t ${task.cpus} \
        ${readsFWD} | \
        $SAMTOOLS  view -Shb -o ${ControlReplicateId}_Control_aligned.bam -
        """
    }
} else if (readsControl != "NO_FILE" && !single_end) {
    process 'BwaControl' {
    // Aligning control reads to reference, sort and index; create BAMs

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/01_preprocessing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(readsFWD), file(readsREV) from control_ch
        set file(RefFasta), file(RefIdx), file(RefDict), file(BwaRef) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.BwaRef])

        output:
        set ControlReplicateId, file("${ControlReplicateId}_Control_aligned.bam") into BwaControl_ch

        script:
        """
        $BWA mem \
        -R "@RG\\tID:${ControlReplicateId}\\tLB:${ControlReplicateId}\\tSM:${ControlReplicateId}\\tPL:ILLUMINA" \
        -M ${RefFasta} \
        -t ${task.cpus} \
        ${readsFWD} \
        ${readsREV} | \
        $SAMTOOLS view -Shb -o ${ControlReplicateId}_Control_aligned.bam -
        """
    }

    process 'MarkDuplicatesControl' {
    // Mark duplicates with Picard

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/01_preprocessing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(bam) from BwaControl_ch

        output:
        set ControlReplicateId, file("${ControlReplicateId}_Control_aligned_sort_mkdp.bam"),
         file("${ControlReplicateId}_Control_aligned_sort_mkdp.bam.bai") into MarkDuplicatesControl0
        set ControlReplicateId, file("${ControlReplicateId}_Control_aligned_sort_mkdp.bam"),
         file("${ControlReplicateId}_Control_aligned_sort_mkdp.bam.bai"),
         file("${ControlReplicateId}_Control_aligned_sort_mkdp.txt") into MarkDuplicatesControl1,
             MarkDuplicatesControl2

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 MarkDuplicatesSpark \
            --tmp-dir ${params.tmpDir} \
            -I ${bam} \
            -O ${ControlReplicateId}_Control_aligned_sort_mkdp.bam \
            -M ${ControlReplicateId}_Control_aligned_sort_mkdp.txt \
            --create-output-bam-index true \
            --read-validation-stringency LENIENT \
            --conf 'spark.executor.cores=${task.cpus}' 2> /dev/stdout
        """
    }

    process 'alignmentMetricsControl' {
    // Generate HS metrics using picard

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/02_QC/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(bam), file(bai) from MarkDuplicatesControl0
        set file(RefFasta), file(RefIdx), file(IntervalsList), file(BaitIntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.IntervalsList, reference.BaitIntervalsList])


        output:
        file("${ControlReplicateId}_Control.HS.metrics.txt")
        file("${ControlReplicateId}_Control.perTarget.coverage.txt")
        file("${ControlReplicateId}_Control.AS.metrics.txt")
        file("${ControlReplicateId}_Control.flagstat.txt")


        script:
        """
        mkdir -p ${params.tmpDir}
        java -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectHsMetrics \
            TMP_DIR=${params.tmpDir} \
            INPUT=${bam} \
            OUTPUT=${ControlReplicateId}_Control.HS.metrics.txt \
            R=${RefFasta} \
            BAIT_INTERVALS=${BaitIntervalsList} \
            TARGET_INTERVALS=${IntervalsList} \
            PER_TARGET_COVERAGE=${ControlReplicateId}_Control.perTarget.coverage.txt && \
        java -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectAlignmentSummaryMetrics \
            TMP_DIR=${params.tmpDir} \
            INPUT=${bam} \
            OUTPUT=${ControlReplicateId}_Control.AS.metrics.txt \
            R=${RefFasta} && \
        $SAMTOOLS flagstat -@${task.cpus} ${bam} > ${ControlReplicateId}_Control.flagstat.txt
        """
    }
}

/*
*********************************************
**             M U T E C T 2               **
*********************************************
*/

process 'BaseRecalApplyTumor' {
/*
 BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
 Recalibration (BQSR)
 ApplyBQSR (GATK4): apply BQSR table to reads
*/

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesTumor1
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(MillsGold), file(MillsGoldIdx), file(DBSNP), file(DBSNPIdx),
        file(KnownIndels), file(KnownIndelsIdx) from Channel.value(
                [database.MillsGold, database.MillsGoldIdx,
                database.DBSNP, database.DBSNPIdx,
                database.KnownIndels, database.KnownIndelsIdx])

    output:
    set TumorReplicateId, file("${TumorReplicateId}_bqsr.table") into BaseRecalibratorTumor
    set TumorReplicateId, file("${TumorReplicateId}_recal4.bam"),
     file("${TumorReplicateId}_recal4.bam.bai") into (ApplyTumor1, ApplyTumor2, ApplyTumor3,
            ApplyTumor4, ApplyTumor5)


    script:
    """
    mkdir -p ${params.tmpDir}

    java -jar $PICARD SetNmMdAndUqTags \
        TMP_DIR=${params.tmpDir} \
        R=${RefFasta} \
        I=${bam} \
        O=fixed.bam \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT && \
    $GATK4 BaseRecalibratorSpark \
        TMP_DIR=${params.tmpDir} \
        -I fixed.bam \
        -R ${RefFasta} \
        -L ${IntervalsList} \
        -O ${TumorReplicateId}_bqsr.table \
        --known-sites ${DBSNP} \
        --known-sites ${KnownIndels} \
        --known-sites ${MillsGold} \
        --conf 'spark.executor.cores=${task.cpus}' && \
     $GATK4 ApplyBQSRSpark \
        -I ${bam} \
        -R ${RefFasta} \
        -L ${IntervalsList} \
        -O ${TumorReplicateId}_recal4.bam \
        --bqsr-recal-file ${TumorReplicateId}_bqsr.table \
        --conf 'spark.executor.cores=${task.cpus}'
    """
}

process 'GetPileupTumor' {
// GetPileupSummaries (GATK4): tabulates pileup metrics for inferring contamination

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
     mode: params.publishDirMode

    input:
    set file(GnomAD), file(GnomADIdx) from Channel.value(
        [database.GnomAD, database.GnomADIdx])
    file(IntervalsList) from Channel.value([reference.IntervalsList])
    set TumorReplicateId, file(bam), file(bai) from ApplyTumor1

    output:
    set TumorReplicateId, file("${TumorReplicateId}_pileup.table") into PileupTumor1,
     PileupTumor2

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 GetPileupSummaries \
        --tmp-dir ${params.tmpDir} \
        -I ${bam} \
        -O ${TumorReplicateId}_pileup.table \
        -L ${IntervalsList} \
        --variant ${GnomAD}
    """
}

process 'AnalyzeCovariates' {
/*
 2nd BaseRecalibrator (GATK4)
 AnalyzeCovariates (GATK4): creates plots to visualize base recalibration results
*/

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/02_QC/", mode: params.publishDirMode

    input:
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(DBSNP), file(DBSNPIdx), file(KnownIndels), file(KnownIndelsIdx),
        file(MillsGold), file(MillsGoldIdx) from Channel.value(
        [database.DBSNP, database.DBSNPIdx, database.KnownIndels,
        database.KnownIndelsIdx, database.MillsGold, database.MillsGoldIdx])
    set TumorReplicateId, file(recalTable) from BaseRecalibratorTumor
    set TumorReplicateId, file(bam), file(bai) from ApplyTumor2

    output:
    set TumorReplicateId, file("${TumorReplicateId}_postbqsr.table") into AnalyzeCovariates

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 BaseRecalibratorSpark \
        --tmp-dir ${params.tmpDir} \
        -I ${bam} \
        -R ${RefFasta} \
        -L ${IntervalsList} \
        -O ${TumorReplicateId}_postbqsr.table \
        --known-sites ${DBSNP} \
        --known-sites ${KnownIndels} \
        --known-sites ${MillsGold} \
        --conf 'spark.executor.cores=${task.cpus}' && \
    $GATK4 AnalyzeCovariates \
        --tmp-dir ${params.tmpDir} \
        -before ${recalTable} \
        -after ${TumorReplicateId}_postbqsr.table \
        -csv ${TumorReplicateId}_BQSR.csv
    """
}

process 'CollectSequencingArtifactMetrics' {
/*
 CollectSequencingArtifactMetrics (Picard): collect metrics to quantify single-base
 sequencing artifacts
*/

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/02_QC/", mode: params.publishDirMode

    input:
    set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict])
    set TumorReplicateId, file(bam), file(bai) from ApplyTumor3

    output:
    file("${TumorReplicateId}.pre_adapter_detail_metrics") into (CollectSequencingArtifactMetrics1, CollectSequencingArtifactMetrics2)

    script:
    """
    mkdir -p ${params.tmpDir}
    java -jar $PICARD CollectSequencingArtifactMetrics \
        TMP_DIR=${params.tmpDir} \
        I=${bam} \
        R=${RefFasta} \
        O=${TumorReplicateId}
    """
}

if (readsControl == "NO_FILE") {
    process 'Mutect2Tumor' {
    // Call somatic SNPs and indels via local re-assembly of haplotypes; only tumor sample

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
         mode: params.publishDirMode

        input:
        set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict])
        set TumorReplicateId, file(Tumorbam), file(Tumorbai),
            file(IntervalsList) from ApplyTumor4.combine(interval_ch1.flatten())

        output:
        file("${IntervalsList}.vcf.gz") into Mutect2VcfTumor
        file("${IntervalsList}.vcf.gz.stats") into Mutect2StatsTumor
        set TumorReplicateId, file("${IntervalsList}.vcf.gz.tbi") into Mutect2IdxTumor

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 Mutect2 \
            --tmp-dir ${params.tmpDir} \
            -R ${RefFasta} \
            -I ${Tumorbam} -tumor ${TumorReplicateId} \
            -L ${IntervalsList} --native-pair-hmm-threads ${task.cpus} \
            -O ${IntervalsList}.vcf.gz
        """
    }

    process 'MergeTumor' {
    // Merge scattered Mutect2 vcfs

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/", mode: params.publishDirMode

        input:
        file(vcf) from Mutect2VcfTumor.collect()
        file(stats) from Mutect2StatsTumor.collect()
        set TumorReplicateId, file(idx) from Mutect2IdxTumor.groupTuple()

        output:
        set TumorReplicateId, file("${TumorReplicateId}_mutect2_raw.vcf.gz"),
        file("${TumorReplicateId}_mutect2_raw.vcf.gz.tbi"),
        file("${TumorReplicateId}_mutect2_raw.vcf.gz.stats") into mutect2Tumor

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 MergeVcfs \
            --tmp-dir ${params.tmpDir} \
            -I ${vcf.join(" -I ")} \
            -O ${TumorReplicateId}_mutect2_raw.vcf.gz

        $GATK4 MergeMutectStats \
            --tmp-dir ${params.tmpDir} \
            --stats ${stats.join(" --stats ")} \
            -O ${TumorReplicateId}_mutect2_raw.vcf.gz.stats
        """
    }

    process 'FilterMutec2Tumor' {
    /*
     CalculateContamination (GATK4): calculate fraction of reads coming from
     cross-sample contamination
     FilterMutectCalls (GATK4): filter somatic SNVs and indels
     FilterByOrientationBias (GATK4): filter variant calls using orientation bias
     SelectVariants (GATK4): select subset of variants from a larger callset
     VariantFiltration (GATK4): filter calls based on INFO and FORMAT annotations
    */

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/", mode: params.publishDirMode

        input:
        set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict])
        set TumorReplicateId, file(pileup) from PileupTumor1
        set TumorReplicateId, file(vcf), file(vcfIdx), file(vcfStats) from mutect2Tumor
        file(preAdapterDetail) from CollectSequencingArtifactMetrics1

        output:
        set TumorReplicateId, file("${TumorReplicateId}_mutect2_final.vcf.gz"),
            file("${TumorReplicateId}_mutect2_final.vcf.gz.tbi") into FilterMutect2Tumor

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 CalculateContamination \
            --tmp-dir ${params.tmpDir} \
            -I ${pileup} \
            -O ${TumorReplicateId}_cont.table && \
        $GATK4 FilterMutectCalls \
            --tmp-dir ${params.tmpDir} \
            -R ${RefFasta} \
            -V ${vcf} \
            --contamination-table ${TumorReplicateId}_cont.table \
            -O ${TumorReplicateId}_oncefiltered.vcf.gz && \
        $GATK4 FilterByOrientationBias \
            --tmp-dir ${params.tmpDir} \
            -V ${TumorReplicateId}_oncefiltered.vcf.gz \
            -P ${preAdapterDetail} \
            -O ${TumorReplicateId}_twicefitlered.vcf.gz && \
        $GATK4 SelectVariants \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_twicefitlered.vcf.gz \
            -R ${RefFasta} \
            --exclude-filtered true \
            --output ${TumorReplicateId}_mutect2_pass.vcf && \
        $GATK4 VariantFiltration \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_mutect2_pass.vcf \
            -R ${RefFasta} \
            --genotype-filter-expression 'g.getAD().1 < ${params.minAD}' \
            --genotype-filter-name "AD.1_${params.minAD}" \
            --output ${TumorReplicateId}_mutect2_final.vcf.gz
        """
    }
}

if (readsControl != 'NO_FILE') {
    process 'BaseRecalApplyControl' {
    /*
     BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
     Recalibration (BQSR)
     ApplyBQSR (GATK4): apply BQSR table to reads
    */

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/03_mutect2/processing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesControl1
        set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
        set file(MillsGold), file(MillsGoldIdx), file(DBSNP), file(DBSNPIdx),
         file(KnownIndels), file(KnownIndelsIdx) from Channel.value(
                    [database.MillsGold, database.MillsGoldIdx,
                    database.DBSNP, database.DBSNPIdx,
                    database.KnownIndels, database.KnownIndelsIdx])

        output:
        set ControlReplicateId,
         file("${ControlReplicateId}_Control_bqsr.table") into BaseRecalibratorControl
        set ControlReplicateId, file("${ControlReplicateId}_Control_recal4.bam"),
         file("${ControlReplicateId}_Control_recal4.bam.bai") into ApplyControl1, ApplyControl2

        script:
        """
        mkdir -p ${params.tmpDir}

        java -jar $PICARD SetNmMdAndUqTags \
            TMP_DIR=${params.tmpDir} \
            R=${RefFasta} \
            I=${bam} \
            O=Control_fixed.bam \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT && \
        $GATK4 BaseRecalibratorSpark \
            --tmp-dir ${params.tmpDir} \
            -I Control_fixed.bam \
            -R ${RefFasta} \
            -L ${IntervalsList} \
            -O ${ControlReplicateId}_Control_bqsr.table \
            --known-sites ${DBSNP} \
            --known-sites ${KnownIndels} \
            --known-sites ${MillsGold} \
            --conf 'spark.executor.cores=${task.cpus}'&& \
         $GATK4 ApplyBQSRSpark \
            --tmp-dir ${params.tmpDir} \
            -I ${bam} \
            -R ${RefFasta} \
            -L ${IntervalsList} \
            -O ${ControlReplicateId}_Control_recal4.bam \
            --bqsr-recal-file ${ControlReplicateId}_Control_bqsr.table \
            --conf 'spark.executor.cores=${task.cpus}'
        """
    }

    process 'GetPileupControl' {
    // GetPileupSummaries (GATK4): tabulates pileup metrics for inferring contamination

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/03_mutect2/processing/",
         mode: params.publishDirMode

        input:
        set file(GnomAD), file(GnomADIdx) from Channel.value(
            [database.GnomAD, database.GnomADIdx])
        file(intervals) from Channel.value([reference.IntervalsList])
        set ControlReplicateId, file(bam), file(bai) from ApplyControl1

        output:
        set ControlReplicateId, file("${ControlReplicateId}_Control_pileup.table") into PileupControl

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 GetPileupSummaries \
            --tmp-dir ${params.tmpDir} \
            -I ${bam} \
            -O ${ControlReplicateId}_Control_pileup.table \
            -L ${intervals} \
            --variant ${GnomAD}
        """
    }

    process 'Mutect2' {
    /* 
     Call somatic SNPs and indels via local re-assembly of haplotypes; tumor sample
     and matched normal sample 
    */

        tag "$TumorreplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
         mode: params.publishDirMode

        input:
        set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict])
        set TumorReplicateId, file(Tumorbam), file(Tumorbai),
            ControlReplicateId, file(Controlbam), file(Controlbai),
            file(intervals) from ApplyTumor5.combine(ApplyControl2).combine(interval_ch2.flatten())

        output:
        file("${intervals}.vcf.gz") into Mutect2Vcf
        file("${intervals}.vcf.gz.stats") into Mutect2Stats
        set TumorReplicateId, ControlReplicateId, file("${intervals}.vcf.gz.tbi") into Mutect2Idx

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 Mutect2 \
            --tmp-dir ${params.tmpDir} \
            -R ${RefFasta} \
            -I ${Tumorbam} -tumor ${TumorReplicateId} \
            -I ${Controlbam} -normal ${ControlReplicateId} \
            -L ${intervals} --native-pair-hmm-threads ${task.cpus} \
            -O ${intervals}.vcf.gz
        """
    }

    process 'Merge' {
    // Merge scattered Mutect2 vcfs

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/",
         mode: params.publishDirMode

        input:
        file(vcf) from Mutect2Vcf.collect()
        file (stats) from Mutect2Stats.collect()
        set TumorReplicateId, ControlReplicateId, file(idx) from Mutect2Idx.groupTuple(by: [0,1])

        output:
        set TumorReplicateId, ControlReplicateId,
         file("${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz"),
             file("${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz.tbi"),
             file("${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz.stats") into mutect2

        script:
        """
        mkdir -p ${params.tmpDir}
        
        $GATK4 MergeVcfs \
            --tmp-dir ${params.tmpDir} \
            -I ${vcf.join(" -I ")} \
            -O ${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz

        $GATK4 MergeMutectStats \
            --tmp-dir ${params.tmpDir} \
            --stats ${stats.join(" --stats ")} \
            -O ${TumorReplicateId}_${ControlReplicateId}_mutect2_raw.vcf.gz.stats
        """
    }

    process 'FilterMutect2' {
    /* 
     CalculateContamination (GATK4): calculate fraction of reads coming from
     cross-sample contamination
     FilterMutectCalls (GATK4): filter somatic SNVs and indels
     FilterByOrientationBias (GATK4): filter variant calls using orientation bias
     SelectVariants (GATK4): select subset of variants from a larger callset
     VariantFiltration (GATK4): filter calls based on INFO and FORMAT annotations
    */

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/",
         mode: params.publishDirMode

        input:
        set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict])
        set TumorReplicateId, file(pileupTumor) from PileupTumor2
        set ControlReplicateId, file(pileupControl) from PileupControl
        set TumorReplicateId, ControlReplicateId, file(vcf), file(vcfIdx),
         file(vcfStats) from mutect2
        file(preAdapterDetail) from CollectSequencingArtifactMetrics2

        output:
        set TumorReplicateId, ControlReplicateId,
            file("${TumorReplicateId}_${ControlReplicateId}_mutect2_final.vcf.gz"),
            file("${TumorReplicateId}_${ControlReplicateId}_mutect2_final.vcf.gz.tbi") into FilterMutect2

        script:
        """
        mkdir -p ${params.tmpDir}

        $GATK4 CalculateContamination \
            --tmp-dir ${params.tmpDir} \
            -I ${pileupTumor} \
            --matched-normal ${pileupControl} \
            -O ${TumorReplicateId}_${ControlReplicateId}_cont.table && \
        $GATK4 FilterMutectCalls \
            --tmp-dir ${params.tmpDir} \
            -R ${RefFasta} \
            -V ${vcf} \
            --contamination-table ${TumorReplicateId}_${ControlReplicateId}_cont.table \
            -O ${TumorReplicateId}_${ControlReplicateId}_oncefiltered.vcf.gz && \
        $GATK4 FilterByOrientationBias \
            --tmp-dir ${params.tmpDir} \
            -V ${TumorReplicateId}_${ControlReplicateId}_oncefiltered.vcf.gz \
            -P ${preAdapterDetail} \
            -O ${TumorReplicateId}_${ControlReplicateId}_twicefitlered.vcf.gz && \
        $GATK4 SelectVariants \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_${ControlReplicateId}_twicefitlered.vcf.gz \
            -R ${RefFasta} \
            --exclude-filtered true \
            --output ${TumorReplicateId}_${ControlReplicateId}_mutect2_pass.vcf && \
        $GATK4 VariantFiltration \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_${ControlReplicateId}_mutect2_pass.vcf \
            -R ${RefFasta} \
            --genotype-filter-expression 'g.getAD().1 < ${params.minAD}' \
            --genotype-filter-name "AD.1_${params.minAD}" \
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
/* 
 RealignerTargetCreator (GATK3): define intervals to target for local realignment
 IndelRealigner (GATK3): perform local realignment of reads around indels
*/

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/05_varscan/processing/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesTumor2
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(KnownIndels), file(KnownIndelsIdx), file(MillsGold),
     file(MillsGoldIdx) from Channel.value([database.KnownIndels, database.KnownIndelsIdx,
            database.MillsGold, database.MillsGoldIdx])

    output:
    set TumorReplicateId, file("${TumorReplicateId}_aligned_sort_mkdp_realign.bam"),
     file("${TumorReplicateId}_aligned_sort_mkdp_realign.bai") into IndelRealignerTumor

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \
        -T RealignerTargetCreator \
        --known ${MillsGold} \
        --known ${KnownIndels} \
        -R ${RefFasta} \
        -L ${IntervalsList} \
        -I ${bam} \
        -o target.list \
        -nt ${task.cpus} && \
    $JAVA8 -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \
        -T IndelRealigner \
        -R ${RefFasta} \
        -L ${IntervalsList} \
        -I ${bam} \
        -targetIntervals target.list \
        -known ${KnownIndels} \
        -known ${MillsGold} \
        -nWayOut _realign.bam && \
    rm target.list
    """
}

process 'FixMateBaseRecalTumor' {
/*
 FixMateInformation (Picard): verify mate-pair information between mates; ensure that
 all mate-pair information is in sync between reads and its pairs
*/

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/05_varscan/processing/",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam), file(bai) from IndelRealignerTumor
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
    set file(DBSNP), file(DBSNPIdx), file(KnownIndels), file(KnownIndelsIdx),
        file(MillsGold), file(MillsGoldIdx) from Channel.value(
        [database.DBSNP, database.DBSNPIdx, database.KnownIndels,
        database.KnownIndelsIdx, database.MillsGold, database.MillsGoldIdx])

    output:
    set TumorReplicateId, file("${TumorReplicateId}_fixmate.bam"),
     file("${TumorReplicateId}_fixmate.bai"),
     file("${TumorReplicateId}_bqsr3.table") into FixMateBaseRecalTumor

    script:
    """
    mkdir -p ${params.tmpDir}

    java -XX:ParallelGCThreads=${task.cpus} -jar $PICARD FixMateInformation \
        TMP_DIR=${params.tmpDir} \
        I=${bam} \
        O=${TumorReplicateId}_fixmate.bam \
        CREATE_INDEX=true && \
    $JAVA8 -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \
        -T BaseRecalibrator \
        -R ${RefFasta} \
        -L ${IntervalsList} \
        -I ${TumorReplicateId}_fixmate.bam \
        -knownSites ${DBSNP} \
        -knownSites ${KnownIndels} \
        -knownSites ${MillsGold} \
        -nct ${task.cpus} \
        -o ${TumorReplicateId}_bqsr3.table
    """
}

process 'PrintReadsTumor' {
// PrintReads (GATK3): write out sequence read data for filtering, merging, subsetting

    tag "$TumorReplicateId"
    publishDir "$params.outputDir/$TumorReplicateId/05_varscan/processing",
     mode: params.publishDirMode

    input:
    set TumorReplicateId, file(bam), file(bai), file(bqsr3) from FixMateBaseRecalTumor
    set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
        [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])

    output:
    set TumorReplicateId, file("${TumorReplicateId}_printreads.bam"),
     file("${TumorReplicateId}_printreads.bai") into PrintReadsTumor1, PrintReadsTumor2,
        PrintReadsTumor3

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

if (readsControl != "NO_FILE") {

    process 'IndelRealignerControl' {
    /*
     RealignerTargetCreator (GATK3): define intervals to target for local realignment
     IndelRealigner (GATK3): perform local realignment of reads around indels
    */

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/05_varscan/processing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(bam), file(bai), file(list) from MarkDuplicatesControl2
        set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
        set file(KnownIndels), file(KnownIndelsIdx), file(MillsGold),
            file(MillsGoldIdx) from Channel.value(
            [database.KnownIndels, database.KnownIndelsIdx,
            database.MillsGold, database.MillsGoldIdx])

        output:
        set ControlReplicateId, file("${ControlReplicateId}_Control_aligned_sort_mkdp_realign.bam"),
         file("${ControlReplicateId}_Control_aligned_sort_mkdp_realign.bai") into IndelRealignerControl

        script:
        """
        mkdir -p ${params.tmpDir}

        $JAVA8 -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \
            -T RealignerTargetCreator \
            --known ${MillsGold} \
            --known ${KnownIndels} \
            -R ${RefFasta} \
            -L ${IntervalsList} \
            -I ${bam} \
            -o target.list \
            -nt ${task.cpus} && \
        $JAVA8 -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \
            -T IndelRealigner \
            -R ${RefFasta} \
            -L ${IntervalsList} \
            -I ${bam} \
            -targetIntervals target.list \
            -known ${KnownIndels} \
            -known ${MillsGold} \
            -nWayOut _realign.bam && \
        rm target.list
        """
    }

    process 'FixMateBaseRecalControl' {
    /*
     FixMateInformation (Picard): verify mate-pair information beteen mates; ensure that
     all mate-pair information is in sync between reads and its pairs
    */

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/05_varscan/processing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(bam), file(bai) from IndelRealignerControl
        set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
        set file(DBSNP), file(DBSNPIdx), file(KnownIndels), file(KnownIndelsIdx),
            file(MillsGold), file(MillsGoldIdx) from Channel.value(
            [database.DBSNP, database.DBSNPIdx, database.KnownIndels, database.KnownIndelsIdx,
            database.MillsGold, database.MillsGoldIdx])

        output:
        set ControlReplicateId, file("${ControlReplicateId}_fixmate.bam"),
            file("${ControlReplicateId}_Control_fixmate.bai"),
             file("${ControlReplicateId}_Control_bqsr3.table") into FixMateBaseRecalControl

        script:
        """
        mkdir -p ${params.tmpDir}
        java -XX:ParallelGCThreads=${task.cpus} -jar $PICARD FixMateInformation \
            TMP_DIR=${params.tmpDir} \
            I=${bam} \
            O=${ControlReplicateId}_fixmate.bam \
            CREATE_INDEX=true && \
            $JAVA8 -jar $GATK3\
            -T BaseRecalibrator \
            -R ${RefFasta} \
            -L ${IntervalsList} \
            -I ${ControlReplicateId}_Control_fixmate.bam \
            -knownSites ${DBSNP} \
            -knownSites ${KnownIndels} \
            -knownSites ${MillsGold} \
            -nct ${task.cpus} \
            -o ${ControlReplicateId}_Control_bqsr3.table
        """
    }

    process 'PrintReadsControl' {
    // PrintReads (GATK3): write out sequence read data for filtering, merging, subsetting

        tag "$ControlReplicateId"
        publishDir "$params.outputDir/$ControlReplicateId/05_varscan/processing/",
         mode: params.publishDirMode

        input:
        set ControlReplicateId, file(bam), file(bai), file(bqsr3) from FixMateBaseRecalControl
        set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])

        output:
        set ControlReplicateId, file("${ControlReplicateId}_Control_printreads.bam"),
         file("${ControlReplicateId}_Control_printreads.bai") into PrintReadsControl1,
            PrintReadsControl2

        script:
        """
        mkdir -p ${params.tmpDir}

        $JAVA8 -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3  \
            -T PrintReads \
            -R ${RefFasta} \
            -L ${IntervalsList} \
            -I ${bam} \
            -BQSR ${bqsr3} \
            -nct ${task.cpus} \
            -o ${ControlReplicateId}_Control_printreads.bam
        """
    }

    process 'ReorderSam' {
    // ReorderSam (Picard): reorders reads to match the contig order in a reference file

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/05_varscan/processing/",
         mode: params.publishDirMode

        input:
        set TumorReplicateId, file(bamTumor), file(baiTumor) from PrintReadsTumor1
        set ControlReplicateId, file(bamControl), file(baiControl) from PrintReadsControl1
        set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict])

        output:
        set TumorReplicateId, file("${TumorReplicateId}_reorder_tumor.bam"),
         file("${TumorReplicateId}_reorder_tumor.bai") into ReorderSamTumor1, ReorderSamTumor2
        set ControlReplicateId, file("${ControlReplicateId}_reorder_control.bam"),
         file("${ControlReplicateId}_reorder_control.bai") into ReorderSamControl

        script:
        """
        mkdir -p ${params.tmpDir}
        java -jar $PICARD ReorderSam \
            TMP_DIR=${params.tmpDir} \
            SD=${RefDict} \
            INPUT=${bamTumor} \
            CREATE_INDEX=true \
            OUTPUT=${TumorReplicateId}_reorder_tumor.bam &
        tumor_process_id=\$!

        java -jar $PICARD ReorderSam \
            TMP_DIR=${params.tmpDir} \
            SD=${RefDict} \
            INPUT=${bamControl} \
            CREATE_INDEX=true \
            OUTPUT=${ControlReplicateId}_reorder_control.bam &
            control_process_id=\$!

        wait $tumor_process_id
        wait $control_process_id
        """
    }

    process 'VarscanSomatic' {
    // somatic (Varscan): calls somatic variants (SNPS and indels)

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/05_varscan/", mode: params.publishDirMode

        input:
        set TumorReplicateId, file(bamTumor), file(baiTumor) from ReorderSamTumor1
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
        --p-value ${params.somatic_pvalue} \
        --somatic-p-value ${params.somatic_somaticpvalue} \
        --strand-filter ${params.strand_filter} && \
        rm -f ${TumorReplicateId}_${ControlReplicateId}_mpileup.fifo
        """
    }

    process 'ProcessSomatic' {
    // Filter variants by somatic status and confidences

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/05_varscan/processing",
         mode: params.publishDirMode

        input:
        set TumorReplicateId, ControlReplicateId, file(snp), file(indel) from VarscanSomatic

        output:
        set TumorReplicateId, ControlReplicateId,
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.hc.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.LOH.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.LOH.hc.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Germline.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Germline.hc.vcf") into ProcessSomaticSNP
        set TumorReplicateId, ControlReplicateId,
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Somatic.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Somatic.hc.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.LOH.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.LOH.hc.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Germline.vcf"),
            file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Germline.hc.vcf") into ProcessSomaticIndel

        script:
        """
        java -jar $VARSCAN processSomatic \
        ${snp} \
        --min-tumor-freq ${params.min_tumor_freq} \
        --max-normal-freq ${params.max_normal_freq} \
        --p-value ${params.processSomatic_pvalue} && \
        java -jar $VARSCAN processSomatic \
        ${indel} \
        --min-tumor-freq ${params.min_tumor_freq} \
        --max-normal-freq ${params.max_normal_freq} \
        --p-value ${params.processSomatic_pvalue}
        """
    }

    process 'FilterVarscan' {
    /*
     AWK-script: calcualtes start-end position of variant
     Bamreadcount: generate metrics at single nucleotide positions for filtering
     fpfilter (Varscan): apply false-positive filter to variants
    */

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/05_varscan/",
         mode: params.publishDirMode

        input:
        set TumorReplicateId, file(bam), file(bai) from ReorderSamTumor2
        set TumorReplicateId, ControlReplicateId, file(snpSomatic), file(snpSomaticHc),
            file(snpLOH), file(snpLOHhc), file(snpGerm), file(snpGemHc) from ProcessSomaticSNP
        set TumorReplicateId, ControlReplicateId, file(indelSomatic), file(indelSomaticHc),
            file(indelLOH), file(indelLOHhc), file(indelGerm),
            file(indelGemHc) from ProcessSomaticIndel
        set file(RefFasta), file(RefIdx), file(RefDict) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict])

        output:
        set TumorReplicateId, ControlReplicateId,
         file("${TumorReplicateId}_${ControlReplicateId}_varscan.snp.Somatic.hc.filtered.vcf") into FilterSnp
        set TumorReplicateId, ControlReplicateId,
         file("${TumorReplicateId}_${ControlReplicateId}_varscan.indel.Somatic.hc.filtered.vcf") into FilterIndel

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

if (readsControl != "NO_FILE") {
    process 'MutectTumorNormal' {
    // Mutect1: calls SNPS from tumor and matched normal sample

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/04_mutect1/",
         mode: params.publishDirMode

        input:
        set TumorReplicateId, file(bamTumor), file(baiTumor) from PrintReadsTumor2
        set ControlReplicateId, file(bamControl), file(baiControl) from PrintReadsControl2
        set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
        set file(DBSNP), file(DBSNPIdx), file(Cosmic), file(CosmicIdx) from Channel.value(
            [database.DBSNP, database.DBSNPIdx, database.Cosmic, database.CosmicIdx])

        output:
        set TumorReplicateId, ControlReplicateId,
            file("${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.vcf.gz"),
            file("${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.vcf.gz.idx"),
            file("${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.stats.txt") into raw
        set TumorReplicateId, ControlReplicateId,
        file("${TumorReplicateId}_${ControlReplicateId}_mutect1_final.vcf.gz"),
        file("${TumorReplicateId}_${ControlReplicateId}_mutect1_final.vcf.gz.tbi") into Mutect1filter

        script:
        """
        mkdir -p ${params.tmpDir}

        $JAVA7 -Djava.io.tmpdir=${params.tmpDir} -jar $MUTECT1 \
            --analysis_type MuTect \
            --reference_sequence ${RefFasta} \
            --cosmic ${Cosmic} \
            --dbsnp ${DBSNP} \
            -L ${IntervalsList} \
            --input_file:normal ${bamControl} \
            --input_file:tumor ${bamTumor} \
            --out ${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.stats.txt \
            --vcf ${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.vcf.gz && \
        $GATK4 SelectVariants \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_${ControlReplicateId}_mutect1.raw.vcf.gz \
            -R ${RefFasta} \
            --exclude-filtered true \
            --output ${TumorReplicateId}_${ControlReplicateId}_mutect1_pass.vcf && \
        $GATK4 VariantFiltration \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_${ControlReplicateId}_mutect1_pass.vcf \
            -R ${RefFasta} \
            --genotype-filter-expression "g.getAD().1 < 2" \
            --genotype-filter-name "AD.1_2" \
            --output ${TumorReplicateId}_${ControlReplicateId}_mutect1_final.vcf.gz
        """
    }
}else {
    process 'MutectTumor' {
    // Mutect1: calls SNPS from tumor sample

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/04_mutect1/", mode: params.publishDirMode

        input:
        set TumorReplicateId, file(bamTumor), file(baiTumor) from PrintReadsTumor3
        set file(RefFasta), file(RefIdx), file(RefDict), file(IntervalsList) from Channel.value(
            [reference.RefFasta, reference.RefIdx, reference.RefDict, reference.IntervalsList])
        set file(DBSNP), file(DBSNPIdx), file(Cosmic), file(CosmicIdx) from Channel.value(
            [database.DBSNP, database.DBSNPIdx, database.Cosmic, database.CosmicIdx])

        output:
        set TumorReplicateId,
            file("${TumorReplicateId}_mutect1.raw.vcf.gz"),
            file("${TumorReplicateId}_mutect1.raw.vcf.gz.idx")
            file("${TumorReplicateId}_mutect1.raw.stats.txt") into Mutect1rawTumor
        set TumorReplicateId,
            file("${TumorReplicateId}_mutect1_final.vcf.gz"),
            file("${TumorReplicateId}_mutect1_final.vcf.gz.tbi") into Mutect1filterTumor

        script:
        """
        mkdir -p ${params.tmpDir}

        $JAVA7 -Djava.io.tmpdir=${params.tmpDir} -jar $MUTECT1 \
            --analysis_type MuTect \
            --reference_sequence ${RefFasta} \
            --cosmic ${Cosmic} \
            --dbsnp ${DBSNP} \
            -L ${IntervalsList} \
            --input_file:tumor ${bamTumor} \
            --out ${TumorReplicateId}_mutect1.raw.stats.txt \
            --vcf ${TumorReplicateId}_mutect1.raw.vcf.gz && \
        $GATK4 SelectVariants \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_mutect1.raw.vcf.gz \
            -R ${RefFasta} \
            --exclude-filtered true \
            --output ${TumorReplicateId}_mutect1_pass.vcf && \
        $GATK4 VariantFiltration \
            --tmp-dir ${params.tmpDir} \
            --variant ${TumorReplicateId}_mutect1_pass.vcf \
            -R ${RefFasta} \
            --genotype-filter-expression "g.getAD().1 < 2" \
            --genotype-filter-name "AD.1_2" \
            --output ${TumorReplicateId}_mutect1_final.vcf.gz
        """
    }
}

if (readsControl != "NO_FILE") {
    process 'Vep' {
    //

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/06_vep/", mode: params.publishDirMode

        input:
        set TumorReplicateId, ControlReplicateId, file(mutect2Vcf),
            file(mutect2Idx) from FilterMutect2
        set TumorReplicateId, ControlReplicateId, file(varscanSnp) from FilterSnp
        set TumorReplicateId, ControlReplicateId, file(varscanIndel) from FilterIndel
        set TumorReplicateId, ControlReplicateId, file(mutect1Vcf),
        file(mutect1Idx) from Mutect1filter
        file(VepFasta) from Channel.value([reference.VepFasta])

        script:
        """
        $PERL $VEP -i ${mutect1Vcf} \
        -o ${TumorReplicateId}_${ControlReplicateId}_mutect1_vep.txt \
        --fork ${params.vep_threads} \
        --stats_file ${TumorReplicateId}_${ControlReplicateId}_mutect1_vep_summary.html \
        --species ${params.vep_species} \
        --assembly ${params.vep_assembly} \
        --offline \
        --dir ${params.vep_dir} \
        --cache --dir_cache ${params.vep_cache} \
        --fasta ${VepFasta} \
        --format "vcf" \
        --tab  && \
        $PERL $VEP -i ${mutect2Vcf} \
        -o ${TumorReplicateId}_${ControlReplicateId}_mutect2_vep.txt \
        --fork ${params.vep_threads} \
        --stats_file ${TumorReplicateId}_${ControlReplicateId}_mutect2_vep_summary.html \
        --species ${params.vep_species} \
        --assembly ${params.vep_assembly} \
        --offline \
        --dir ${params.vep_dir} \
        --cache --dir_cache ${params.vep_cache} \
        --fasta ${VepFasta} \
        --format "vcf" \
        --tab  && \
        $PERL $VEP -i ${varscanSnp} \
        -o ${TumorReplicateId}_${ControlReplicateId}_varscanSnp_vep.txt \
        --fork ${params.vep_threads} \
        --stats_file ${TumorReplicateId}_${ControlReplicateId}_varscanSnp_vep_summary.html \
        --species ${params.vep_species} \
        --assembly ${params.vep_assembly} \
        --offline \
        --dir ${params.vep_dir} \
        --cache --dir_cache ${params.vep_cache} \
        --format "vcf" \
        --tab  && \
        $PERL $VEP -i ${varscanIndel} \
        -o ${TumorReplicateId}_${ControlReplicateId}_varscanIndel_vep.txt \
        --fork ${params.vep_threads} \
        --stats_file ${TumorReplicateId}_${ControlReplicateId}_varscanIndel_vep_summary.html \
        --species ${params.vep_species} \
        --assembly ${params.vep_assembly} \
        --offline \
        --dir ${params.vep_dir} \
        --cache --dir_cache ${params.vep_cache} \
        --format "vcf" \
        --tab && \
        variant_sort.sh  \
        -1 ${TumorReplicateId}_${ControlReplicateId}_mutect1_vep.txt \
        -2 ${TumorReplicateId}_${ControlReplicateId}_mutect2_vep.txt \
        -S ${TumorReplicateId}_${ControlReplicateId}_varscanSnp_vep.txt \
        -I ${TumorReplicateId}_${ControlReplicateId}_varscanIndel_vep.txt \
        -t ${TumorReplicateId} \
        -c ${ControlReplicateId}
        """
    }
} else {
    process 'VepTumor' {
    //

        tag "$TumorReplicateId"
        publishDir "$params.outputDir/$TumorReplicateId/06_vep/", mode: params.publishDirMode

        input:
        set TumorReplicateId, file(mutect2Vcf),
            file(mutect2Idx) from FilterMutect2Tumor
        set TumorReplicateId, file(mutect1Vcf),
        file(mutect1Idx) from Mutect1filterTumor
        file(VepFasta) from Channel.value([reference.VepFasta])

        script:
        """
        $PERL $VEP -i ${mutect1Vcf} \
        -o ${TumorReplicateId}_mutect1_vep.txt \
        --fork ${params.vep_threads} \
        --stats_file ${TumorReplicateId}_mutect1_vep_summary.html \
        --species ${params.vep_species} \
        --assembly ${params.vep_assembly} \
        --offline \
        --dir ${params.vep_dir} \
        --cache --dir_cache ${params.vep_cache} \
        --fasta ${VepFasta} \
        --format "vcf" \
        --tab  && \
        $PERL $VEP -i ${mutect2Vcf} \
        -o ${TumorReplicateId}_mutect2_vep.txt \
        --fork ${params.vep_threads} \
        --stats_file ${TumorReplicateId}_mutect2_vep_summary.html \
        --species ${params.vep_species} \
        --assembly ${params.vep_assembly} \
        --offline \
        --dir ${params.vep_dir} \
        --cache --dir_cache ${params.vep_cache} \
        --fasta ${VepFasta} \
        --format "vcf" \
        --tab  && \
        variant_sort.sh \
        -1 ${TumorReplicateId}_mutect1_vep.txt \
        -2 ${TumorReplicateId}_mutect2_vep.txt \
        -t ${TumorReplicateId}
        """
    }
}

/*
________________________________________________________________________________

                            F U N C T I O N S
________________________________________________________________________________

*/

def mkTmpDir() {
    myTmpDir = file(params.tmpDir)
    result = myTmpDir.mkdirs()
    println result ? "tmpDir created: $myTmpDir" : "Cannot create directory: $myTmpDir"
}

def checkParamReturnFileReferences(item) {
    params."${item}" = params.references."${item}"
    return file(params."${item}")
}

def checkParamReturnFileDatabases(item) {
    params."${item}" = params.databases."${item}"
    return file(params."${item}")
}

def defineReference() {
    if (params.references.size() != 8) exit 1, """
    ERROR: Not all References needed found in configuration
    Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
    """
    return [
        'RefFasta'   : checkParamReturnFileReferences("RefFasta"),
        'RefIdx'       : checkParamReturnFileReferences("RefIdx"),
        'RefDict'    : checkParamReturnFileReferences("RefDict"),
        'BwaRef'     : checkParamReturnFileReferences("BwaRef"),
        'IntervalsList' : checkParamReturnFileReferences("IntervalsList"),
        'BaitIntervalsList' : checkParamReturnFileReferences("BaitIntervalsList"),
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
    'MillsGold'   : checkParamReturnFileDatabases("MillsGold"),
    'MillsGoldIdx'   : checkParamReturnFileDatabases("MillsGoldIdx"),
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

def helpMessage() {
    log.info ""
    log.info "----------------------------"
    log.info "--        U S A G E       --"
    log.info "----------------------------"
    log.info ""
    log.info ' nextflow run wes.nf "--readsTumor|--batchFile" "[--readsControl]" "--IntervalsList" "--IntervalsBed" [--single_end]'
    log.info ""
    log.info "-------------------------------------------------------------------------"
    log.info ""
    log.info "Single-end reads:"
    log.info "---------------------------"
    log.info "--single_end \t\t sets parameter to TRUE (default false)"
    log.info ""
    log.info " Mandatory arguments:"
    log.info " --------------------"
    log.info "--readsTumor \t\t reads_{1,2}.fastq \t\t paired-end reads; FASTQ files (can be zipped)"
    log.info "   or"
    log.info "--batchFile"
    log.info "CSV-file, paired-end T/N reads:"
    log.info "sampleId,readsTumorFWD,readsTumorREV,readsControlFWD,readsControlREV,group"
    log.info "sample1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,Control1_reads_1.fastq,Control1_reads_2.fastq,group1"
    log.info "sample2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,Control2_reads_1.fastq,Control2_reads_2.fastq,group1"
    log.info "..."
    log.info "sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,ControlN_reads_1.fastq,ControlN_reads_2.fastq,groupX"
    log.info ""
    log.info "CSV-file, single-end T only reads:"
    log.info "sampleId,readsTumorFWD,readsTumorREV,readsControlFWD,readsControlREV,group"
    log.info "sample1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,NO_FILE,NO_FILE,group1"
    log.info "sample2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,NO_FILE,NO_FILE,group1"
    log.info "..."
    log.info "sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,NO_FILE,,NO_FILE,groupX"
    log.info ""
    log.info "FASTQ files (can be zipped), if single-end reads are used put NO_FILE instead of *_reads_2.fastq in the REV fields"
    log.info ""
    log.info "--IntervalsList \t intervals.list \t\t interval.list file for targeting"
    log.info "--IntervalsBed \t\t intervals.bed \t\t\t interval.bed file for targeting"
    log.info ""
    log.info " Optional argument:"
    log.info " ------------------"
    log.info "--readsControl \t\t reads_{1,2}.fastq \t\t paired-end reads; FASTQ file (can be zipped)"
    log.info "--sampleName \t\t  sample name. If not specified samples will be named according to the fastq filenames."
    log.info ""
    log.info " All references, databases, software should be edited in the nextflow.config file"
    log.info ""
    log.info " Mandatory references:"
    log.info " ---------------------"
    log.info " RefFasta \t\t reference.fa \t\t\t Reference Genome; FASTA file"
    log.info " RefIdx \t\t\t reference.fai \t\t\t Referene Genome Index, FAI file"
    log.info " RefDict \t\t reference.dict \t\t Reference Genome Dictionary, DICT file"
    log.info " BwaRef \t\t\t reference.{amb,sa,ann,pac,bwt}  Reference Genome prepared for BWA mem"
    log.info " VepFasta \t\t reference.toplevel.fa \t\t Reference genome; ENSEMBL toplevel"
    log.info ""
    log.info " Mandatory databases:"
    log.info " --------------------"
    log.info " MillsGold/Idx \t\t Mills_and_1000G_gold_standard.indels.vcf/idx \t Gold standard Indels database, VCF file, IDX File"
    log.info " Cosmic/Idx \t\t CosmicCodingMuts.vcf \t\t\t\t Cosmic conding mutations, VCF file, IDX file"
    log.info " DBSNP/Idx \t\t Homo_sapiens_assembly.dbsnp.vcf/idx \t\t SNPS, microsatellites, and small-scale insertions and deletions, VCF file, IDX file"
    log.info " GnomAD/Idx \t\t small_exac_common_3.vcf/idx \t\t\t exonix sites only for contamination estimation from GATK, VCF file, IDX file"
    log.info " KnownIdenls/Idx \t Homo_sapiens_assembly.known_indels.vcf/idx \t Known Indels from GATK resource Bundle, VCF file, IDX file"
    log.info ""
    log.info " Required software:"
    log.info " ------------------"
    log.info " JAVA7 \t\t\t Version 1.7"
    log.info " JAVA8 \t\t\t Version 1.8"
    log.info " BWA \t\t\t Version 0.7.17"
    log.info " SAMTOOLS \t\t Version 1.9"
    log.info " PICARD \t\t\t Version 2.16.0"
    log.info " GATK3 \t\t\t Version 3.8-0"
    log.info " GATK4 \t\t\t Version 4.1.3.0"
    log.info " VARSCAN \t\t Version 2.4.3"
    log.info " MUTECT1 \t\t Version 1.1.7"
    log.info " BAMREADCOUNT \t\t Version 0.8.0"
    log.info " VEP \t\t\t Version 2.0"
    log.info "-------------------------------------------------------------------------"
}
