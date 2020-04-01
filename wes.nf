#!/usr/bin/env nextflow

log.info ""
log.info " NEXTFLOW ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
log.info "-------------------------------------------------------------------------"
log.info "         W H O L E - E X O M E   P I P E L I N E             "
log.info "-------------------------------------------------------------------------"
log.info ""
log.info " Finding SNPs and Indels in tumor only or tumor + matched normal samples"
log.info " using 5 different caller: \n \t * MuTect1 \n \t * MuTect2 \n \t * VarScan2 \n \t * Manta \n \t * Strelka"
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
if  (params.readsNormal != "NO_FILE") log.info " Reads Tumor: \t\t " + params.readsTumor
if  (params.readsNormal != "NO_FILE") log.info " Reads Normal: \t\t " + params.readsNormal
log.info ""
log.info "Please check --help for further instruction"
log.info "-------------------------------------------------------------------------"

/*
________________________________________________________________________________

                            C O N F I G U R A T I O N
________________________________________________________________________________
*/
if (params.help) exit 0, helpMessage()

// switch for enable/disable processes (debug/devel only: use if(params.RUNTHIS) { .... })
params.RUNTHIS = false

// default is not to process a batchfile
params.batchFile = false

// set single_end variable to supplied param
single_end = params.single_end
single_end_RNA = params.single_end_RNA

/*--------------------------------------------------
  For workflow summary
---------------------------------------------------*/
// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if ( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ) {
    custom_runName = workflow.runName
}

// Summary
def summary = [:]
summary['Pipeline Name']                                   = 'icbi/wes'
summary['Pipeline Version']                                = workflow.manifest.version
if(params.batchFile) summary['Batch file']                 = params.batchFile
if(params.readsNormal != "NO_FILE") summary['Reads normal fastq files'] = params.readsNormal
if(params.readsTumor != "NO_FILE") summary['Reads tumor fastq files']   = params.readsTumor
summary['Read length']                   = params.readLength
summary['Baits bed file']                = params.references.BaitsBed
summary['Regions bed file']              = params.references.RegionsBed
summary['Fasta Ref']                     = params.references.RefFasta
summary['Fasta Index']                   = params.references.RefIdx
summary['Fasta dict ']                   = params.references.RefDict
summary['BWA Index']                     = params.references.BwaRef
summary['VEP Fasta']                     = params.references.VepFasta
summary['MillsGold']                     = params.databases.MillsGold
summary['MillsGoldIdx']                  = params.databases.MillsGoldIdx
summary['hcSNPS1000G']                   = params.databases.hcSNPS1000G
summary['hcSNPS1000GIdx']                = params.databases.hcSNPS1000GIdx
summary['HapMap']                        = params.databases.HapMap
summary['HapMapIdx']                     = params.databases.HapMapIdx
summary['Cosmic']                        = params.databases.Cosmic
summary['CosmicIdx']                     = params.databases.CosmicIdx
summary['DBSNP']                         = params.databases.DBSNP
summary['DBSNPIdx']                      = params.databases.DBSNPIdx
summary['GnomAD']                        = params.databases.GnomAD
summary['GnomADIdx']                     = params.databases.GnomADIdx
summary['GnomADfull']                    = params.databases.GnomADfull
summary['GnomADfullIdx']                 = params.databases.GnomADfullIdx
summary['KnownIndels']                   = params.databases.KnownIndels
summary['KnownIndelsIdx']                = params.databases.KnownIndelsIdx
summary['priority variant Caller']       = params.priorityCaller
summary['Mutect 1 and 2 minAD']          = params.minAD
summary['VarScan min_cov']               = params.min_cov
summary['VarScan min_cov_tumor']         = params.min_cov_tumor
summary['VarScan min_cov_normal']        = params.min_cov_normal
summary['VarScan min_freq_for_hom']      = params.min_freq_for_hom
summary['VarScan tumor_purity']          = params.tumor_purity
summary['VarScan somatic_pvalue']        = params.somatic_pvalue
summary['VarScan somatic_somaticpvalue'] = params.somatic_somaticpvalue
summary['VarScan strand_filter']         = params.strand_filter
summary['VarScan processSomatic_pvalue'] = params.processSomatic_pvalue
summary['VarScan max_normal_freq']       = params.max_normal_freq
summary['VarScan min_tumor_freq']        = params.min_tumor_freq
summary['VarScan min_map_cov']           = params.min_map_cov
summary['VarScan min_base_q']            = params.min_base_q
summary['VEP assembly']                  = params.vep_assembly
summary['VEP species']                   = params.vep_species
summary['VEP options']                   = params.vep_options
summary['Number of scatters']            = params.scatter_count
summary['Max Memory']                    = params.max_memory
summary['Max CPUs']                      = params.max_cpus
summary['Max Time']                      = params.max_time
summary['Output dir']                    = params.outputDir
summary['Working dir']                   = workflow.workDir
summary['TMP dir']                       = params.tmpDir
summary['Current home']                  = "$HOME"
summary['Current user']                  = "$USER"
summary['Current path']                  = "$PWD"
summary['JAVA_Xmx']                      = params.JAVA_Xmx
summary['Picard maxRecordsInRam']        = params.maxRecordsInRam
summary['Script dir']                    = workflow.projectDir
summary['Config Profile']                = workflow.profile


if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "-------------------------------------------------------------------------"


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'icbi-wes-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'icbi/wes Workflow Summary'
    section_href: 'https://gitlab.i-med.ac.at/icbi-lab/pipelines/wes-nf'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}
// End Summary

// did not get a CSV batch file: just run a single sample
if (! params.batchFile) {
        // create channel with tumorSampleName/reads file set
        if (params.readsTumor != "NO_FILE") {
            // Sample name and reads file is passed via cmd line options
            // Sample name to use, if not given uses fastq file simpleName
            params.tumorSampleName  = "undefined"
            params.normalSampleName = "undefined"
            if(!params.single_end) {
                tumorSampleName = params.tumorSampleName != "undefined" ? params.tumorSampleName : file(params.readsTumor)[0].simpleName
            } else {
                tumorSampleName = params.tumorSampleName != "undefined" ? params.tumorSampleName : file(params.readsTumor).simpleName
            }

            Channel
                   .fromFilePairs(params.readsTumor)
                   .map { reads -> tuple(tumorSampleName, reads[1][0], reads[1][1], "None") }
                   .into { raw_reads_tumor_ch;
                           fastqc_reads_tumor_ch }
            Channel
                   .fromFilePairs(params.readsRNAseq)
                   .map { reads -> tuple(tumorSampleName, reads[1][0], reads[1][1], "None") }
                   .into { raw_reads_tumor_neofuse_ch; fastqc_readsRNAseq_ch }
            
        } else  {
            exit 1, "No tumor sample defined"
        }

        if (params.readsNormal != "NO_FILE") {
            if(!params.single_end) {
                normalSampleName = params.normalSampleName != "undefined" ? params.normalSampleName : file(params.readsNormal)[0].simpleName
            } else {
                normalSampleName = params.normalSampleName != "undefined" ? params.normalSampleName : file(params.readsNormal).simpleName
            }
            Channel
                    .fromFilePairs(params.readsNormal)
                    .map { reads -> tuple(normalSampleName, tumorSampleName, reads[1][0], reads[1][1], "None") }
                    .into { raw_reads_normal_ch; fastqc_reads_normal_ch }
        }  else  {
            exit 1, "No normal sample defined"
        }

} else {
        // batchfile ()= csv with sampleId and T reads [N reads] [and group]) was provided
        // create channel with all sampleId/reads file sets from the batch file
        
        // check if reverse reads are specified, if not set up single end processing
        // check if Normal reads are specified, if not set up exit with error
        // attention: if one of the samples is no-Normal or single end, all others
        // will be handled as such. You might want to process mixed sample data as
        // separate batches.
        batchCSV = file(params.batchFile).splitCsv(header:true)
        for ( row in batchCSV ) {
            if (row.readsTumorREV == "None") {
                single_end = true
            }
            if (row.readsNormalFWD == "None") {
                exit 1, "No normal sample defined for " + row.readsTumorFWD
            }
            if (row.readsRNAseqREV == "None") {
                single_end_RNA = true
            }          
        }

        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.tumorSampleName,
                                    row.normalSampleName,
                                    file(row.readsTumorFWD),
                                    file(row.readsTumorREV),
                                    row.group) }
                .into { raw_reads_tumor_ch;
                        fastqc_reads_tumor_ch }
        
        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.tumorSampleName,
                                    row.normalSampleName,
                                    file(row.readsNormalFWD),
                                    file(row.readsNormalREV),
                                    row.group) }
                .into { raw_reads_normal_ch; fastqc_reads_normal_ch }

        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.tumorSampleName,
                                    row.normalSampleName,
                                    file(row.readsRNAseqFWD),
                                    file(row.readsRNAseqREV)) }
                .into { raw_reads_tumor_neofuse_ch; fastqc_readsRNAseq_ch }

}


reference = defineReference()
database = defineDatabases()
mkTmpDir()

scatter_count = Channel.from(params.scatter_count)
padding = params.readLength + 100

FASTQC        = file(params.FASTQC)
FASTP        = file(params.FASTP)
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
BGZIP         = file(params.BGZIP)
TABIX         = file(params.TABIX)
BCFTOOLS      = file(params.BCFTOOLS)
YARA		  = file(params.YARA)
PYTHON		  = file(params.PYTHON)
OPTITYPE	  = file(params.OPTITYPE)
HLAHD		  = file(params.HLAHD)
MIXCR		  = file(params.MIXCR)

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
if (params.WES) {
    process 'RegionsBedToIntervalList' {

        tag 'RegionsBedToIntervalList'

        publishDir "$params.outputDir/00_prepare_Intervals/", mode: params.publishDirMode

        input:
        set(
            file(RefDict),
            file(RegionsBed)
        ) from Channel.value(
            [ reference.RefDict,
            reference.RegionsBed ]
        )

        output:
        file(
            "${RegionsBed.baseName}.interval_list"
        ) into (
            RegionsBedToIntervalList_out_ch0,
            RegionsBedToIntervalList_out_ch1,
            RegionsBedToIntervalList_out_ch2
        )

        script:
        """
        $JAVA8 ${params.JAVA_Xmx} -jar $PICARD BedToIntervalList \\
            I=${RegionsBed} \\
            O=${RegionsBed.baseName}.interval_list \\
            SD=$RefDict
        """
    }

    process 'BaitsBedToIntervalList' {

        tag 'BaitsBedToIntervalList'

        publishDir "$params.outputDir/00_prepare_Intervals/", mode: params.publishDirMode

        input:
        set(
            file(RefDict),
            file(BaitsBed)
        ) from Channel.value(
            [ reference.RefDict,
            reference.BaitsBed ]
        )

        output:
        file(
            "${BaitsBed.baseName}.interval_list"
        ) into (
            BaitsBedToIntervalList_out_ch0,
            BaitsBedToIntervalList_out_ch1
        )

        script:
        """
        $JAVA8 ${params.JAVA_Xmx} -jar $PICARD BedToIntervalList \\
            I=${BaitsBed} \\
            O=${BaitsBed.baseName}.interval_list \\
            SD=$RefDict
        """
    }
} else {
    RegionsBedToIntervalList_out_ch0 = Channel.empty()
    RegionsBedToIntervalList_out_ch1 = Channel.empty()
    RegionsBedToIntervalList_out_ch2 = Channel.empty()
    BaitsBedToIntervalList_out_ch0 = Channel.empty()
    BaitsBedToIntervalList_out_ch1 = Channel.empty()
}

process 'preprocessIntervalList' {

    tag 'preprocessIntervalList'

    publishDir "$params.outputDir/00_prepare_Intervals/", mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )
    file(interval_list) from RegionsBedToIntervalList_out_ch0

    output:
    file(
        "${interval_list.baseName}_merged_padded.interval_list"
    ) into (
        preprocessIntervalList_out_ch0,
        preprocessIntervalList_out_ch1,
        preprocessIntervalList_out_ch2,
        preprocessIntervalList_out_ch3,
        preprocessIntervalList_out_ch4,
        preprocessIntervalList_out_ch5,
        preprocessIntervalList_out_ch6,
        preprocessIntervalList_out_ch7,
        preprocessIntervalList_out_ch8
    )

    script:
    if(params.WES)
        """
        $GATK4 PreprocessIntervals \\
            -R $RefFasta \\
            -L ${interval_list} \\
            --bin-length 0 \\
            --padding ${padding} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            -O ${interval_list.baseName}_merged_padded.interval_list
        """
    else
        """
        $GATK4 PreprocessIntervals \\
            -R $RefFasta \\
            --bin-length 1000 \\
            --padding 0 \\
            -O ${interval_list.baseName}_merged_padded.interval_list
        """
}

process 'SplitIntervals' {
// Splitting interval file in 20(default) files for scattering Mutect2

    tag "SplitIntervals"

    publishDir "$params.outputDir/00_prepare_Intervals/SplitIntervals/", mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch0

    val x from scatter_count

    output:
    file(
        "${IntervalName}/*-scattered.interval_list"
    ) into (
        SplitIntervals_out_ch0,
        SplitIntervals_out_ch1,
        SplitIntervals_out_ch2,
        SplitIntervals_out_ch3,
        SplitIntervals_out_ch4,
        SplitIntervals_out_ch5,
        SplitIntervals_out_ch6,
    )
    val("${IntervalName}") into SplitIntervals_out_ch0_Name

    script:
    IntervalName = IntervalsList.baseName
    """
    mkdir -p ${params.tmpDir}

    $GATK4 SplitIntervals \\
        --tmp-dir ${params.tmpDir} \\
        -R ${RefFasta}  \\
        -scatter ${x} \\
        --interval-merging-rule ALL \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
        -L ${IntervalsList} \\
        -O ${IntervalName}

    """
}

process 'IntervalListToBed' {
    // convert padded interval list to Bed file (used by varscan)
    // generate a padded tabix indexed region BED file for strelka

    tag 'BedFromIntervalList'

    publishDir "$params.outputDir/00_prepare_Intervals/", mode: params.publishDirMode

    input:
        file(paddedIntervalList) from preprocessIntervalList_out_ch1

    output:
    tuple(
        file("${paddedIntervalList.baseName}.bed.gz"),
        file("${paddedIntervalList.baseName}.bed.gz.tbi")
    ) into (
        RegionsBedToTabix_out_ch0,
        RegionsBedToTabix_out_ch1
    )

    script:
    """
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD IntervalListToBed \\
        I=${paddedIntervalList} \\
        O=${paddedIntervalList.baseName}.bed
    
    ${BGZIP} -c ${paddedIntervalList.baseName}.bed > ${paddedIntervalList.baseName}.bed.gz &&
    ${TABIX} -p bed ${paddedIntervalList.baseName}.bed.gz
    """
}

process 'ScatteredIntervalListToBed' {
    // convert scattered padded interval list to Bed file (used by varscan)

    tag 'ScatteredIntervalListToBed'

    publishDir "$params.outputDir/00_prepare_Intervals/SplitIntervals/${IntervalName}", mode: params.publishDirMode

    input:
    set(
        val(IntervalName),
        file(IntervalsList)
    ) from SplitIntervals_out_ch0_Name
        .combine(
            SplitIntervals_out_ch0.flatten()
        )


    output:
    file(
        "*.bed"
    ) into (
        ScatteredIntervalListToBed_out_ch0
    )

    script:
    """
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD IntervalListToBed \\
        I=${IntervalsList} \\
        O=${IntervalsList.baseName}.bed
    """
}

// FastQC
process FastQC {
    tag "$TumorReplicateId"

    publishDir "${params.outputDir}/$TumorReplicateId/02_QC/",
        mode: params.publishDirMode,
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    if (!single_end && !single_end_RNA) {
        cpus = 6
    } else if (!single_end && single_end_RNA) {
        cpus = 5
    } else if (single_end && !single_end_RNA) {
        cpus = 4
    } else if (single_end && single_end_RNA) {
        cpus = 3
    }


    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(tumor_readsFWD),
        file(tumor_readsREV),
        sampleGroup,      // unused so far
    ) from fastqc_reads_tumor_ch

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(normal_readsFWD),
        file(normal_readsREV),
        sampleGroup,      // unused so far
    ) from fastqc_reads_normal_ch

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(readsRNAseq_FWD),
        file(readsRNAseq_REV)
    ) from fastqc_readsRNAseq_ch

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("*_fastqc*")
    ) into ch_fastqc // multiQC

    script:
    tumor_readsFWD_simpleName = tumor_readsFWD.getSimpleName()
    normal_readsFWD_simpleName = normal_readsFWD.getSimpleName()
    tumor_readsFWD_ext = tumor_readsFWD.getExtension()
    normal_readsFWD_ext = normal_readsFWD.getExtension()

    readsRNAseq_FWD_simpleName = readsRNAseq_FWD.getSimpleName()
    readsRNAseq_FWD_ext = readsRNAseq_FWD.getExtension()

    tumor_readsFWD_ext = (tumor_readsFWD_ext == "gz") ? "fastq.gz" : tumor_readsFWD_ext
    normal_readsFWD_ext = (normal_readsFWD_ext == "gz") ? "fastq.gz" : normal_readsFWD_ext

    readsRNAseq_FWD_ext = (readsRNAseq_FWD_ext == "gz") ? "fastq.gz" : readsRNAseq_FWD_ext


    if (! single_end) {
        tumor_readsREV_simpleName = tumor_readsREV.getSimpleName()
        normal_readsREV_simpleName = normal_readsREV.getSimpleName()
        tumor_readsREV_ext = tumor_readsREV.getExtension()
        normal_readsREV_ext = normal_readsREV.getExtension()

        tumor_readsREV_ext = (tumor_readsREV_ext == "gz") ? "fastq.gz" : tumor_readsREV_ext
        normal_readsREV_ext = (normal_readsREV_ext == "gz") ? "fastq.gz" : normal_readsREV_ext
    }
    if (! single_end_RNA) {
        readsRNAseq_REV_simpleName = readsRNAseq_REV.getSimpleName()
        readsRNAseq_REV_ext = readsRNAseq_REV.getExtension()

        readsRNAseq_REV_ext = (readsRNAseq_REV_ext == "gz") ? "fastq.gz" : readsRNAseq_REV_ext
    }

    if (single_end && single_end_RNA)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}_R1.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}_R1.${normal_readsFWD_ext}
        ln -s $readsRNAseq_FWD ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}

        fastqc --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext} \\
            ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}
        """
    else if (!single_end && !single_end_RNA)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}_R1.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}_R1.${normal_readsFWD_ext}
        ln -s $tumor_readsREV ${TumorReplicateId}_R2.${tumor_readsREV_ext}
        ln -s $normal_readsREV ${NormalReplicateId}_R2.${normal_readsREV_ext}

        ln -s $readsRNAseq_FWD ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}
        ln -s $readsRNAseq_REV ${TumorReplicateId}_RNA_R2.${readsRNAseq_REV_ext}

        $FASTQC --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} ${TumorReplicateId}_R2.${tumor_readsREV_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext} ${NormalReplicateId}_R2.${normal_readsREV_ext} \\
            ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext} ${TumorReplicateId}_RNA_R2.${readsRNAseq_REV_ext}
        """
    else if (!single_end && single_end_RNA)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}_R1.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}_R1.${normal_readsFWD_ext}
        ln -s $tumor_readsREV ${TumorReplicateId}_R2.${tumor_readsREV_ext}
        ln -s $normal_readsREV ${NormalReplicateId}_R2.${normal_readsREV_ext}

        ln -s $readsRNAseq_FWD ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}

        $FASTQC --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} ${TumorReplicateId}_R2.${tumor_readsREV_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext} ${NormalReplicateId}_R2.${normal_readsREV_ext} \\
            ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}
        """
    else if (single_end && !single_end_RNA)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}.${normal_readsFWD_ext}

        ln -s $readsRNAseq_FWD ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}
        ln -s $readsRNAseq_REV ${TumorReplicateId}_RNA_R2.${readsRNAseq_REV_ext}

        $FASTQC --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext} \\
            ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext} ${TumorReplicateId}_RNA_R2.${readsRNAseq_REV_ext}
        """

}

// adapter trimming Tumor
if (params.trim_adapters) {
    process fastp_tumor {

        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(tumor_readsFWD),
            file(tumor_readsREV),
            sampleGroup,      // unused so far
        ) from raw_reads_tumor_ch

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("${TumorReplicateId}_trimmed_R1.fastq.gz"),
            file("${trimmedReads_2}"),
            sampleGroup
        ) into (
            reads_tumor_ch,
            reads_tumor_mixcr_DNA_ch,
            fastqc_reads_tumor_trimmed_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*.json")
        ) into ch_fastp_tumor // multiQC


        script:
        trimmedReads_2 = (single_end) ? val("NO_FILE") : TumorReplicateId + "_trimmed_R2.fastq.gz"

        def fastpAdapter = ''
        if(params.adapterSeqFile != false) {
            val adapterSeqFile = Channel.fromPath(params.adapterSeqFile)
            fastpAdapter = "--adapter_fasta $adapterSeqFile"
        } else {
            if(params.adapterSeq != false) {
                adapterSeq   = Channel.value(params.adapterSeq)
                fastpAdapter = "--adapter_sequence " + adapterSeq.getVal()
            
                if(params.adapterSeqR2 != false) {
                    adapterSeqR2   = Channel.value(params.adapterSeqR2)
                    fastpAdapter += " --adapter_sequence_r2 " + adapterSeqR2.getVal()
                }
            }
        }

        if(single_end)
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${tumor_readsFWD} \\
                --out1 ${TumorReplicateId}_trimmed_R1.fastq.gz \\
                --json ${TumorReplicateId}_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            """
        else
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${tumor_readsFWD} \\
                --in2 ${tumor_readsREV} \\
                --out1 ${TumorReplicateId}_trimmed_R1.fastq.gz \\
                --out2 ${TumorReplicateId}_trimmed_R2.fastq.gz \\
                --json ${TumorReplicateId}_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            """
    }

    process fastp_normal {

        tag "$NormalReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(normal_readsFWD),
            file(normal_readsREV),
            sampleGroup,      // unused so far
        ) from raw_reads_normal_ch

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("${NormalReplicateId}_trimmed_R1.fastq.gz"),
            file("${trimmedReads_2}"),
            sampleGroup
        ) into (
            reads_normal_ch,
            reads_normal_mixcr_DNA_ch,
            fastqc_reads_normal_trimmed_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*.json")
        ) into ch_fastp_normal // multiQC


        script:
        trimmedReads_2 = (single_end) ? val("NO_FILE") : NormalReplicateId + "_trimmed_R2.fastq.gz"

        def fastpAdapter = ''
        if(params.adapterSeqFile != false) {
            val adapterSeqFile = Channel.fromPath(params.adapterSeqFile)
            fastpAdapter = "--adapter_fasta $adapterSeqFile"
        } else {
            if(params.adapterSeq != false) {
                adapterSeq   = Channel.value(params.adapterSeq)
                fastpAdapter = "--adapter_sequence " + adapterSeq.getVal()
            
                if(params.adapterSeqR2 != false) {
                    adapterSeqR2   = Channel.value(params.adapterSeqR2)
                    fastpAdapter += " --adapter_sequence_r2 " + adapterSeqR2.getVal()
                }
            }
        }

        if(single_end)
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${normal_readsFWD} \\
                --out1 ${NormalReplicateId}_trimmed_R1.fastq.gz \\
                --json ${NormalReplicateId}_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            """
        else
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${normal_readsFWD} \\
                --in2 ${normal_readsREV} \\
                --out1 ${NormalReplicateId}_trimmed_R1.fastq.gz \\
                --out2 ${NormalReplicateId}_trimmed_R2.fastq.gz \\
                --json ${NormalReplicateId}_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            """
    }


    // FastQC after adapter trimming
    process FastQC_trimmed {
        tag "$TumorReplicateId"

        publishDir "${params.outputDir}/$TumorReplicateId/02_QC/",
            mode: params.publishDirMode,
            saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

        if (single_end) {
            cpus = 2
        } else {
            cpus = 4
        }

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(tumor_readsFWD),
            file(tumor_readsREV),
            sampleGroup,      // unused so far
        ) from fastqc_reads_tumor_trimmed_ch

        set(
            TumorReplicateId,
            NormalReplicateId,
            file(normal_readsFWD),
            file(normal_readsREV),
            sampleGroup,      // unused so far
        ) from fastqc_reads_normal_trimmed_ch


        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*_fastqc*")
        ) into ch_fastqc_trimmed // multiQC

        script:
        if (single_end)
            """
            fastqc --quiet --threads ${task.cpus} \\
                ${tumor_readsFWD} \\
                ${normal_readsFWD}
            """
        else
            """
            $FASTQC --quiet --threads ${task.cpus} \\
                ${tumor_readsFWD} ${tumor_readsREV} \\
                ${normal_readsFWD} ${normal_readsREV}
            """
    }
} else { // no adapter trimming
    process mk_readChannelsFromRaw {
        tag "$TumorReplicateId"

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(tumor_readsFWD),
            file(tumor_readsREV),
            sampleGroup,      // unused so far
        ) from raw_reads_tumor_ch

        set(
            TumorReplicateId,
            NormalReplicateId,
            file(normal_readsFWD),
            file(normal_readsREV),
            sampleGroup,      // unused so far
        ) from raw_reads_normal_ch

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(tumor_readsFWD),
            file(tumor_readsREV),
            sampleGroup,      // unused so far
        ) into (
            reads_tumor_ch,
            reads_tumor_mixcr_DNA_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(normal_readsFWD),
            file(normal_readsREV),
            sampleGroup,      // unused so far
        ) into (
            reads_normal_ch,
            reads_normal_mixcr_DNA_ch
        )

        set(
            TumorReplicateId,
            NormalReplicateId,
            ""
        ) into (
            ch_fastqc_trimmed,
            ch_fastp_tumor,
            ch_fastp_normal
        )

        // do nothing
        """
        """
    }    
}

// adapter trimming RNAseq
if (params.trim_adapters_RNAseq) {
    process fastp_RNAseq {
        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(readsRNAseq_FWD),
            file(readsRNAseq_REV),
        ) from raw_reads_tumor_neofuse_ch

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("${TumorReplicateId}_RNA_trimmed_R1.fastq.gz"),
            file("${trimmedReads_2}")
        ) into (
            reads_tumor_neofuse_ch,
            reads_tumor_mixcr_RNA_ch,
            fastqc_readsRNAseq_trimmed_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*.json")
        ) into ch_fastp_RNAseq // multiQC

        script:
        trimmedReads_2 = (single_end_RNA) ? val("NO_FILE") : TumorReplicateId + "_RNA_trimmed_R2.fastq.gz"


        def fastpAdapter = ''
        if(params.adapterSeqFileRNAseq != false) {
            val adapterSeqFile = Channel.fromPath(params.adapterSeqFileRNAseq)
            fastpAdapter = "--adapter_fasta $adapterSeqFile"
        } else {
            if(params.adapterSeqRNAseq != false) {
                adapterSeq   = Channel.value(params.adapterSeqRNAseq)
                fastpAdapter = "--adapter_sequence " + adapterSeq.getVal()
            
                if(params.adapterSeqR2RNAseq != false) {
                    adapterSeqR2   = Channel.value(params.adapterSeqR2RNAseq)
                    fastpAdapter += " --adapter_sequence_r2 " + adapterSeqR2.getVal()
                }
            }
        }

        if(single_end)
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${readsRNAseq_FWD} \\
                --out1 ${TumorReplicateId}_RNA_trimmed_R1.fastq.gz \\
                --json ${TumorReplicateId}_RNA_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            """
        else
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${readsRNAseq_FWD} \\
                --in2 ${readsRNAseq_REV} \\
                --out1 ${TumorReplicateId}_RNA_trimmed_R1.fastq.gz \\
                --out2 ${TumorReplicateId}_RNA_trimmed_R2.fastq.gz \\
                --json ${TumorReplicateId}_RNA_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            """
    }

    // FastQC after RNAseq adapter trimming
    process FastQC_trimmed_RNAseq {
        tag "$TumorReplicateId"

        publishDir "${params.outputDir}/$TumorReplicateId/02_QC/",
            mode: params.publishDirMode,
            saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

        if (single_end_RNA) {
            cpus = 1
        } else {
            cpus = 2
        }

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(readsRNAseq_FWD),
            file(readsRNAseq_REV)
        ) from fastqc_readsRNAseq_trimmed_ch

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*_fastqc*")
        ) into ch_fastqc_trimmed_RNAseq // multiQC

        script:
        if (single_end_RNA)
            """
            fastqc --quiet --threads ${task.cpus} ${tumor_readsFWD}
            """
        else
            """
            $FASTQC --quiet --threads ${task.cpus} \\
                ${readsRNAseq_FWD} ${readsRNAseq_REV}
            """
    }

} else { // no adapter trimming for RNAseq
    process mk_readRNAseqChannelsFromRaw {
        tag "$TumorReplicateId"

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(readsRNAseq_FWD),
            file(readsRNAseq_REV),
        ) from raw_reads_tumor_neofuse_ch

        empty = []

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(readsRNAseq_FWD),
            file(readsRNAseq_REV),
        ) into (
            reads_tumor_neofuse_ch,
            reads_tumor_mixcr_RNA_ch
        )

        set(
            TumorReplicateId,
            NormalReplicateId,
            ""
        ) into (
            fastqc_readsRNAseq_trimmed_ch,
            ch_fastp_RNAseq,
            ch_fastqc_trimmed_RNAseq
        ) 

        // do nothing
        """
        """    
    }
}


/// start processing reads
process 'BwaTumor' {
// Aligning tumor reads to reference, sort and index; create BAMs

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(readsFWD),
        file(readsREV),
        sampleGroup,      // unused so far
    ) from reads_tumor_ch
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
        file(BwaRef)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.BwaRef ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_aligned.bam")
    ) into BwaTumor_out_ch0

    script:
    if (single_end)
        """
        $BWA mem \\
            -R "@RG\\tID:${TumorReplicateId}\\tLB:${TumorReplicateId}\\tSM:${TumorReplicateId}\\tPL:ILLUMINA" \\
            -M ${RefFasta} \\
            -t ${task.cpus} \\
            ${readsFWD} | \\
            $SAMTOOLS view -Shb -o ${TumorReplicateId}_aligned.bam -
        """
    else
        """
        $BWA mem \\
            -R "@RG\\tID:${TumorReplicateId}\\tLB:${TumorReplicateId}\\tSM:${TumorReplicateId}\\tPL:ILLUMINA" \\
            -M ${RefFasta} \\
            -t ${task.cpus} \\
            ${readsFWD} \\
            ${readsREV}  | \\
            $SAMTOOLS view -Shb -o ${TumorReplicateId}_aligned.bam -
        """
}


process 'MarkDuplicatesTumor' {
// Mark duplicates with Picard

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam)
    ) from BwaTumor_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_aligned_sort_mkdp.bam"),
        file("${TumorReplicateId}_aligned_sort_mkdp.bam.bai"),
        file("${TumorReplicateId}_aligned_sort_mkdp.txt")
    ) into (
        MarkDuplicatesTumor_out_ch0,
        MarkDuplicatesTumor_out_ch1,
        MarkDuplicatesTumor_out_ch2,
        MarkDuplicatesTumor_out_ch3 // mhc_extract -> hld-hd, optitype
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_aligned_sort_mkdp.txt")
    ) into MarkDuplicatesTumor_out_ch4 // multiQC

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 MarkDuplicatesSpark \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -O ${TumorReplicateId}_aligned_sort_mkdp.bam \\
        -M ${TumorReplicateId}_aligned_sort_mkdp.txt \\
        --create-output-bam-index true \\
        --read-validation-stringency LENIENT \\
        --spark-master local[${task.cpus}] \\
        --conf 'spark.executor.cores=${task.cpus}' \\
        --conf 'spark.local.dir=${params.tmpDir}' 2> /dev/stdout
    """
}

if(params.WES) {
    process 'alignmentMetricsTumor' {
    // Generate HS metrics using picard

        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/02_QC/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(bam),
            file(bai),
            _
        ) from MarkDuplicatesTumor_out_ch0
        
        set(
            file(RefFasta),
            file(RefIdx)
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx ]
        )

        file(BaitIntervalsList) from BaitsBedToIntervalList_out_ch0
        file(IntervalsList) from RegionsBedToIntervalList_out_ch1


        output:
        set(TumorReplicateId,
            NormalReplicateId,
            file("${TumorReplicateId}.*.txt")
        ) into alignmentMetricsTumor_ch // multiQC

        script:
        """
        mkdir -p ${params.tmpDir}
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectHsMetrics \\
            TMP_DIR=${params.tmpDir} \\
            INPUT=${bam} \\
            OUTPUT=${TumorReplicateId}.HS.metrics.txt \\
            R=${RefFasta} \\
            BAIT_INTERVALS=${BaitIntervalsList} \\
            TARGET_INTERVALS=${IntervalsList} \\
            PER_TARGET_COVERAGE=${TumorReplicateId}.perTarget.coverage.txt && \\
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectAlignmentSummaryMetrics \\
            TMP_DIR=${params.tmpDir} \\
            INPUT=${bam} \\
            OUTPUT=${TumorReplicateId}.AS.metrics.txt \\
            R=${RefFasta} &&
        $SAMTOOLS flagstat -@${task.cpus} ${bam} > ${TumorReplicateId}.flagstat.txt
        """
    }
} else {
    // bogus channel for multiqc
    process mk_bogus_alignmentMetricsTumor_ch {
        tag "$TumorReplicateId"

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            _,
            _
        ) from MarkDuplicatesTumor_out_ch0


        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            ""
        ) into alignmentMetricsTumor_ch

        // do nothing
        """
        """    
    }    
}

process 'BwaNormal' {
// Aligning Normal reads to reference, sort and index; create BAMs

    tag "$NormalReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(readsFWD),
        file(readsREV),
        sampleGroup      // unused so far
    ) from reads_normal_ch

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
        file(BwaRef)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.BwaRef ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_aligned.bam")
    ) into BwaNormal_out_ch0

    script:
    if(single_end)
        """
        $BWA mem \\
            -R "@RG\\tID:${NormalReplicateId}\\tLB:${NormalReplicateId}\\tSM:${NormalReplicateId}\\tPL:ILLUMINA" \\
            -M ${RefFasta} \\
            -t ${task.cpus} \\
            ${readsFWD} | \\
        $SAMTOOLS  view -Shb -o ${NormalReplicateId}_aligned.bam -
        """
    else
        """
        $BWA mem \\
            -R "@RG\\tID:${NormalReplicateId}\\tLB:${NormalReplicateId}\\tSM:${NormalReplicateId}\\tPL:ILLUMINA" \\
            -M ${RefFasta} \\
            -t ${task.cpus} \\
            ${readsFWD} \\
            ${readsREV} | \\
        $SAMTOOLS view -Shb -o ${NormalReplicateId}_aligned.bam -
    """
}


process 'MarkDuplicatesNormal' {
// Mark duplicates with Picard

    tag "$NormalReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam)
    ) from BwaNormal_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_aligned_sort_mkdp.bam"),
        file("${NormalReplicateId}_aligned_sort_mkdp.bam.bai"),
        file("${NormalReplicateId}_aligned_sort_mkdp.txt")
    ) into (
        MarkDuplicatesNormal_out_ch0,
        MarkDuplicatesNormal_out_ch1,
        MarkDuplicatesNormal_out_ch2
    )
    
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_aligned_sort_mkdp.txt")
    ) into MarkDuplicatesNormal_out_ch3 // multiQC

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 MarkDuplicatesSpark \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -O ${NormalReplicateId}_aligned_sort_mkdp.bam \\
        -M ${NormalReplicateId}_aligned_sort_mkdp.txt \\
        --create-output-bam-index true \\
        --read-validation-stringency LENIENT \\
        --spark-master local[${task.cpus}] \\
        --conf 'spark.executor.cores=${task.cpus}' \\
        --conf 'spark.local.dir=${params.tmpDir}' 2> /dev/stdout
    """
}

if (params.WES) {
    process 'alignmentMetricsNormal' {
    // Generate HS metrics using picard

        tag "$NormalReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/02_QC/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(bam),
            file(bai),
            _
        ) from MarkDuplicatesNormal_out_ch0

        set(
            file(RefFasta),
            file(RefIdx)
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx ]
        )

        file(BaitIntervalsList) from BaitsBedToIntervalList_out_ch1
        file(IntervalsList) from RegionsBedToIntervalList_out_ch2


        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("${NormalReplicateId}.*.txt")
        ) into alignmentMetricsNormal_ch // multiQC


        script:
        """
        mkdir -p ${params.tmpDir}
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectHsMetrics \\
            TMP_DIR=${params.tmpDir} \\
            INPUT=${bam} \\
            OUTPUT=${NormalReplicateId}.HS.metrics.txt \\
            R=${RefFasta} \\
            BAIT_INTERVALS=${BaitIntervalsList} \\
            TARGET_INTERVALS=${IntervalsList} \\
            PER_TARGET_COVERAGE=${NormalReplicateId}.perTarget.coverage.txt && \\
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectAlignmentSummaryMetrics \\
            TMP_DIR=${params.tmpDir} \\
            INPUT=${bam} \\
            OUTPUT=${NormalReplicateId}.AS.metrics.txt \\
            R=${RefFasta} && \\
        $SAMTOOLS flagstat -@${task.cpus} ${bam} > ${NormalReplicateId}.flagstat.txt
        """
    }
} else {
    // bogus channel for multiqc
    process mk_bogus_alignmentMetricsNormal_ch {
        tag "$TumorReplicateId"

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            _,
            _
        ) from MarkDuplicatesTumor_out_ch0


        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            ""
        ) into alignmentMetricsTumor_ch
        
        // do nothing
        """
        """    
    }    
}


/*
*********************************************
**             M U T E C T 2               **
*********************************************
*/

process 'BaseRecalTumorGATK4' {
/*
 BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
 Recalibration (BQSR)
 ApplyBQSR (GATK4): apply BQSR table to reads
*/

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai),
        file(list)
    ) from MarkDuplicatesTumor_out_ch1

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch2

    set(
        file(MillsGold),
        file(MillsGoldIdx),
        file(DBSNP),
        file(DBSNPIdx),
        file(KnownIndels),
        file(KnownIndelsIdx)
    ) from Channel.value(
        [ database.MillsGold,
          database.MillsGoldIdx,
          database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_bqsr.table")
    ) into BaseRecalTumorGATK4_out_ch0

    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_recal4.bam"),
        file("${TumorReplicateId}_recal4.bam.bai")
    ) into (
        BaseRecalTumorGATK4_out_ch1,
        BaseRecalTumorGATK4_out_ch2,
        BaseRecalTumorGATK4_out_ch3,
        BaseRecalTumorGATK4_out_ch4,
        BaseRecalTumorGATK4_out_ch5
    )


    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD SetNmMdAndUqTags \\
        TMP_DIR=${params.tmpDir} \\
        R=${RefFasta} \\
        I=${bam} \\
        O=fixed.bam \\
        CREATE_INDEX=true \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam} \\
        VALIDATION_STRINGENCY=LENIENT && \\
    # TODO: re-add later when more stable. Commented out do to crashing randomly
    # $GATK4 BaseRecalibratorSpark \\
    #    --java-options '${params.JAVA_Xmx_spark}' \\
    #    --tmp-dir ${params.tmpDir} \\
    #    -I fixed.bam \\
    #    -R ${RefFasta} \\
    #    -L ${IntervalsList} \\
    #    -O ${TumorReplicateId}_bqsr.table \\
    #    --known-sites ${DBSNP} \\
    #    --known-sites ${KnownIndels} \\
    #    --known-sites ${MillsGold} \\
    #    --spark-master local[${task.cpus}] \\
    #    --conf 'spark.executor.cores=${task.cpus}' \\
    #    --conf 'spark.local.dir=${params.tmpDir}' && \\
    $GATK4 BaseRecalibrator \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I fixed.bam \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${TumorReplicateId}_bqsr.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold} && \\
    $GATK4 ApplyBQSRSpark \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${TumorReplicateId}_recal4.bam \\
        --bqsr-recal-file ${TumorReplicateId}_bqsr.table \\
        --spark-master local[${task.cpus}] \\
        --conf 'spark.executor.cores=${task.cpus}' \\
        --conf 'spark.local.dir=${params.tmpDir}'
    """
}

process 'GetPileupTumor' {
// GetPileupSummaries (GATK4): tabulates pileup metrics for inferring contamination

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
        mode: params.publishDirMode

    input:
    set(
        file(GnomAD),
        file(GnomADIdx)
    ) from Channel.value(
        [ database.GnomAD,
          database.GnomADIdx ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch3
    
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from BaseRecalTumorGATK4_out_ch1

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_pileup.table")
    ) into (
        GetPileupTumor_out_ch0,
        GetPileupTumor_out_ch1
    )

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 GetPileupSummaries \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -O ${TumorReplicateId}_pileup.table \\
        -L ${IntervalsList} \\
        --variant ${GnomAD}
    """
}

process 'AnalyzeCovariates' {
/*
 2nd BaseRecalibrator (GATK4)
 AnalyzeCovariates (GATK4): creates plots to visualize base recalibration results
*/

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/02_QC/",
        mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch4

    set(
        file(DBSNP), 
        file(DBSNPIdx),
        file(KnownIndels),
        file(KnownIndelsIdx),
        file(MillsGold),
        file(MillsGoldIdx)
    ) from Channel.value(
        [ database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx,
          database.MillsGold,
          database.MillsGoldIdx ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(recalTable)
    ) from BaseRecalTumorGATK4_out_ch0

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from BaseRecalTumorGATK4_out_ch2

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_postbqsr.table")
    ) into AnalyzeCovariates_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    # TODO: re-add later when more stable. Commented out do to crashing randomly
    # $GATK4 BaseRecalibratorSpark \\
    #    --java-options '${params.JAVA_Xmx_spark}' \\
    #    --tmp-dir ${params.tmpDir} \\
    #    -I ${bam} \\
    #    -R ${RefFasta} \\
    #    -L ${IntervalsList} \\
    #    -O ${TumorReplicateId}_postbqsr.table \\
    #    --known-sites ${DBSNP} \\
    #    --known-sites ${KnownIndels} \\
    #    --known-sites ${MillsGold} \\
    #    --spark-master local[${task.cpus}] \\
    #    --conf 'spark.executor.cores=${task.cpus}' \\
    #    --conf 'spark.local.dir=${params.tmpDir}' && \\
    $GATK4 BaseRecalibrator \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${TumorReplicateId}_postbqsr.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold} && \\
    $GATK4 AnalyzeCovariates \\
        --tmp-dir ${params.tmpDir} \\
        -before ${recalTable} \\
        -after ${TumorReplicateId}_postbqsr.table \\
        -csv ${TumorReplicateId}_BQSR.csv
    """
}

process 'CollectSequencingArtifactMetrics' {
/*
CollectSequencingArtifactMetrics (Picard): collect metrics to quantify single-base
sequencing artifacts.
It has to be done only for paired end
CollectSequencingArtifactMetrics and FilterByOrientationBias has been deprecated
in favor of our new orientation bias workflow that uses LearnReadOrientationModel
see "https://github.com/broadinstitute/gatk/issues/5553"

If single-end reads are used, do nothing, just create an empty file!!!
*/

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/02_QC/",
        mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from BaseRecalTumorGATK4_out_ch3

    output:
    file(
        "${TumorReplicateId}.pre_adapter_detail_metrics"
    ) into (
        CollectSequencingArtifactMetrics_out_ch0
    )

    script:
    if(single_end)
        """
        touch ${TumorReplicateId}.pre_adapter_detail_metrics
        """
    else
        """
        mkdir -p ${params.tmpDir}

        $JAVA8 ${params.JAVA_Xmx} -jar $PICARD CollectSequencingArtifactMetrics \\
            TMP_DIR=${params.tmpDir} \\
            I=${bam} \\
            R=${RefFasta} \\
            O=${TumorReplicateId}
        """
}


process 'BaseRecalNormalGATK4' {
/*
    BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
    Recalibration (BQSR)
    ApplyBQSR (GATK4): apply BQSR table to reads
*/

    tag "$NormalReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam), file(bai),
        file(list)
    ) from MarkDuplicatesNormal_out_ch1

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict
        ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch5

    set(
        file(MillsGold),
        file(MillsGoldIdx),
        file(DBSNP),
        file(DBSNPIdx),
        file(KnownIndels),
        file(KnownIndelsIdx)
    ) from Channel.value(
        [ database.MillsGold,
          database.MillsGoldIdx,
          database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_bqsr.table")
    ) into BaseRecalNormal_out_ch0

    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_recal4.bam"),
        file("${NormalReplicateId}_recal4.bam.bai")
    ) into (
        BaseRecalNormal_out_ch1,
        BaseRecalNormal_out_ch2,
        BaseRecalNormal_out_ch3  // for HaploTypeCaller
    )

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD SetNmMdAndUqTags \\
        TMP_DIR=${params.tmpDir} \\
        R=${RefFasta} \\
        I=${bam} \\
        O=Normal_fixed.bam \\
        CREATE_INDEX=true \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam} \\
        VALIDATION_STRINGENCY=LENIENT && \\
    # TODO: re-add later when more stable. Commented out do to crashing randomly
    # $GATK4 BaseRecalibratorSpark \\
    #    --java-options '${params.JAVA_Xmx_spark}' \\
    #    --tmp-dir ${params.tmpDir} \\
    #    -I Normal_fixed.bam \\
    #    -R ${RefFasta} \\
    #    -L ${IntervalsList} \\
    #    -O ${NormalReplicateId}_bqsr.table \\
    #    --known-sites ${DBSNP} \\
    #    --known-sites ${KnownIndels} \\
    #    --known-sites ${MillsGold} \\
    #    --spark-master local[${task.cpus}] \\
    #    --conf 'spark.executor.cores=${task.cpus}' \\
    #    --conf 'spark.local.dir=${params.tmpDir}' && \\
    $GATK4 BaseRecalibrator \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I Normal_fixed.bam \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${NormalReplicateId}_bqsr.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold} && \\
    $GATK4 ApplyBQSRSpark \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${NormalReplicateId}_recal4.bam \\
        --bqsr-recal-file ${NormalReplicateId}_bqsr.table \\
        --spark-master local[${task.cpus}] \\
        --conf 'spark.executor.cores=${task.cpus}' \\
        --conf 'spark.local.dir=${params.tmpDir}'
    """
}

process 'GetPileupNormal' {
// GetPileupSummaries (GATK4): tabulates pileup metrics for inferring contamination

    tag "$NormalReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
        mode: params.publishDirMode

    input:
    set(
        file(GnomAD),
        file(GnomADIdx)
    ) from Channel.value(
        [ database.GnomAD,
          database.GnomADIdx ]
    )

    file(intervals) from preprocessIntervalList_out_ch6

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from BaseRecalNormal_out_ch1

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_pileup.table")
    ) into GetPileupNormal_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 GetPileupSummaries \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -O ${NormalReplicateId}_pileup.table \\
        -L ${intervals} \\
        --variant ${GnomAD}
    """
}

process 'Mutect2' {
/* 
    Call somatic SNPs and indels via local re-assembly of haplotypes; tumor sample
    and matched normal sample 
*/

    tag "$TumorReplicateId"

    // publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
    //     mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
        file(gnomADfull),
        file(gnomADfullIdx)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          database.GnomADfull,
          database.GnomADfullIdx ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Tumorbam),
        file(Tumorbai),
        _,                    // unused TumorReplicateId from BaseRecalNormal_out_ch2
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from BaseRecalTumorGATK4_out_ch5
        .combine(BaseRecalNormal_out_ch2, by: 0)
        .combine(
            SplitIntervals_out_ch2.flatten()
        )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${intervals}.vcf.gz"),
        file("${TumorReplicateId}_${intervals}.vcf.gz.stats"),
        file("${TumorReplicateId}_${intervals}.vcf.gz.tbi")
    ) into Mutect2_out_ch0


    script:
    if(params.mutect2ponFile != false) {
        val mutect2ponFile = Channel.fromPath(params.mutect2ponFile)
        panel_of_normals = "--panel-of-normals ${mutect2ponFile}"
    } else {
        panel_of_normals = ""
    }
    """
    mkdir -p ${params.tmpDir}

    $GATK4 Mutect2 \\
        --tmp-dir ${params.tmpDir} \\
        -R ${RefFasta} \\
        -I ${Tumorbam} -tumor ${TumorReplicateId} \\
        -I ${Normalbam} -normal ${NormalReplicateId} \\
        --germline-resource ${gnomADfull} \\
        ${panel_of_normals} \\
        -L ${intervals} \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O ${TumorReplicateId}_${intervals}.vcf.gz
    """
}

process 'gatherMutect2VCFs' {
// Merge scattered Mutect2 vcfs

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(vcf),
        file(stats),
        file(idx)
    ) from Mutect2_out_ch0
        .groupTuple(by: [0, 1])


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz.tbi"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz.stats")
    ) into gatherMutect2VCFs_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}
    
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${vcf.join(" I=")} \\
        O=${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz

    $GATK4 MergeMutectStats \\
        --tmp-dir ${params.tmpDir} \\
        --stats ${stats.join(" --stats ")} \\
        -O ${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz.stats
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
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(pileupTumor),
        _,
        file(pileupNormal),
        _,
        file(vcf),
        file(vcfIdx),
        file(vcfStats)
    ) from GetPileupTumor_out_ch1
    .combine(GetPileupNormal_out_ch0, by :0)
    .combine(gatherMutect2VCFs_out_ch0, by :0)

    file(preAdapterDetail) from CollectSequencingArtifactMetrics_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        val("mutect2"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_final.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_final.vcf.gz.tbi")
    ) into (
        FilterMutect2_out_ch0,
        FilterMutect2_out_ch1
    )

    script:
    if (single_end)
        """
        mkdir -p ${params.tmpDir}

        $GATK4 CalculateContamination \\
            --tmp-dir ${params.tmpDir} \\
            -I ${pileupTumor} \\
            --matched-normal ${pileupNormal} \\
            -O ${TumorReplicateId}_${NormalReplicateId}_cont.table && \\
        $GATK4 FilterMutectCalls \\
            --tmp-dir ${params.tmpDir} \\
            -R ${RefFasta} \\
            -V ${vcf} \\
            --contamination-table ${TumorReplicateId}_${NormalReplicateId}_cont.table \\
            -O ${TumorReplicateId}_${NormalReplicateId}_oncefiltered.vcf.gz && \\
        $GATK4 SelectVariants \\
            --tmp-dir ${params.tmpDir} \\
            --variant ${TumorReplicateId}_${NormalReplicateId}_oncefiltered.vcf.gz \\
            -R ${RefFasta} \\
            --exclude-filtered true \\
            --select 'vc.getGenotype(\"${TumorReplicateId}\").getAD().1 >= ${params.minAD}' \\
            --output ${TumorReplicateId}_${NormalReplicateId}_mutect2_final.vcf
        """
    else
        """
        mkdir -p ${params.tmpDir}

        $GATK4 CalculateContamination \\
            --tmp-dir ${params.tmpDir} \\
            -I ${pileupTumor} \\
            --matched-normal ${pileupNormal} \\
            -O ${TumorReplicateId}_${NormalReplicateId}_cont.table && \\
        $GATK4 FilterMutectCalls \\
            --tmp-dir ${params.tmpDir} \\
            -R ${RefFasta} \\
            -V ${vcf} \\
            --contamination-table ${TumorReplicateId}_${NormalReplicateId}_cont.table \\
            -O ${TumorReplicateId}_${NormalReplicateId}_oncefiltered.vcf.gz && \\
        $GATK4 FilterByOrientationBias \\
            --tmp-dir ${params.tmpDir} \\
            -V ${TumorReplicateId}_${NormalReplicateId}_oncefiltered.vcf.gz \\
            -P ${preAdapterDetail} \\
            -O ${TumorReplicateId}_${NormalReplicateId}_twicefitlered.vcf.gz && \\
        $GATK4 SelectVariants \\
            --tmp-dir ${params.tmpDir} \\
            --variant ${TumorReplicateId}_${NormalReplicateId}_twicefitlered.vcf.gz \\
            -R ${RefFasta} \\
            --exclude-filtered true \\
            --select 'vc.getGenotype(\"${TumorReplicateId}\").getAD().1 >= ${params.minAD}' \\
            --output ${TumorReplicateId}_${NormalReplicateId}_mutect2_final.vcf.gz
        """
}


// HaploTypeCaller
process 'HaploTypeCaller' {
/* 
    Call germline SNPs and indels via local re-assembly of haplotypes; normal sample
    germline variants are needed for generating phased vcfs for pVACtools
*/

    tag "$TumorReplicateId"

    // publishDir "$params.outputDir/$TumorReplicateId/03_haplotypeCaller/processing",
    // mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
        file(DBSNP),
        file(DBSNPIdx)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          database.DBSNP,
          database.DBSNPIdx ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from BaseRecalNormal_out_ch3
        .combine(
            SplitIntervals_out_ch5.flatten()
        )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_germline_${intervals}.vcf.gz"),
        file("${NormalReplicateId}_germline_${intervals}.vcf.gz.tbi"),
        file(Normalbam),
        file(Normalbai)
    ) into (
        HaploTypeCaller_out_ch0
    )


    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 --java-options ${params.JAVA_Xmx} HaplotypeCaller \\
        --tmp-dir ${params.tmpDir} \\
        -R ${RefFasta} \\
        -I ${Normalbam} \\
        -L ${intervals} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --dbsnp ${DBSNP} \\
        -O ${NormalReplicateId}_germline_${intervals}.vcf.gz
    """
}


// Run a Convolutional Neural Net to filter annotated germline variants
process 'CNNScoreVariants' {
/* 
    Run a Convolutional Neural Net to filter annotated germline variants; normal sample
    germline variants are needed for generating phased vcfs for pVACtools
*/

    // TODO: deal with this smarter
    conda 'assets/gatkcondaenv.yml'
    // conda '/data/projects/2019/ADSI/Exome_01/src/gatk-4.1.4.1_conda'
    // conda 'bioconda::gatk4-spark=4.1.4.1'

    tag "$TumorReplicateId"

    // publishDir "$params.outputDir/$TumorReplicateId/03_haplotypeCaller/processing",
    // mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(raw_germline_vcf),
        file(raw_germline_vcf_tbi),
        file(Normalbam),
        file(Normalbai)
    ) from HaploTypeCaller_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${raw_germline_vcf.baseName}_CNNScored.vcf.gz"),
        file("${raw_germline_vcf.baseName}_CNNScored.vcf.gz.tbi")
    ) into CNNScoreVariants_out_ch0


    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 CNNScoreVariants \\
        --tmp-dir ${params.tmpDir} \\
        -R ${RefFasta} \\
        -I ${Normalbam} \\
        -V ${raw_germline_vcf} \\
        -tensor-type read_tensor \\
        --inter-op-threads ${task.cpus} \\
        --intra-op-threads ${task.cpus} \\
        -O ${raw_germline_vcf.baseName}_CNNScored.vcf.gz
    """
}


process 'MergeHaploTypeCallerGermlineVCF' {
// Merge scattered filtered germline vcfs

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/04_haplotypeCaller/processing",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(filtered_germline_vcf),
        file(filtered_germline_vcf_tbi)
    ) from CNNScoreVariants_out_ch0
        .groupTuple(by: [0, 1])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_germline_CNNscored.vcf.gz"),
        file("${NormalReplicateId}_germline_CNNscored.vcf.gz.tbi"),
    ) into MergeHaploTypeCallerGermlineVCF_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}
    
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${filtered_germline_vcf.join(" I=")} \\
        O=${NormalReplicateId}_germline_CNNscored.vcf.gz
    """
}

// Apply a Convolutional Neural Net to filter annotated germline variants
process 'FilterGermlineVariantTranches' {
/* 
    Apply a Convolutional Neural Net to filter annotated germline variants; normal sample
    germline variants are needed for generating phased vcfs for pVACtools
*/

    // TODO: deal with this smarter
    conda 'assets/gatkcondaenv.yml'
    // conda '/data/projects/2019/ADSI/Exome_01/src/gatk-4.1.4.1_conda'

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/04_haplotypeCaller/",
        mode: params.publishDirMode

    input:
    set(
        file(MillsGold),
        file(MillsGoldIdx),
        file(HapMap),
        file(HapMapIdx),
        file(hcSNPS1000G),
        file(hcSNPS1000GIdx)
    ) from Channel.value(
        [ database.MillsGold,
          database.MillsGoldIdx,
          database.HapMap,
          database.HapMapIdx,
          database.hcSNPS1000G,
          database.hcSNPS1000GIdx ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(scored_germline_vcf),
        file(scored_germline_vcf_tbi)
    ) from MergeHaploTypeCallerGermlineVCF_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${scored_germline_vcf.simpleName}_Filtered.vcf.gz"),
        file("${scored_germline_vcf.simpleName}_Filtered.vcf.gz.tbi")
    ) into FilterGermlineVariantTranches_out_ch0


    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 FilterVariantTranches \\
        --tmp-dir ${params.tmpDir} \\
        -V ${scored_germline_vcf} \\
        --resource ${hcSNPS1000G} \\
        --resource ${HapMap} \\
        --resource ${MillsGold} \\
        --info-key CNN_2D \\
        --snp-tranche 99.95 \\
        --indel-tranche 99.4 \\
        --invalidate-previous-filters \\
        -O ${scored_germline_vcf.simpleName}_Filtered.vcf.gz
    """
}

/////////////////////////////////////////


// END HTC



/*
*********************************************
**             V A R S C A N               **
*********************************************
*/

process 'IndelRealignerTumorIntervals' {
/* 
 RealignerTargetCreator (GATK3): define intervals to target for local realignment
 IndelRealigner (GATK3): perform local realignment of reads around indels
*/

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai),
        file(list)
    ) from MarkDuplicatesTumor_out_ch2

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        file(KnownIndels),
        file(KnownIndelsIdx),
        file(MillsGold),
        file(MillsGoldIdx)
    ) from Channel.value(
        [ database.KnownIndels,
          database.KnownIndelsIdx,
          database.MillsGold,
          database.MillsGoldIdx ]
    )

    each file(interval) from SplitIntervals_out_ch3.flatten()

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_aligned_sort_mkdp_realign_${interval}.bam"),
        file("${TumorReplicateId}_aligned_sort_mkdp_realign_${interval}.bai")
    ) into IndelRealignerTumorIntervals_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \\
        -T RealignerTargetCreator \\
        --known ${MillsGold} \\
        --known ${KnownIndels} \\
        -R ${RefFasta} \\
        -L ${interval} \\
        -I ${bam} \\
        -o ${interval}_target.list \\
        -nt ${task.cpus} && \\
    $JAVA8 -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \\
        -T IndelRealigner \\
        -R ${RefFasta} \\
        -L ${interval} \\
        -I ${bam} \\
        -targetIntervals ${interval}_target.list \\
        -known ${KnownIndels} \\
        -known ${MillsGold} \\
        -nWayOut _realign_${interval}.bam && \\
    rm ${interval}_target.list
    """
}

process 'GatherRealignedBamFilesTumor' {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/processing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from IndelRealignerTumorIntervals_out_ch0
        .toSortedList({a, b -> a[2].baseName <=> b[2].baseName})
        .flatten()
        .collate(4)
        .groupTuple(by: [0,1])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_aligned_sort_mkdp_realign.bam"),
        file("${TumorReplicateId}_aligned_sort_mkdp_realign.bai") 
    ) into GatherRealignedBamFilesTumor_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 -XX:ParallelGCThreads=${task.cpus} ${params.JAVA_Xmx} -jar $PICARD GatherBamFiles \\
        TMP_DIR=${params.tmpDir} \\
        I=${bam.join(" I=")} \\
        O=${TumorReplicateId}_aligned_sort_mkdp_realign.bam \\
        CREATE_INDEX=true \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam}
    """
}

process 'BaseRecalTumorGATK3' {
/*
 FixMateInformation (Picard): verify mate-pair information between mates; ensure that
 all mate-pair information is in sync between reads and its pairs
*/

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/processing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from GatherRealignedBamFilesTumor_out_ch0

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch7

    set(
        file(DBSNP),
        file(DBSNPIdx),
        file(KnownIndels),
        file(KnownIndelsIdx),
        file(MillsGold),
        file(MillsGoldIdx)
    ) from Channel.value(
        [ database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx,
          database.MillsGold,
          database.MillsGoldIdx ]
    )

    output:
    file("${TumorReplicateId}_bqsr4.table")
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_recal4.bam"),
        file("${TumorReplicateId}_recal4.bam.bai"),
    ) into (
        BaseRecalTumorGATK3_out_ch0,
        BaseRecalTumorGATK3_out_ch1,
        BaseRecalTumorGATK3_out_ch2,
        BaseRecalTumorGATK3_out_ch3,
        BaseRecalTumorGATK3_out_ch4,
        BaseRecalTumorGATK3_out_ch5
    )


    ///
    script:
    """
    mkdir -p ${params.tmpDir}

    # TODO: re-add later when more stable. Commented out do to crashing randomly
    # $GATK4 BaseRecalibratorSpark \\
    #    --java-options '${params.JAVA_Xmx_spark}' \\
    #    --tmp-dir ${params.tmpDir} \\
    #    -I ${bam} \\
    #    -R ${RefFasta} \\
    #    -L ${IntervalsList} \\
    #    -O ${TumorReplicateId}_bqsr4.table \\
    #    --known-sites ${DBSNP} \\
    #    --known-sites ${KnownIndels} \\
    #    --known-sites ${MillsGold} \\
    #    --spark-master local[${task.cpus}] \\
    #    --conf 'spark.executor.cores=${task.cpus}' \\
    #    --conf 'spark.local.dir=${params.tmpDir}' && \\
    $GATK4 BaseRecalibrator \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${TumorReplicateId}_bqsr4.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold} && \\
    $GATK4 ApplyBQSRSpark \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${TumorReplicateId}_recal4.bam \\
        --bqsr-recal-file ${TumorReplicateId}_bqsr4.table \\
        --spark-master local[${task.cpus}] \\
        --conf 'spark.executor.cores=${task.cpus}' \\
        --conf 'spark.local.dir=${params.tmpDir}'
    """
}

process 'IndelRealignerNormalIntervals' {
/* 
RealignerTargetCreator (GATK3): define intervals to target for local realignment
IndelRealigner (GATK3): perform local realignment of reads around indels
*/

    tag "$NormalReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai),
        file(list)
    ) from MarkDuplicatesNormal_out_ch2

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        file(KnownIndels),
        file(KnownIndelsIdx),
        file(MillsGold),
        file(MillsGoldIdx)
    ) from Channel.value(
        [ database.KnownIndels,
          database.KnownIndelsIdx,
          database.MillsGold,
          database.MillsGoldIdx ]
    )

    each file(interval) from SplitIntervals_out_ch4.flatten()

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_aligned_sort_mkdp_realign_${interval}.bam"),
        file("${NormalReplicateId}_aligned_sort_mkdp_realign_${interval}.bai")
    ) into IndelRealignerNormalIntervals_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \\
        -T RealignerTargetCreator \\
        --known ${MillsGold} \\
        --known ${KnownIndels} \\
        -R ${RefFasta} \\
        -L ${interval} \\
        -I ${bam} \\
        -o ${interval}_target.list \\
        -nt ${task.cpus} && \\
    $JAVA8 -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \\
        -T IndelRealigner \\
        -R ${RefFasta} \\
        -L ${interval} \\
        -I ${bam} \\
        -targetIntervals ${interval}_target.list \\
        -known ${KnownIndels} \\
        -known ${MillsGold} \\
        -nWayOut _realign_${interval}.bam && \\
    rm ${interval}_target.list
    """
}

process 'GatherRealignedBamFilesNormal' {

    tag "$NormalReplicateId"

    publishDir "$params.outputDir/${TumorReplicateId}/03_varscan/processing/",
        mode: params.publishDirMode

    // publishDir "$params.outputDir/${TumorReplicateId[0]}/05_varscan/processing/",
    //     mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from IndelRealignerNormalIntervals_out_ch0
        .toSortedList({a, b -> a[2].baseName <=> b[2].baseName})
        .flatten()
        .collate(4)
        .groupTuple(by: [0,1])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_aligned_sort_mkdp_realign.bam"),
        file("${NormalReplicateId}_aligned_sort_mkdp_realign.bai") 
    ) into GatherRealignedBamFilesNormal_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 -XX:ParallelGCThreads=${task.cpus} ${params.JAVA_Xmx} -jar $PICARD GatherBamFiles \\
        TMP_DIR=${params.tmpDir} \\
        I=${bam.join(" I=")} \\
        O=${NormalReplicateId}_aligned_sort_mkdp_realign.bam \\
        CREATE_INDEX=true \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam}
    """
}

process 'BaseRecalNormalGATK3' {
/*
FixMateInformation (Picard): verify mate-pair information between mates; ensure that
all mate-pair information is in sync between reads and its pairs
*/

    tag "$NormalReplicateId"
    // There was a reason for this?
    publishDir "$params.outputDir/${TumorReplicateId}/03_varscan/processing/",
        mode: params.publishDirMode
    // publishDir "$params.outputDir/${TumorReplicateId[0]}/05_varscan/processing/",
    //     mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai)
    ) from GatherRealignedBamFilesNormal_out_ch0

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    file(IntervalsList) from preprocessIntervalList_out_ch8

    set(
        file(DBSNP),
        file(DBSNPIdx),
        file(KnownIndels),
        file(KnownIndelsIdx),
        file(MillsGold),
        file(MillsGoldIdx)
    ) from Channel.value(
        [ database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx,
          database.MillsGold,
          database.MillsGoldIdx ]
    )

    output:
    file("${NormalReplicateId}_bqsr4.table")
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${NormalReplicateId}_recal4.bam"),
        file("${NormalReplicateId}_recal4.bam.bai")
    ) into (
        BaseRecalNormalGATK3_out_ch0,
        BaseRecalNormalGATK3_out_ch1,
        BaseRecalNormalGATK3_out_ch2,
        BaseRecalNormalGATK3_out_ch3
    )


    ///
    script:
    """
    mkdir -p ${params.tmpDir}

    # TODO: re-add later when more stable. Commented out do to crashing randomly
    # $GATK4 BaseRecalibratorSpark \\
    #    --java-options '${params.JAVA_Xmx_spark}' \\
    #    --tmp-dir ${params.tmpDir} \\
    #    -I ${bam} \\
    #    -R ${RefFasta} \\
    #    -L ${IntervalsList} \\
    #    -O ${NormalReplicateId}_bqsr4.table \\
    #    --known-sites ${DBSNP} \\
    #    --known-sites ${KnownIndels} \\
    #    --known-sites ${MillsGold} \\
    #    --spark-master local[${task.cpus}] \\
    #    --conf 'spark.executor.cores=${task.cpus}' \\
    #    --conf 'spark.local.dir=${params.tmpDir}' && \\
    $GATK4 BaseRecalibrator \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${NormalReplicateId}_bqsr4.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold} && \\
    $GATK4 ApplyBQSRSpark \\
        --java-options '${params.JAVA_Xmx_spark}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${IntervalsList} \\
        -O ${NormalReplicateId}_recal4.bam \\
        --bqsr-recal-file ${NormalReplicateId}_bqsr4.table \\
        --spark-master local[${task.cpus}] \\
        --conf 'spark.executor.cores=${task.cpus}' \\
        --conf 'spark.local.dir=${params.tmpDir}'
    """
}

process 'VarscanSomaticScattered' {

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Tumorbam),
        file(Tumorbai),
        _,                    // unused TumorReplicateId from BaseRecalNormal_out_ch2
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from BaseRecalTumorGATK3_out_ch0
        .combine(BaseRecalNormalGATK3_out_ch0, by: 0)
        .combine(
            ScatteredIntervalListToBed_out_ch0.flatten()
        )

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan.snp.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan.indel.vcf")
    ) into VarscanSomaticScattered_out_ch0

    script:
    // awk filters at the end needed, found at least one occurence of "W" in Ref field of
    // varscan vcf (? wtf). Non ACGT seems to cause MergeVCF (picard) crashing
    """
    mkfifo ${TumorReplicateId}_${NormalReplicateId}_${intervals}_mpileup.fifo
    $SAMTOOLS mpileup \\
        -q 1 \\
        -f ${RefFasta} \\
        -l ${intervals} \\
        ${Normalbam} ${Tumorbam} > ${TumorReplicateId}_${NormalReplicateId}_${intervals}_mpileup.fifo &
    $JAVA8 ${params.JAVA_Xmx} -jar $VARSCAN somatic \\
        ${TumorReplicateId}_${NormalReplicateId}_${intervals}_mpileup.fifo \\
        ${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan_tmp \\
        --output-vcf 1 \\
        --mpileup 1 \\
        --min-coverage ${params.min_cov} \\
        --min-coverage-normal ${params.min_cov_normal} \\
        --min-coverage-tumor ${params.min_cov_tumor} \\
        --min-freq-for-hom ${params.min_freq_for_hom} \\
        --tumor-purity ${params.tumor_purity} \\
        --p-value ${params.somatic_pvalue} \\
        --somatic-p-value ${params.somatic_somaticpvalue} \\
        --strand-filter ${params.strand_filter} && \\
    rm -f ${TumorReplicateId}_${NormalReplicateId}_${intervals}_mpileup.fifo && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]/) { print } } else { print } }' \\
        ${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan_tmp.snp.vcf \\
        > ${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan.snp.vcf && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]+/) { print } } else { print } }' \\
        ${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan_tmp.indel.vcf \\
        > ${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan.indel.vcf

    rm -f ${TumorReplicateId}_${NormalReplicateId}_${intervals}_varscan_tmp.*
    """


}

process 'gatherVarscanVCFs' {
// Merge scattered Varscan vcfs

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/",
        mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(snp_vcf),
        file(indel_vcf)
    ) from VarscanSomaticScattered_out_ch0
        .toSortedList({a, b -> a[2].baseName <=> b[2].baseName})
        .flatten()
        .collate(4)
        .groupTuple(by: [0,1])


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.vcf")
    ) into gatherVarscanVCFs_out_ch0
 
    script:
    """
    mkdir -p ${params.tmpDir}
    
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${snp_vcf.join(" I=")} \\
        O=${TumorReplicateId}_${NormalReplicateId}_varscan.snp.vcf \\
        SEQUENCE_DICTIONARY=${RefDict}

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${indel_vcf.join(" I=")} \\
        O=${TumorReplicateId}_${NormalReplicateId}_varscan.indel.vcf \\
        SEQUENCE_DICTIONARY=${RefDict}

    """
}

process 'ProcessVarscan' {
// Filter variants by somatic status and confidences

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/processing",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(snp),
        file(indel)
    ) from gatherVarscanVCFs_out_ch0 // VarscanSomatic_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.Somatic.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.Somatic.hc.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.LOH.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.LOH.hc.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.Germline.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.Germline.hc.vcf")
    ) into ProcessVarscanSNP_out_ch0

    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.Somatic.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.Somatic.hc.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.LOH.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.LOH.hc.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.Germline.vcf"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.Germline.hc.vcf")
    ) into ProcessVarscanIndel_out_ch0

    script:
    """
    $JAVA8 ${params.JAVA_Xmx} -jar $VARSCAN processSomatic \\
        ${snp} \\
        --min-tumor-freq ${params.min_tumor_freq} \\
        --max-normal-freq ${params.max_normal_freq} \\
        --p-value ${params.processSomatic_pvalue} && \\
    $JAVA8 ${params.JAVA_Xmx} -jar $VARSCAN processSomatic \\
        ${indel} \\
        --min-tumor-freq ${params.min_tumor_freq} \\
        --max-normal-freq ${params.max_normal_freq} \\
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

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/processing",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bam),
        file(bai),
        _,
        file(snpSomatic),
        file(snpSomaticHc),
        file(snpLOH),
        file(snpLOHhc),
        file(snpGerm),
        file(snpGemHc),
        _,
        file(indelSomatic),
        file(indelSomaticHc),
        file(indelLOH),
        file(indelLOHhc),
        file(indelGerm),
        file(indelGemHc)
    ) from BaseRecalTumorGATK3_out_ch1
        .combine(ProcessVarscanSNP_out_ch0, by :0)
        .combine(ProcessVarscanIndel_out_ch0, by :0)

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.snp.Somatic.hc.filtered.vcf")
    ) into FilterVarscanSnp_out_ch0

    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.indel.Somatic.hc.filtered.vcf")
    ) into FilterVarscanIndel_out_ch0

    script:
    """
    cat ${snpSomaticHc} | \\
    awk '{if (!/^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' | \\
    $BAMREADCOUNT \\
        -q${params.min_map_cov} \\
        -b${params.min_base_q} \\
        -w1 \\
        -l /dev/stdin \\
        -f ${RefFasta} \\
        ${bam} | \\
    $JAVA8 ${params.JAVA_Xmx} -jar $VARSCAN fpfilter \\
        ${snpSomaticHc} \\
        /dev/stdin \\
        --output-file ${TumorReplicateId}_${NormalReplicateId}_varscan.snp.Somatic.hc.filtered.vcf && \\
    cat ${indelSomaticHc} | \\
    awk '{if (! /^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' | \\
    $BAMREADCOUNT \\
        -q${params.min_map_cov} \\
        -b${params.min_base_q} \\
        -w1 \\
        -l /dev/stdin \\
        -f ${RefFasta} ${bam} | \\
    $JAVA8 ${params.JAVA_Xmx} -jar $VARSCAN fpfilter \\
        ${indelSomaticHc} \\
        /dev/stdin \\
        --output-file ${TumorReplicateId}_${NormalReplicateId}_varscan.indel.Somatic.hc.filtered.vcf

    """
}

process 'MergeAndRenameSamplesInVarscanVCF' {
/*
    1. Merge filtered SNPS and INDELs from VarScan
    2. Rename the sample names (TUMOR/NORMAL) from varscan vcfs to the real samplenames
*/

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/",
        mode: params.publishDirMode

    input:
    file(RefDict) from Channel.value(reference.RefDict)

    set(
        TumorReplicateId,
        NormalReplicateId,
        VarScanSNP_VCF,
        VarScanINDEL_VCF
    ) from FilterVarscanSnp_out_ch0
        .combine(FilterVarscanIndel_out_ch0, by: [0, 1])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        val("varscan"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.Somatic.hc.filtered.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_varscan.Somatic.hc.filtered.vcf.gz.tbi")
    ) into (
        MergeAndRenameSamplesInVarscanVCF_out_ch0,
        MergeAndRenameSamplesInVarscanVCF_out_ch1
    ) 

    script:
    """

    $BGZIP -c ${VarScanSNP_VCF} > ${VarScanSNP_VCF}.gz
    $TABIX -p vcf ${VarScanSNP_VCF}.gz
    $BGZIP -c ${VarScanINDEL_VCF} > ${VarScanINDEL_VCF}.gz
    $TABIX -p vcf ${VarScanINDEL_VCF}.gz

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${VarScanSNP_VCF}.gz \\
        I=${VarScanINDEL_VCF}.gz \\
        O=${TumorReplicateId}_varscan_combined.vcf.gz \\
        SEQUENCE_DICTIONARY=${RefDict}

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD SortVcf \\
        TMP_DIR=${params.tmpDir} \\
        I=${TumorReplicateId}_varscan_combined.vcf.gz \\
        O=${TumorReplicateId}_varscan_combined_sorted.vcf.gz \\
        SEQUENCE_DICTIONARY=${RefDict}

    # rename samples in varscan vcf
    printf "TUMOR ${TumorReplicateId}\nNORMAL ${NormalReplicateId}\n" > vcf_rename_${TumorReplicateId}_${NormalReplicateId}_tmp

    $BCFTOOLS reheader \\
        -s vcf_rename_${TumorReplicateId}_${NormalReplicateId}_tmp \\
        ${TumorReplicateId}_varscan_combined_sorted.vcf.gz \\
        > ${TumorReplicateId}_${NormalReplicateId}_varscan.Somatic.hc.filtered.vcf.gz

    $TABIX -p vcf ${TumorReplicateId}_${NormalReplicateId}_varscan.Somatic.hc.filtered.vcf.gz
    rm -f vcf_rename_${TumorReplicateId}_${NormalReplicateId}_tmp

    """

}

/*
*********************************************
**             M U T E C T 1               **
*********************************************
*/

process 'Mutect1scattered' {
// Mutect1: calls SNPS from tumor and matched normal sample

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Tumorbam),
        file(Tumorbai),
        _,                    // unused TumorReplicateId from BaseRecalNormal_out_ch2
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from BaseRecalTumorGATK3_out_ch3
        .combine(BaseRecalNormalGATK3_out_ch1, by: 0)
        .combine(
            SplitIntervals_out_ch6.flatten()
        )

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )


    set(
        file(DBSNP),
        file(DBSNPIdx),
        file(Cosmic),
        file(CosmicIdx)
    ) from Channel.value(
        [ database.DBSNP,
          database.DBSNPIdx,
          database.Cosmic,
          database.CosmicIdx ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${intervals}.raw.vcf.gz"),
        file("${TumorReplicateId}_${intervals}.raw.stats.txt"),
        file("${TumorReplicateId}_${intervals}.raw.vcf.gz.idx")
    ) into Mutect1scattered_out_ch0

    script:
    """
    mkdir -p ${params.tmpDir}

    $JAVA7 ${params.JAVA_Xmx} -Djava.io.tmpdir=${params.tmpDir} -jar $MUTECT1 \\
        --analysis_type MuTect \\
        --reference_sequence ${RefFasta} \\
        --cosmic ${Cosmic} \\
        --dbsnp ${DBSNP} \\
        -L ${intervals} \\
        --input_file:normal ${Normalbam} \\
        --input_file:tumor ${Tumorbam} \\
        --out ${TumorReplicateId}_${intervals}.raw.stats.txt \\
        --vcf ${TumorReplicateId}_${intervals}.raw.vcf.gz
    """
}

process 'gatherMutect1VCFs' {
// Merge scattered Mutect1 vcfs

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect1/",
        mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(vcf),
        file(stats),
        file(idx)
    ) from Mutect1scattered_out_ch0
        .toSortedList({a, b -> a[2].baseName <=> b[2].baseName})
        .flatten()
        .collate(5)
        .groupTuple(by: [0,1])


    output:
    file("${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.vcf.gz")
    file("${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.vcf.gz.tbi")
    file("${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.stats.txt")

    set(
        TumorReplicateId,
        NormalReplicateId,
        val("mutect1"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect1_final.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect1_final.vcf.gz.tbi")
    ) into (
        Mutect1_out_ch0,
        Mutect1_out_ch1
    )


    script:
    """
    mkdir -p ${params.tmpDir}
    
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${vcf.join(" I=")} \\
        O=${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.vcf.gz

    $GATK4 SelectVariants \\
        --tmp-dir ${params.tmpDir} \\
        --variant ${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.vcf.gz \\
        -R ${RefFasta} \\
        --exclude-filtered true \\
        --select 'vc.getGenotype(\"${TumorReplicateId}\").getAD().1 >= ${params.minAD}' \\
        --output ${TumorReplicateId}_${NormalReplicateId}_mutect1_final.vcf.gz


    head -2 ${stats[0]} > ${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.stats.txt
    tail -q -n +3 ${stats.join(" ")} >> ${TumorReplicateId}_${NormalReplicateId}_mutect1_raw.stats.txt
    """
}


// Strelka2 and Manta
process 'MantaSomaticIndels' {
    conda 'bioconda::manta=1.6.0'

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_manta_somatic",
        mode: params.publishDirMode
    
    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Tumorbam),
        file(Tumorbai),
        _,                    // unused TumorReplicateId from manta_inputNormal
        file(Normalbam),
        file(Normalbai),
        file(RegionsBedGz),
        file(RegionsBedGzTbi)
    ) from BaseRecalTumorGATK3_out_ch4
        .combine(BaseRecalNormalGATK3_out_ch2, by:0)
        .combine(RegionsBedToTabix_out_ch1)

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_candidateSmallIndels.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_candidateSmallIndels.vcf.gz.tbi")
    ) into MantaSomaticIndels_out_ch0

    file("${TumorReplicateId}_${NormalReplicateId}_diploidSV.vcf.gz")
    file("${TumorReplicateId}_${NormalReplicateId}_diploidSV.vcf.gz.tbi")
    file("${TumorReplicateId}_${NormalReplicateId}_candidateSV.vcf.gz")
    file("${TumorReplicateId}_${NormalReplicateId}_candidateSV.vcf.gz.tbi")
    file("${TumorReplicateId}_${NormalReplicateId}_candidateSmallIndels.vcf.gz")
    file("${TumorReplicateId}_${NormalReplicateId}_candidateSmallIndels.vcf.gz.tbi")
    file("${TumorReplicateId}_${NormalReplicateId}_svCandidateGenerationStats.tsv")

    script:
    """
    configManta.py --tumorBam ${Tumorbam} --normalBam  ${Normalbam} \\
        --referenceFasta ${RefFasta} \\
        --runDir manta_${TumorReplicateId} --callRegions ${RegionsBedGz} --exome &&
    manta_${TumorReplicateId}/runWorkflow.py -m local -j  ${task.cpus} && 
    cp manta_${TumorReplicateId}/results/variants/diploidSV.vcf.gz ${TumorReplicateId}_${NormalReplicateId}_diploidSV.vcf.gz
    cp manta_${TumorReplicateId}/results/variants/diploidSV.vcf.gz.tbi ${TumorReplicateId}_${NormalReplicateId}_diploidSV.vcf.gz.tbi
    cp manta_${TumorReplicateId}/results/variants/candidateSV.vcf.gz ${TumorReplicateId}_${NormalReplicateId}_candidateSV.vcf.gz
    cp manta_${TumorReplicateId}/results/variants/candidateSV.vcf.gz.tbi ${TumorReplicateId}_${NormalReplicateId}_candidateSV.vcf.gz.tbi
    cp manta_${TumorReplicateId}/results/variants/candidateSmallIndels.vcf.gz ${TumorReplicateId}_${NormalReplicateId}_candidateSmallIndels.vcf.gz
    cp manta_${TumorReplicateId}/results/variants/candidateSmallIndels.vcf.gz.tbi ${TumorReplicateId}_${NormalReplicateId}_candidateSmallIndels.vcf.gz.tbi
    cp manta_${TumorReplicateId}/results/stats/svCandidateGenerationStats.tsv ${TumorReplicateId}_${NormalReplicateId}_svCandidateGenerationStats.tsv
    """
}


process StrelkaSomatic {
    conda 'bioconda::strelka=2.9.10'
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_strelka_somatic/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Tumorbam),
        file(Tumorbai),
        _,                    // unused NormalReplicateId from strelka_inputNormal
        file(Normalbam),
        file(Normalbai),
        _,                    // unused NormalReplicateId from MantaSomaticIndels_out_ch0
        file(manta_indel),
        file(manta_indel_tbi),
        file(RegionsBedGz),
        file(RegionsBedGzTbi)
    ) from BaseRecalTumorGATK3_out_ch5
        .combine(BaseRecalNormalGATK3_out_ch3, by:0)
        .combine(MantaSomaticIndels_out_ch0, by:0)
        .combine(RegionsBedToTabix_out_ch0)

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    output:
    tuple(
        TumorReplicateId,
        NormalReplicateId,
        val("strelka"),
        file("${TumorReplicateId}_${NormalReplicateId}_strelka_somatic_final.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_strelka_somatic_final.vcf.gz.tbi")
    ) into (
        StrelkaSomatic_out_ch0,
        StrelkaSomatic_out_ch1
    )
    file("${TumorReplicateId}_${NormalReplicateId}_strelka_combined_somatic.vcf.gz")
    file("${TumorReplicateId}_${NormalReplicateId}_runStats.tsv")
    file("${TumorReplicateId}_${NormalReplicateId}_runStats.xml")
    
    script:
    """
    configureStrelkaSomaticWorkflow.py --tumorBam ${Tumorbam} --normalBam  ${Normalbam} \\
        --referenceFasta ${RefFasta} \\
        --indelCandidates ${manta_indel} \\
        --runDir strelka_${TumorReplicateId} --callRegions ${RegionsBedGz} --exome &&
    strelka_${TumorReplicateId}/runWorkflow.py -m local -j ${task.cpus} &&
    cp strelka_${TumorReplicateId}/results/variants/somatic.indels.vcf.gz ${TumorReplicateId}_${NormalReplicateId}_somatic.indels.vcf.gz
    cp strelka_${TumorReplicateId}/results/variants/somatic.indels.vcf.gz.tbi ${TumorReplicateId}_${NormalReplicateId}_somatic.indels.vcf.gz.tbi
    cp strelka_${TumorReplicateId}/results/variants/somatic.snvs.vcf.gz ${TumorReplicateId}_${NormalReplicateId}_somatic.snvs.vcf.gz
    cp strelka_${TumorReplicateId}/results/variants/somatic.snvs.vcf.gz.tbi ${TumorReplicateId}_${NormalReplicateId}_somatic.snvs.vcf.gz.tbi
    cp strelka_${TumorReplicateId}/results/stats/runStats.tsv ${TumorReplicateId}_${NormalReplicateId}_runStats.tsv
    cp strelka_${TumorReplicateId}/results/stats/runStats.xml ${TumorReplicateId}_${NormalReplicateId}_runStats.xml

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${TumorReplicateId}_${NormalReplicateId}_somatic.snvs.vcf.gz \\
        I=${TumorReplicateId}_${NormalReplicateId}_somatic.indels.vcf.gz \\
        O=${TumorReplicateId}_${NormalReplicateId}_strelka_combined.vcf.gz \\
        SEQUENCE_DICTIONARY=${RefDict}

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD SortVcf \\
        TMP_DIR=${params.tmpDir} \\
        I=${TumorReplicateId}_${NormalReplicateId}_strelka_combined.vcf.gz \\
        O=${TumorReplicateId}_${NormalReplicateId}_strelka_combined_sorted.vcf.gz \\
        SEQUENCE_DICTIONARY=${RefDict}

    # rename samples in varscan vcf
    printf "TUMOR ${TumorReplicateId}\nNORMAL ${NormalReplicateId}\n" > vcf_rename_${TumorReplicateId}_${NormalReplicateId}_tmp

    $BCFTOOLS reheader \\
        -s vcf_rename_${TumorReplicateId}_${NormalReplicateId}_tmp \\
        ${TumorReplicateId}_${NormalReplicateId}_strelka_combined_sorted.vcf.gz \\
        > ${TumorReplicateId}_${NormalReplicateId}_strelka_combined_somatic.vcf.gz

    $TABIX -p vcf ${TumorReplicateId}_${NormalReplicateId}_strelka_combined_somatic.vcf.gz
    rm -f vcf_rename_${TumorReplicateId}_${NormalReplicateId}_tmp

    $GATK4 SelectVariants \\
        --tmp-dir ${params.tmpDir} \\
        --variant ${TumorReplicateId}_${NormalReplicateId}_strelka_combined_somatic.vcf.gz \\
        -R ${RefFasta} \\
        --exclude-filtered true \\
        --output ${TumorReplicateId}_${NormalReplicateId}_strelka_somatic_final.vcf.gz


    """

}


// END Strelka2 and Manta

process 'mkHCsomaticVCF' {
/*
    Creates a VCF that is based on the priority caller (e.g. mutect2) vcf but contains only variants
    that are confirmed by any of the two confirming callers (e..g. mutect1, varscan)
*/

    // TODO: deal with this smarter
    conda 'assets/gatkcondaenv.yml'
    // conda '/data/projects/2019/ADSI/Exome_01/src/gatk-4.1.4.1_conda'

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/05_hcVCF/",
        mode: params.publishDirMode

    input:
    file(RefDict) from Channel.value(reference.RefDict)

    set(
        TumorReplicateId,
        NormalReplicateId,
        _,
        Mutect2_VCF,
        Mutect2_Idx,
        _,
        Mutect1_VCF,
        Mutect1_Idx,
        _,
        VarScan_VCF,
        VarScan_Idx,
        _,
        Strelka_VCF,
        Strelka_IDX
    ) from FilterMutect2_out_ch1
        .combine(MergeAndRenameSamplesInVarscanVCF_out_ch1, by: [0, 1])
        .combine(Mutect1_out_ch1, by: [0, 1])
        .combine(StrelkaSomatic_out_ch0, by: [0, 1])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        val("hc"),
        file("${TumorReplicateId}_${NormalReplicateId}_Somatic.hc.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_Somatic.hc.vcf.gz.tbi")
    ) into (
        mkHCsomaticVCF_out_ch0,
        mkHCsomaticVCF_out_ch1,
        mkHCsomaticVCF_out_ch2
    ) 

    script:
    """
    make_hc_vcf.py --priority ${params.priorityCaller} \\
        --m2_vcf ${Mutect2_VCF} \\
        --m1_vcf ${Mutect1_VCF} \\
        --vs_vcf ${VarScan_VCF} \\
        --st_vcf ${Strelka_VCF} \\
        --out_vcf ${TumorReplicateId}_${NormalReplicateId}_Somatic.hc.vcf \\
        --out_single_vcf ${TumorReplicateId}_${NormalReplicateId}_Somatic.single.vcf

    $BGZIP -c ${TumorReplicateId}_${NormalReplicateId}_Somatic.hc.vcf > ${TumorReplicateId}_${NormalReplicateId}_Somatic.hc.vcf.gz
    $TABIX -p vcf ${TumorReplicateId}_${NormalReplicateId}_Somatic.hc.vcf.gz
    """

}

process 'VepTab' {
// Variant Effect Prediction: using ensembl vep

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/06_vep/",
        mode: params.publishDirMode

    input:

    set(
        TumorReplicateId,
        NormalReplicateId,
        CallerName,
        file(Vcf),
        file(Idx),
    ) from FilterMutect2_out_ch0
        .concat(MergeAndRenameSamplesInVarscanVCF_out_ch0)
        .concat(Mutect1_out_ch0)
        .concat(StrelkaSomatic_out_ch1)
        .concat(mkHCsomaticVCF_out_ch0)
        .flatten()
        .collate(5)

    output:
    file("${TumorReplicateId}_${NormalReplicateId}_${CallerName}_vep.txt")
    file("${TumorReplicateId}_${NormalReplicateId}_${CallerName}_vep_summary.html")
    
    script:
    """
    $PERL $VEP -i ${Vcf} \\
        -o ${TumorReplicateId}_${NormalReplicateId}_${CallerName}_vep.txt \\
        --fork ${task.cpus} \\
        --stats_file ${TumorReplicateId}_${NormalReplicateId}_${CallerName}_vep_summary.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --dir ${params.vep_dir} \\
        --cache --dir_cache ${params.vep_cache} \\
        --fasta ${params.VepFasta} \\
        --format "vcf" \\
        ${params.vep_options} \\
        --tab
    """
}    

// CREATE phased VCF
process 'mkPhasedVCF' {
/* 
    make phased vcf for pVACseq using tumor and germliine variants:
    based on https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/proximal_vcf.html
*/

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/07_PhasedVCF",
        mode: params.publishDirMode

    input:
    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    set(
        TumorReplicateId,
        NormalReplicateId,
        file(tumorBAM),
        file(tumorBAI),
        file(germlineVCF),
        file(germlineVCFidx),
        _,
        file(tumorVCF),
        file(tumorVCFidx)
    ) from BaseRecalTumorGATK3_out_ch2
        .combine(FilterGermlineVariantTranches_out_ch0, by: [0,1])
        .combine(mkHCsomaticVCF_out_ch1, by: [0,1])             // uses confirmed mutect2 variants

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_vep_phased.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_vep_phased.vcf.gz.tbi"),
        file("${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf.gz.tbi")
    ) into (
        mkPhasedVCF_out_ch0, mkPhasedVCF_out_pVACseqch0
    )
    file("${TumorReplicateId}_${NormalReplicateId}_tumor_reference.fa")
    file("${TumorReplicateId}_${NormalReplicateId}_tumor_mutated.fa")

    script:
    """
    mkdir -p ${params.tmpDir}

    $GATK4 --java-options ${params.JAVA_Xmx} SelectVariants \\
        --tmp-dir ${params.tmpDir} \\
        -R ${RefFasta} \\
        -V ${tumorVCF} \\
        --sample-name ${TumorReplicateId} \\
        -O ${TumorReplicateId}_${NormalReplicateId}_tumor.vcf.gz

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD RenameSampleInVcf \\
        TMP_DIR=${params.tmpDir} \\
        I=${germlineVCF} \\
        NEW_SAMPLE_NAME=${TumorReplicateId} \\
        O=${NormalReplicateId}_germlineVAR_rename2tumorID.vcf.gz

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD MergeVcfs \\
        TMP_DIR=${params.tmpDir} \\
        I=${TumorReplicateId}_${NormalReplicateId}_tumor.vcf.gz \\
        I=${NormalReplicateId}_germlineVAR_rename2tumorID.vcf.gz \\
        O=${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined.vcf.gz

    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD SortVcf \\
        TMP_DIR=${params.tmpDir} \\
        I=${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined.vcf.gz \\
        O=${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted.vcf.gz \\
        SEQUENCE_DICTIONARY=${RefDict}

    $PERL $VEP -i ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted.vcf.gz \\
        -o ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep.vcf \\
        --fork ${task.cpus} \\
        --stats_file ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep_summary.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --cache \\
        --dir ${params.vep_dir} \\
        --dir_cache ${params.vep_cache} \\
        --fasta ${params.VepFasta} \\
        --pick --plugin Downstream --plugin Wildtype \\
        --symbol --terms SO --transcript_version --tsl \\
        --vcf

    $PERL $VEP -i ${tumorVCF} \\
        -o ${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf \\
        --fork ${task.cpus} \\
        --stats_file ${TumorReplicateId}_${NormalReplicateId}_tumor_vep_summary.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --cache \\
        --dir ${params.vep_dir} \\
        --dir_cache ${params.vep_cache} \\
        --fasta ${params.VepFasta} \\
        --pick --plugin Downstream --plugin Wildtype \\
        --plugin ProteinSeqs,${TumorReplicateId}_${NormalReplicateId}_tumor_reference.fa,${TumorReplicateId}_${NormalReplicateId}_tumor_mutated.fa \\
        --symbol --terms SO --transcript_version --tsl \\
        --vcf

    $BGZIP -c ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep.vcf \\
        > ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep.vcf.gz
    
    $TABIX -p vcf ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep.vcf.gz

    $BGZIP -c ${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf \\
        > ${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf.gz
    
    $TABIX -p vcf ${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf.gz


    $JAVA8 -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${params.tmpDir} -jar $GATK3 \\
        -T ReadBackedPhasing \\
        -R ${RefFasta} \\
        -I ${tumorBAM} \\
        -V ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep.vcf.gz \\
        -L ${TumorReplicateId}_${NormalReplicateId}_germlineVAR_combined_sorted_vep.vcf.gz \\
        -o ${TumorReplicateId}_${NormalReplicateId}_vep_phased.vcf.gz 
    """
}
// END CREATE phased VCF

// HLA TYPING

process 'mhc_extract' {
    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(tumor_BAM_aligned_sort_mkdp),
        file(tumor_BAI_aligned_sort_mkdp),
        _
    ) from MarkDuplicatesTumor_out_ch3


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${mhcReads_1}"),
        file("${mhcReads_2}"),
        _,      // unused so far
    ) into (
        reads_tumor_hla_ch,
        reads_tumor_hlaHD_ch
    )

    script:
    mhc_region = params.HLA_HD_genome_version ? params.MHC_genomic_region[ params.HLA_HD_genome_version ].region ?: false : false
    
    if (!mhc_region) {
        exit 1, "MHC region not found for genome version: ${params.HLA_HD_genome_version}" 
    }

    mhcReads_1 = (single_end) ? TumorReplicateId + "_reads_mhc.fastq.gz" : TumorReplicateId + "_reads_mhc_R1.fastq.gz"
    mhcReads_2 = (single_end) ? val("NO_FILE") : TumorReplicateId + "_reads_mhc_R2.fastq.gz"

    if(single_end)
        """
        mkfifo unmapped_bam
        mkfifo mhc_mapped_bam
        mkfifo R.fastq

        $SAMTOOLS  view -@4 -h -b -u -f 4 ${tumor_BAM_aligned_sort_mkdp} > unmapped_bam &
        $SAMTOOLS  view -@4 -h -b -u ${tumor_BAM_aligned_sort_mkdp} ${mhc_region} > mhc_mapped_bam &

        $SAMTOOLS merge -@4 -u - mhc_mapped_bam unmapped_bam | \\
            $SAMTOOLS sort -@4 -n - | \\
            $SAMTOOLS fastq -@2 -0 R.fastq \\
            -i - &
        $PERL -ple 'if ((\$. % 4) == 1) { s/\$/ 1:N:0:NNNNNNNN/; }' R.fastq | gzip -1 > ${TumorReplicateId}_reads_mhc.fastq.gz

        wait

        rm -f unmapped_bam mhc_mapped_bam R.fastq
        """
    else
        """
        mkfifo unmapped_bam
        mkfifo mhc_mapped_bam
        mkfifo R1.fastq
        mkfifo R2.fastq

        $SAMTOOLS  view -@4 -h -b -u -f 4 ${tumor_BAM_aligned_sort_mkdp} > unmapped_bam &
        $SAMTOOLS  view -@4 -h -b -u ${tumor_BAM_aligned_sort_mkdp} ${mhc_region} > mhc_mapped_bam &

        $SAMTOOLS merge -@4 -u - mhc_mapped_bam unmapped_bam | \\
            $SAMTOOLS sort -@4 -n - | \\
            $SAMTOOLS fastq -@2 -1 R1.fastq -2 R2.fastq -s /dev/null -0 /dev/null \\
            -i - &
        $PERL -ple 'if ((\$. % 4) == 1) { s/\$/ 1:N:0:NNNNNNNN/; }' R1.fastq | gzip -1 > ${TumorReplicateId}_reads_mhc_R1.fastq.gz &
        $PERL -ple 'if ((\$. % 4) == 1) { s/\$/ 2:N:0:NNNNNNNN/; }' R2.fastq | gzip -1 > ${TumorReplicateId}_reads_mhc_R2.fastq.gz &

        wait

        rm -f unmapped_bam mhc_mapped_bam R1.fastq R2.fastq
        """
}


/*
*********************************************
**             O P T I T Y P E             **
*********************************************
*/

/*
 * Preparation Step - Pre-mapping against HLA
 *
 * In order to avoid the internal usage of RazerS from within OptiType when
 * the input files are of type `fastq`, we perform a pre-mapping step
 * here with the `yara` mapper, and map against the HLA reference only.
 *
 */

process 'pre_map_hla' {
    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(readsFWD),
        file(readsREV),
        sampleGroup,      // unused so far
    ) from reads_tumor_hla_ch
    val yaraIdx from Channel.value(reference.YaraIndex)

    output:
    set (
        TumorReplicateId,
        file("mapped_{1,2}.bam")
    ) into fished_reads
    
    script:
    if(single_end) {
        yara_cpus = ((task.cpus - 2).compareTo(2) == -1) ? 2 : (task.cpus - 2)
        samtools_cpus = ((task.cpus - yara_cpus).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus)
    } else {
        yara_cpus = ((task.cpus - 6).compareTo(2) == -1) ? 2 : (task.cpus - 6)
        samtools_cpus =  (((task.cpus - yara_cpus).div(3)).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus).div(3)
    }
    if (single_end)
        """
        $YARA -e 3 -t $yara_cpus -f bam ${yaraIdx} ${readsFWD} | \\
            $SAMTOOLS view -@ $samtools_cpus -h -F 4 -b1 - > mapped_1.bam
        """
    else
        """
        mkfifo R1 R2
        $YARA -e 3 -t $yara_cpus -f bam ${yaraIdx} ${readsFWD} ${readsREV} | \\
            $SAMTOOLS view -@ $samtools_cpus -h -F 4 -b1 | \\
            tee R1 R2 > /dev/null &
        $SAMTOOLS view -@ $samtools_cpus -h -f 0x40 -b1 R1 > mapped_1.bam &
        $SAMTOOLS view -@ $samtools_cpus -h -f 0x80 -b1 R2 > mapped_2.bam &
        wait
        rm -f R1 R2
        """
}

/*
 * STEP 2 - Run Optitype
 *
 * This is the major process, that formulates the IP and calls the selected
 * IP solver.
 *
 * Ouput formats: <still to enter>
 */

process 'OptiType' {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/08_OptiType/",
        mode: params.publishDirMode

    input:
    set (
        TumorReplicateId,
        file(reads)
    ) from fished_reads.collect()

    output:
    set (
        TumorReplicateId,
        file("**/*_result.tsv")
    ) into optitype_output
    file("**")
    // file("**/*_result.tsv") into optitype_output

    script:
    """
    $PYTHON $OPTITYPE -i ${reads} -e 1 -b 0.009 --dna -o .
    """
}

/*
*********************************************
**             H L A - H D                 **
*********************************************
*/

process 'run_hla_hd' {

    tag "$TumorReplicateId"

    conda 'assets/hlahdenv.yml'

    publishDir "$params.outputDir/$TumorReplicateId/09_HLA_HD/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        _,
        file(readsFWD),
        file(readsREV),
        sampleGroup,      // unused so far
    ) from reads_tumor_hlaHD_ch
    val frData from Channel.value(reference.HLAHDFreqData)
    file gSplit from Channel.value(reference.HLAHDGeneSplit)
    val dict from Channel.value(reference.HLAHDDict)

    output:
    set (
        TumorReplicateId,
        file("**/*_final.result.txt")
    ) into hlahd_output

    script:
    hlahd_p = Channel.value(params.HLAHD_PATH).getVal()

    if (single_end)
        """
        export PATH=\$PATH:$hlahd_p
        $HLAHD -t ${task.cpus} \\
            -m 50 \\
            -f ${frData} ${readsFWD} ${readsFWD} \\
            ${gSplit} ${dict} $TumorReplicateId .
        """
    else
        """
        export PATH=\$PATH:$hlahd_p
        $HLAHD -t ${task.cpus} \\
            -m 50 \\
            -f ${frData} ${readsFWD} ${readsREV} \\
            ${gSplit} ${dict} $TumorReplicateId .
        """
}
// END HLA TYPING

// NeoAntigen predictions

/*
*********************************************
**      N E O F U S E / P V A C S E Q      **
*********************************************
*/

/*
Prediction of gene fusion neoantigens with Neofuse and calculation of TPM values
*/

process Neofuse {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/10_NeoFuse/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        readRNAFWD,
        readRNAREV
    ) from reads_tumor_neofuse_ch
    file STARidx from file(reference.STARidx)
    file RefFASTA from file(reference.RefFASTA)
    file AnnoFile from file(reference.AnnoFile)

    output:
    set (
        TumorReplicateId,
        file("**/*.tpm.txt") // TODO: consider changing wildcards to expected dirnames and filenames
    ) into tpm_file
    set (
        TumorReplicateId,
        file("./${TumorReplicateId}/STAR/${TumorReplicateId}.Aligned.sortedByCoord.out.bam"),
        file("./${TumorReplicateId}/STAR/${TumorReplicateId}.Aligned.sortedByCoord.out.bam.bai")
    ) into star_bam_file
    path("${TumorReplicateId}")


    script:
    if(single_end_RNA)
        """
        NeoFuse_single -1 ${readRNAFWD} \\
            -d ${TumorReplicateId} \\
            -o . \\
            -m ${params.pepMin_length} \\
            -M ${params.pepMax_length} \\
            -n ${task.cpus} \\
            -t ${params.IC50_Threshold} \\
            -T ${params.rank} \\
            -c ${params.conf_lvl} \\
            -s ${STARidx} \\
            -g ${RefFASTA} \\
            -a ${AnnoFile} \\
            -N ${params.netMHCpan}
        """
    else
        """
        NeoFuse_single -1 ${readRNAFWD} -2 ${readRNAREV} \\
            -d ${TumorReplicateId} \\
            -o . \\
            -m ${params.pepMin_length} \\
            -M ${params.pepMax_length} \\
            -n ${task.cpus} \\
            -t ${params.IC50_Threshold} \\
            -T ${params.rank} \\
            -c ${params.conf_lvl} \\
            -s ${STARidx} \\
            -g ${RefFASTA} \\
            -a ${AnnoFile} \\
            -N ${params.netMHCpan}
        """		
}

/*
Add the gene ID (required by vcf-expression-annotator) to the TPM file
*/
process add_geneID {
    
    tag "$TumorReplicateId"

    input:
    set (
        TumorReplicateId,
        tpm
    ) from tpm_file
    // file tpm from tpm_file
    file AnnoFile from file(reference.AnnoFile)

    output:
    set (
        TumorReplicateId,
        file("*.tpm_final.txt")
    ) into final_file

    script:
    """
    NameToID.py -i ${tpm} -a ${AnnoFile} -o .
    """
}

/*
Add gene expression info to the VEP annotated, phased VCF file
*/

process gene_annotator {

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        _,
        _,
        vep_somatic_vcf_gz,
        vep_somatic_vcf_gz_tbi,
        final_tpm,
        RNA_bam,
        RNA_bai,
    ) from mkPhasedVCF_out_ch0
        .combine(final_file, by: 0)
        .combine(star_bam_file, by: 0)
    
    file(RefFasta) from file(reference.RefFASTA)
    file(RefIdx) from file(reference.RefIDX)

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_vep_somatic_gx.vcf.gz"),
        file("${TumorReplicateId}_vep_somatic_gx.vcf.gz.tbi")
    ) into vcf_vep_ex_gz

    script:
    """
    vcf-expression-annotator \\
        -i GeneID \\
        -e TPM \\
        -s ${TumorReplicateId} \\
        ${vep_somatic_vcf_gz} ${final_tpm} \\
        custom gene \\
        -o ./${TumorReplicateId}_vep_somatic_gx_tmp.vcf
    bgzip -f ${TumorReplicateId}_vep_somatic_gx_tmp.vcf
    tabix -p vcf ${TumorReplicateId}_vep_somatic_gx_tmp.vcf.gz

    vt decompose \\
        -s ${TumorReplicateId}_vep_somatic_gx_tmp.vcf.gz \\
        -o ${TumorReplicateId}_vep_somatic_gx_dec_tmp.vcf.gz

    bam_readcount_helper.py \\
        ${TumorReplicateId}_vep_somatic_gx_dec_tmp.vcf.gz \\
        ${TumorReplicateId} \\
        ${RefFasta} \\
        ${RNA_bam} \\
        ./

    vcf-readcount-annotator \\
        -s ${TumorReplicateId} \\
        -t snv \\
        -o ${TumorReplicateId}_vep_somatic_gx_dec_snv_rc_tmp.vcf \\
        ${TumorReplicateId}_vep_somatic_gx_dec_tmp.vcf.gz \\
        ${TumorReplicateId}_bam_readcount_snv.tsv \\
        RNA

    vcf-readcount-annotator \\
        -s ${TumorReplicateId} \\
        -t indel \\
        -o ${TumorReplicateId}_vep_somatic_gx.vcf \\
        ${TumorReplicateId}_vep_somatic_gx_dec_snv_rc_tmp.vcf \\
        ${TumorReplicateId}_bam_readcount_indel.tsv \\
        RNA

    bgzip -f ${TumorReplicateId}_vep_somatic_gx.vcf
    tabix -p vcf ${TumorReplicateId}_vep_somatic_gx.vcf.gz

    """
}

/*
Get the HLA types from OptiType and HLA-HD ouput as a "\n" seperated list.
To be used as input for pVACseq
*/

process get_vhla {
    tag "$TumorReplicateId"

    input:
    set (
        TumorReplicateId,
        opti_out,
        hlahd_out
    ) from optitype_output
        .combine(hlahd_output, by: 0)
    // file(opti_out) from optitype_output
    // file(hlahd_out) from hlahd_output
    
    output:
    set (
        TumorReplicateId,
        file("${TumorReplicateId}_hlas.txt")
    ) into hlas
    // file("${TumorReplicateId}_hlas.txt")
    // stdout hlas

    script:
    """
    HLA_parser.py \\
        --opti_out ${opti_out} \\
        --hlahd_out ${hlahd_out} \\
        --ref_hlas ${params.valid_HLAs} \\
        > ./${TumorReplicateId}_hlas.txt
    """
}

/*
Run pVACseq
*/

process 'pVACseq' {
    tag "$TumorReplicateId"
    cache false
    // remove publish when done with debuging
    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        vep_phased_vcf_gz,
        vep_phased_vcf_gz_tbi,
        _,
        _,
        _,
        anno_vcf,
        anno_vcf_tbi,
        hla_types
    ) from mkPhasedVCF_out_pVACseqch0
        .combine(vcf_vep_ex_gz, by: 0)
        .combine(hlas.splitText(), by: 0)

    output:
    set(
        TumorReplicateId,
        file("**/MHC_Class_I/*.filtered.tsv"),
        file("**/MHC_Class_I/*.all_epitopes.tsv")
    ) optional true into mhcI_out_f

    set(
        TumorReplicateId,
        file("**/MHC_Class_II/*.filtered.tsv"),
        file("**/MHC_Class_II/*.all_epitopes.tsv")
    ) optional true into mhcII_out_f


    script:
    hla_type = (hla_types - ~/\n/)
    NetChop = params.use_NetChop ? "--net-chop-method cterm" : ""
    NetMHCstab = params.use_NetMHCstab ? "--netmhc-stab" : ""
    """
    pvacseq run \\
        --iedb-install-directory /opt/iedb \\
        -t ${task.cpus} \\
        -p ${vep_phased_vcf_gz} \\
        -e1 ${params.mhci_epitope_len} \\
        -e2 ${params.mhcii_epitope_len} \\
        --normal-sample-name ${NormalReplicateId} \\
        ${NetChop} \\
        ${NetMHCstab} \\
        ${anno_vcf} ${TumorReplicateId} ${hla_type} ${params.baff_tools} ./${TumorReplicateId}_${hla_type}
    """
}

process concat_mhcI_files {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/MCH_Class_I/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        '*filtered.tsv',
        '*.all_epitopes.tsv'
    ) from mhcI_out_f.groupTuple(by: 0)

    output:
    set(
        TumorReplicateId,
        file("*_MHCI_filtered.tsv")
    ) into (MHCI_final_ranked, MHCI_final_immunogenicity)
    file("${TumorReplicateId}_MHCI_all_epitopes.tsv")

    script:
    """
    sed -e '2,\${/^Chromosome/d' -e '}' *filtered.tsv > ${TumorReplicateId}_MHCI_filtered.tsv
    sed -e '2,\${/^Chromosome/d' -e '}' *.all_epitopes.tsv > ${TumorReplicateId}_MHCI_all_epitopes.tsv
    """
}

process concat_mhcII_files {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/MCH_Class_II/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        '*filtered.tsv',
        '*.all_epitopes.tsv'
    ) from mhcII_out_f.groupTuple(by: 0)

    output:
    set(
        TumorReplicateId,
        file("*_MHCII_filtered.tsv")
    ) into MHCII_final_ranked
    file("${TumorReplicateId}_MHCII_all_epitopes.tsv")

    script:
    """
    sed -e '2,\${/^Chromosome/d' -e '}' *filtered.tsv > ${TumorReplicateId}_MHCII_filtered.tsv
    sed -e '2,\${/^Chromosome/d' -e '}' *.all_epitopes.tsv > ${TumorReplicateId}_MHCII_all_epitopes.tsv
    """
}


process ranked_reports {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/",
        mode: params.publishDirMode

    input:
    set (
        TumorReplicateId,
        pvacseq_mhcI_file,
        pvacseq_mhcII_file
    ) from MHCI_final_ranked
        .combine(MHCII_final_ranked, by: 0)
    

    output:
    file("**/*_MHCI_filtered.condensed.ranked.tsv")
    file("**/*_MHCII_filtered.condensed.ranked.tsv")

    script:
    """
    mkdir ./MCH_Class_I/
    pvacseq generate_condensed_ranked_report \\
        -m lowest \\
        $pvacseq_mhcI_file \\
        ./MCH_Class_I/${TumorReplicateId}_MHCI_filtered.condensed.ranked.tsv
    mkdir ./MCH_Class_II/
    pvacseq generate_condensed_ranked_report \\
        -m lowest $pvacseq_mhcII_file \\
        ./MCH_Class_II/${TumorReplicateId}_MHCII_filtered.condensed.ranked.tsv
    """
}

// mixMHC2pred
// boilerplate code to be finished (Didi)
// /*

// process 'pVACtools_generate_protein_seq' {
//     tag "$TumorReplicateId"

//     // publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/",
//     //     mode: params.publishDirMode

//     input:
//     set(
//         TumorReplicateId,
//         NormalReplicateId,
//         _,
//         _,
//         _,
//         _,
//         _,
//         anno_vcf,
//         anno_vcf_tbi
//     ) from mkPhasedVCF_out_pVACseq_ch1

//     output:
//     set(
//         TumorReplicateId,
//         file("**/MHC_Class_I/*.filtered.tsv"),
//         file("**/MHC_Class_I/*.all_epitopes.tsv")
//     ) optional true into mhcI_out_f
//     set(
//         TumorReplicateId,
//         file("**/MHC_Class_II/*.filtered.tsv"),
//         file("**/MHC_Class_II/*.all_epitopes.tsv")
//     ) optional true into mhcII_out_f

//     script:
//     hla_type = (hla_types - ~/\n/)
//     NetChop = params.use_NetChop ? "--net-chop-method cterm" : ""
//     NetMHCstab = params.use_NetMHCstab ? "--netmhc-stab" : ""
//     """
//     pvacseq run \\
//         --iedb-install-directory /opt/iedb \\
//         -t ${task.cpus} \\
//         -p ${vep_phased_vcf_gz} \\
//         -e ${params.epitope_len} \\
//         --normal-sample-name ${NormalReplicateId} \\
//         ${NetChop} \\
//         ${NetMHCstab} \\
//         ${anno_vcf} ${TumorReplicateId} ${hla_type} ${params.baff_tools} ./${TumorReplicateId}_${hla_type}
//     """
// }

/*
process mixMHC2pred {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/12_mixMHC2pred",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        file(mutatedProtein),
        alleles // to be implemented
    ) from mkPhasedVCF_out_mixMHCpred_ch0
        .combine(hlas_ch2.splitText(), by: 0) // to be created

    output:
    set(
        TumorReplicateId,
        file("${TumorReplicateId}_mixMHC2pred.tsv"),
    )

    script:
    """
    pepChopper.py -i ${mutatedProtein} -o chopped_peps.txt
    MixMHC2pred -i chopped_peps.txt -o ${TumorReplicateId}_mixMHC2pred.tsv -a ${alleles}
    """
}
*/

// process immunogenicity_scoring {
//     tag "$TumorReplicateId"

//     publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/MCH_Class_I/",
//         mode: params.publishDirMode

//     input:
// 	set (
// 		TumorReplicateId,
// 		mhCI_tag_immunogenicity
// 	) from MHCI_final_immunogenicity
//     // val(TumorReplicateId) from mhCI_tag_immunogenicity
//     // file pvacseq_file from MHCI_final_immunogenicity

//     output:
//     file("*_immunogenicity.tsv") optional true

//     script:
//     """
//     get_epitopes.py \
//         --pvacseq_out $pvacseq_file \\
//         --sample_id $TumorReplicateId \\
//         --output ./${TumorReplicateId}_epitopes.tsv
//     NeoAg_immunogenicity_predicition_GBM.R \\
//         ./${TumorReplicateId}_epitopes.tsv ./${TumorReplicateId}_immunogenicity.tsv
//     """
// }

if(params.TCR) {
    process mixcr_DNA_tumor {
        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/13_TCRs",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            readFWD,
            readREV,
            _
        ) from reads_tumor_mixcr_DNA_ch

        output:
        set(
            TumorReplicateId,
            file("${TumorReplicateId}_mixcr_DNA.*.txt"),
        )

        script:
        reads = (single_end) ? readFWD : readFWD + " " + readREV
        """
        $MIXCR analyze shotgun \\
            --align "--threads ${task.cpus}" \\
            --species hs \\
            --starting-material dna \\
            --only-productive \\
            $reads \\
            ${TumorReplicateId}_mixcr_DNA
        """
    }

    process mixcr_DNA_normal {
        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/13_TCRs",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            readFWD,
            readREV,
            _
        ) from reads_normal_mixcr_DNA_ch

        output:
        set(
            NormalReplicateId,
            file("${NormalReplicateId}_mixcr_DNA.*.txt"),
        )

        script:
        reads = (single_end) ? readFWD : readFWD + " " + readREV
        """
        $MIXCR analyze shotgun \\
            --align "--threads ${task.cpus}" \\
            --species hs \\
            --starting-material dna \\
            --only-productive \\
            $reads \\
            ${NormalReplicateId}_mixcr_DNA
        """
    }

    process mixcr_RNA {
        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/13_TCRs",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            readRNAFWD,
            readRNAREV
        ) from reads_tumor_mixcr_RNA_ch

        output:
        set(
            TumorReplicateId,
            file("${TumorReplicateId}_mixcr_RNA.*.txt"),
        )

        script:
        readsRNA = (single_end_RNA) ? readRNAFWD : readRNAFWD + " " + readRNAREV
        """
        $MIXCR analyze shotgun \\
            --align "--threads ${task.cpus}" \\
            --species hs \\
            --starting-material rna \\
            --only-productive \\
            $readsRNA \\
            ${TumorReplicateId}_mixcr_RNA
        """
    }
}

/*
***********************************
*  Generate final multiQC output  *
***********************************
*/
process multiQC {

    publishDir "${params.outputDir}/$TumorReplicateId/02_QC", mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("*"),
        file("*"),
        file("*"),
        file("*"),
        file("*"),
        file("*"),
        file("*"),
        file("*"),
        file("*"),
        file("*")
    )   from ch_fastqc
            .combine(ch_fastp_tumor, by: [0,1])
            .combine(ch_fastp_normal, by: [0,1])
            .combine(ch_fastqc_trimmed, by: [0,1])
            .combine(ch_fastp_RNAseq, by: [0,1])
            .combine(ch_fastqc_trimmed_RNAseq, by: [0,1])
            .combine(MarkDuplicatesTumor_out_ch4, by: [0,1])
            .combine(alignmentMetricsTumor_ch, by: [0,1])
            .combine(MarkDuplicatesNormal_out_ch3, by: [0,1])
            .combine(alignmentMetricsNormal_ch, by: [0,1])

    output:
    file("multiqc_data/*")
    file("multiqc_report.html")

    script:
    """
    multiqc .
    """

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
    if(params.WES) {
        if (params.references.size() != 15) exit 1, """
        ERROR: Not all References needed found in configuration
        Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
        """
        return [
            'RefFasta'          : checkParamReturnFileReferences("RefFasta"),
            'RefIdx'            : checkParamReturnFileReferences("RefIdx"),
            'RefDict'           : checkParamReturnFileReferences("RefDict"),
            'BwaRef'            : checkParamReturnFileReferences("BwaRef"),
            'VepFasta'          : checkParamReturnFileReferences("VepFasta"),
            'BaitsBed'          : checkParamReturnFileReferences("BaitsBed"),
            'RegionsBed'        : checkParamReturnFileReferences("RegionsBed"),
            'YaraIndex'        	: checkParamReturnFileReferences("YaraIndex"),
            'HLAHDFreqData'     : checkParamReturnFileReferences("HLAHDFreqData"),
            'HLAHDGeneSplit'    : checkParamReturnFileReferences("HLAHDGeneSplit"),
            'HLAHDDict'       	: checkParamReturnFileReferences("HLAHDDict"),
            'STARidx'       	: checkParamReturnFileReferences("STARidx"),
            'AnnoFile'        	: checkParamReturnFileReferences("AnnoFile"),
            'RefFASTA'        	: checkParamReturnFileReferences("RefFASTA"),
            'RefIDX'        	: checkParamReturnFileReferences("RefIDX")
        ]
    } else {
        if (params.references.size() != 15) exit 1, """
        ERROR: Not all References needed found in configuration
        Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
        """
        return [
            'RefFasta'          : checkParamReturnFileReferences("RefFasta"),
            'RefIdx'            : checkParamReturnFileReferences("RefIdx"),
            'RefDict'           : checkParamReturnFileReferences("RefDict"),
            'BwaRef'            : checkParamReturnFileReferences("BwaRef"),
            'VepFasta'          : checkParamReturnFileReferences("VepFasta"),
            'YaraIndex'        	: checkParamReturnFileReferences("YaraIndex"),
            'HLAHDFreqData'     : checkParamReturnFileReferences("HLAHDFreqData"),
            'HLAHDGeneSplit'    : checkParamReturnFileReferences("HLAHDGeneSplit"),
            'HLAHDDict'       	: checkParamReturnFileReferences("HLAHDDict"),
            'STARidx'       	: checkParamReturnFileReferences("STARidx"),
            'AnnoFile'        	: checkParamReturnFileReferences("AnnoFile"),
            'RefFASTA'        	: checkParamReturnFileReferences("RefFASTA"),
            'RefIDX'        	: checkParamReturnFileReferences("RefIDX")
        ]
    }
}

def defineDatabases() {
    if (params.databases.size() != 16) exit 1, """
    ERROR: Not all Databases needed found in configuration
    Please check if Mills_and_1000G_gold_standard, CosmicCodingMuts, DBSNP, GnomAD, and knownIndels are given.
    """
    return [
        'MillsGold'      : checkParamReturnFileDatabases("MillsGold"),
        'MillsGoldIdx'   : checkParamReturnFileDatabases("MillsGoldIdx"),
        'hcSNPS1000G'    : checkParamReturnFileDatabases("hcSNPS1000G"),
        'hcSNPS1000GIdx' : checkParamReturnFileDatabases("hcSNPS1000GIdx"),
        'HapMap'         : checkParamReturnFileDatabases("HapMap"),
        'HapMapIdx'      : checkParamReturnFileDatabases("HapMapIdx"),
        'Cosmic'         : checkParamReturnFileDatabases("Cosmic"),
        'CosmicIdx'      : checkParamReturnFileDatabases("CosmicIdx"),
        'DBSNP'          : checkParamReturnFileDatabases("DBSNP"),
        'DBSNPIdx'       : checkParamReturnFileDatabases("DBSNPIdx"),
        'GnomAD'         : checkParamReturnFileDatabases("GnomAD"),
        'GnomADIdx'      : checkParamReturnFileDatabases("GnomADIdx"),
        'GnomADfull'     : checkParamReturnFileDatabases("GnomADfull"),
        'GnomADfullIdx'  : checkParamReturnFileDatabases("GnomADfullIdx"),
        'KnownIndels'    : checkParamReturnFileDatabases("KnownIndels"),
        'KnownIndelsIdx' : checkParamReturnFileDatabases("KnownIndelsIdx")
    ]
}

def helpMessage() {
    log.info ""
    log.info "----------------------------"
    log.info "--        U S A G E       --"
    log.info "----------------------------"
    log.info ""
    log.info ' nextflow run wes.nf "--readsTumor|--batchFile" "[--readsNormal]" "--RegionsBed" "--BaitsBed" [--single_end]'
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
    log.info "tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,group"
    log.info "tumor1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,normal1,Normal1_reads_1.fastq,Normal1_reads_2.fastq,group1"
    log.info "tumor2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,normal2,Normal2_reads_1.fastq,Normal2_reads_2.fastq,group1"
    log.info "..."
    log.info "sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,NormalN_reads_1.fastq,NormalN_reads_2.fastq,groupX"
    log.info ""
    log.info "CSV-file, single-end T only reads:"
    log.info "tumorSampleName,readsTumorFWD,readsTumorREV,readsNormalFWD,readsNormalREV,group"
    log.info "tumor1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,noname,NO_FILE,NO_FILE,group1"
    log.info "tumor2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,noname,NO_FILE,NO_FILE,group1"
    log.info "..."
    log.info "sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,NO_FILE,,NO_FILE,groupX"
    log.info ""
    log.info "FASTQ files (can be zipped), if single-end reads are used put NO_FILE instead of *_reads_2.fastq in the REV fields"
    log.info ""
    log.info "--RegionsBed \t\t regions.bed \t\t\t regions.bed file for Exon targets"
    log.info "--BaitsBed \t\t baits.bed \t\t\t baits.bed file for Exon baits"
    log.info ""
    log.info " Optional argument:"
    log.info " ------------------"
    log.info "--readsNormal \t\t reads_{1,2}.fastq \t\t paired-end reads; FASTQ file (can be zipped)"
    log.info "--tumorSampleName \t\t  tumor sample name. If not specified samples will be named according to the fastq filenames."
    log.info "--normalSampleName \t\t  normal sample name. If not specified samples will be named according to the fastq filenames."
    log.info "--trim_adapters \t\t  If true Illumina universal adpter (AGATCGGAAGAG) will be trimmed from reads unless adapter seqs are provided."
    log.info "--adapterSeq \t\t  String of atapter sequence (see --trim_adapers)."
    log.info "--adapterSeqFile \t\t  Fasta file with atapter sequences (see --trim_adapers)."
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
    log.info " hcSNPS1000G/Idx \t\t 1000G_phase1.snps.high_confidence.hg38.vcf.gz/idx \t high confidence SNPS from 1000G project, VCF file, IDX File"
    log.info " HapMap/Idx \t\t hapmap_3.3.hg38.vcf.gz/idx \t HapMap for germline filtration of HaploType caller, VCF file, IDX File"
    log.info " Cosmic/Idx \t\t CosmicCodingMuts.vcf \t\t\t\t Cosmic conding mutations, VCF file, IDX file"
    log.info " DBSNP/Idx \t\t Homo_sapiens_assembly.dbsnp.vcf/idx \t\t SNPS, microsatellites, and small-scale insertions and deletions, VCF file, IDX file"
    log.info " GnomAD/Idx \t\t small_exac_common_3.vcf/idx \t\t\t exonix sites only for contamination estimation from GATK, VCF file, IDX file"
    log.info " GnomADfull/Idx \t\t af-only-gnomad.hg38.vcf.gz/idx \t\t\t for mutect2, VCF file, IDX file"
    log.info " KnownIdenls/Idx \t Homo_sapiens_assembly.known_indels.vcf/idx \t Known Indels from GATK resource Bundle, VCF file, IDX file"
    log.info ""
    log.info " Required software:"
    log.info " ------------------"
    log.info " JAVA7 \t\t\t Version 1.7"
    log.info " JAVA8 \t\t\t Version 1.8"
    log.info " BWA \t\t\t Version 0.7.17"
    log.info " SAMTOOLS \t\t Version 1.9"
    log.info " PICARD \t\t\t Version 2.21.4"
    log.info " GATK3 \t\t\t Version 3.8-0"
    log.info " GATK4 \t\t\t Version 4.1.4.1"
    log.info " VARSCAN \t\t Version 2.4.3"
    log.info " MUTECT1 \t\t Version 1.1.7"
    log.info " BAMREADCOUNT \t\t Version 0.8.0"
    log.info " VEP \t\t\t Version 2.0"
    log.info "-------------------------------------------------------------------------"
}

// workflow complete
workflow.onComplete {
    // Set up the e-mail variables
    def subject = "[icbi/wes] Successful: $workflow.runName"
    if(!workflow.success){
        subject = "[icbi/wes] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
            if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[icbi/wes] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << email_txt
            log.info "[icbi/wes] Sent summary e-mail to $params.email (mail)"
        }
    }

  // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outputDir}/Documentation/" )
    if( !output_d.exists() ) {
        output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[icbi/wes] Pipeline Complete! You can find your results in $baseDir/${params.outputDir}"
}
