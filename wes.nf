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
if  (params.readsTumor != "NO_FILE") log.info " Reads Tumor: \t\t " + params.readsTumor
if  (params.readsNormal != "NO_FILE") log.info " Reads Normal: \t\t " + params.readsNormal
if  (params.readsRNAseq != "NO_FILE") log.info " Reads RNAseq: \t\t " + params.readsRNAseq
if  (params.customHLA != "NO_FILE") log.info " Custom HLA file: \t\t " + params.customHLA
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

// initialize RNAseq presence
have_RNAseq = false

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
if(params.customHLA != "NO_FILE") summary['Custom HLA file']   = params.customHLA
summary['Gender']                        = params.gender
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
summary['VarScan min_map_q']             = params.min_map_q
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

        if (params.readsRNAseq) {
            Channel
                    .fromFilePairs(params.readsRNAseq)
                    .map { reads -> tuple(tumorSampleName, reads[1][0], reads[1][1], "None") }
                    .into { raw_reads_tumor_neofuse_ch; fastqc_readsRNAseq_ch }
            have_RNAseq = true
        } else {
            Channel
                .empty()
                .into { raw_reads_tumor_neofuse_ch; fastqc_readsRNAseq_ch }
            Channel
                .of(tuple(
                    tumorSampleName,
                    "NO_FILE"
                ))
                .set { optitype_RNA_output }

            have_RNAseq = false
        }
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
    // user supplied HLA types (default: NO_FILE, will be checked in get_vhla)
    Channel
            .of(tuple(
                    tumorSampleName,
                    file(params.customHLA)
                ))
            .set { custom_hlas_ch }

    genderMap = [:]

    if (params.gender in ["XX", "XY", "Female", "Male"]) {
        genderMap[row.tumorSampleName] = params.gender
    } else {
        exit 1, "Gender should be one of: XX, XY, Female, Male, got: " + params.gender
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

    pe_dna_count = 0
    se_dna_count = 0
    pe_rna_count = 0
    se_rna_count = 0

    genderMap = [:]

    for ( row in batchCSV ) {
        if(row.gender && row.gender != "None") {
            if (params.gender in ["XX", "XY", "Female", "Male"]) {
               genderMap[row.tumorSampleName] = (row.gender == "Female" || row.gender == "XX") ? "XX" : "XY"
            } else {
                exit 1, "Gender should be one of: XX, XY, Female, Male, got: " + params.gender
            }
        } else {
            println("WARNING: gender not specified assuming: XY")
            genderMap[row.tumorSampleName] = "XY"
        }

        if (row.readsTumorREV == "None") {
            single_end = true
            se_dna_count++
        } else {
            pe_dna_count++
        }

        if (row.readsNormalFWD == "None") {
            exit 1, "No normal sample defined for " + row.readsTumorFWD
        }

        if (row.readsRNAseqFWD == "None") {
            have_RNAseq = false
        } else {
            have_RNAseq = true
            if (row.readsRNAseqREV == "None") {
                single_end_RNA = true
                se_rna_count++
            } else {
                pe_rna_count++
            }
        }
    }

    if (pe_dna_count != 0 && se_dna_count != 0) {
        exit 1, "Please do not mix pe and se DNA read samples in batch file. Create a separate batch file for se and pe DNA samples"
    }

    if (pe_rna_count != 0 && se_rna_count != 0) {
        exit 1, "Please do not mix pe and se RNA read samples in batch file. Create a separate batch file for se and pe RNA samples"
    }

    Channel
            .fromPath(params.batchFile)
            .splitCsv(header:true)
            .map { row -> tuple(row.tumorSampleName,
                                row.normalSampleName,
                                file(row.readsTumorFWD),
                                file((row.readsTumorREV == "None") ? "NO_FILE_REV_T" : row.readsTumorREV),
                                row.group) }
            .into { raw_reads_tumor_ch;
                    fastqc_reads_tumor_ch }

    Channel
            .fromPath(params.batchFile)
            .splitCsv(header:true)
            .map { row -> tuple(row.tumorSampleName,
                                row.normalSampleName,
                                file(row.readsNormalFWD),
                                file((row.readsNormalREV == "None") ? "NO_FILE_REV_N" : row.readsNormalREV),
                                row.group) }
            .into { raw_reads_normal_ch; fastqc_reads_normal_ch }

    if (have_RNAseq) {
        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.tumorSampleName,
                                    row.normalSampleName,
                                    file(row.readsRNAseqFWD),
                                    file(row.readsRNAseqREV)) }
                .into { raw_reads_tumor_neofuse_ch; fastqc_readsRNAseq_ch }
    } else {
        Channel
            .empty()
            .into { raw_reads_tumor_neofuse_ch; fastqc_readsRNAseq_ch }

        Channel
                .fromPath(params.batchFile)
                .splitCsv(header:true)
                .map { row -> tuple(row.tumorSampleName,
                                    file("NO_FILE")) }
                .set { optitype_RNA_output }
    }

    // user supplied HLA types (default: NO_FILE, will be checked in get_vhla)
    Channel
            .fromPath(params.batchFile)
            .splitCsv(header:true)
            .map { row -> tuple(row.tumorSampleName,
                                file((row.HLAfile == "None") ? "NO_FILE" : row.HLAfile)) }
            .set { custom_hlas_ch }
}

reference = defineReference()
database = defineDatabases()
mkTmpDir()

scatter_count = Channel.from(params.scatter_count)
padding = params.readLength + 100

FASTQC        = file(params.FASTQC)
FASTP         = file(params.FASTP)
BWA           = file(params.BWA)
VARSCAN       = file(params.VARSCAN)
GATK4         = file(params.GATK4)
GATK3         = file(params.GATK3)
MUTECT1       = file(params.MUTECT1)
SAMTOOLS      = file(params.SAMTOOLS)
SAMBAMBA      = file(params.SAMBAMBA)
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
MiXMHC2PRED   = file(params.MiXMHC2PRED)
ALLELECOUNT   = file(params.ALLELECOUNT)
FREEC         = file(params.FREEC)
RSCRIPT       = file(params.RSCRIPT)

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
        $JAVA8 ${params.JAVA_Xmx} -jar $PICARD ScatterIntervalsByNs \\
            REFERENCE=$RefFasta \\
            OUTPUT_TYPE=ACGT \\
            OUTPUT=${interval_list.baseName}_merged_padded.interval_list
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
        SplitIntervals_out_scatterBaseRecalTumorGATK4_ch,
        SplitIntervals_out_scatterTumorGATK4applyBQSRS_ch,
        SplitIntervals_out_scatterBaseRecalNormalGATK4_ch,
        SplitIntervals_out_scatterNormalGATK4applyBQSRS_ch,
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
        ScatteredIntervalListToBed_out_ch0,
        ScatteredIntervalListToBed_out_ch1
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

    if (have_RNAseq) {
        if (!single_end && !single_end_RNA) {
            cpus = 6
        } else if (!single_end && single_end_RNA) {
            cpus = 5
        } else if (single_end && !single_end_RNA) {
            cpus = 4
        } else if (single_end && single_end_RNA) {
            cpus = 3
        }
    } else {
        if (!single_end) {
            cpus = 4
        } else {
            cpus = 2
        }
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

    tumor_readsFWD_ext = (tumor_readsFWD_ext == "gz") ? "fastq.gz" : tumor_readsFWD_ext
    normal_readsFWD_ext = (normal_readsFWD_ext == "gz") ? "fastq.gz" : normal_readsFWD_ext

    if (! single_end) {
        tumor_readsREV_simpleName = tumor_readsREV.getSimpleName()
        normal_readsREV_simpleName = normal_readsREV.getSimpleName()
        tumor_readsREV_ext = tumor_readsREV.getExtension()
        normal_readsREV_ext = normal_readsREV.getExtension()

        tumor_readsREV_ext = (tumor_readsREV_ext == "gz") ? "fastq.gz" : tumor_readsREV_ext
        normal_readsREV_ext = (normal_readsREV_ext == "gz") ? "fastq.gz" : normal_readsREV_ext
    }


    if (have_RNAseq) {
        readsRNAseq_FWD_simpleName = readsRNAseq_FWD.getSimpleName()
        readsRNAseq_FWD_ext = readsRNAseq_FWD.getExtension()
        readsRNAseq_FWD_ext = (readsRNAseq_FWD_ext == "gz") ? "fastq.gz" : readsRNAseq_FWD_ext

        if (! single_end_RNA) {
            readsRNAseq_REV_simpleName = readsRNAseq_REV.getSimpleName()
            readsRNAseq_REV_ext = readsRNAseq_REV.getExtension()

            readsRNAseq_REV_ext = (readsRNAseq_REV_ext == "gz") ? "fastq.gz" : readsRNAseq_REV_ext
        }
    }

    if (single_end && single_end_RNA && have_RNAseq)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}_R1.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}_R1.${normal_readsFWD_ext}
        ln -s $readsRNAseq_FWD ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}

        fastqc --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext} \\
            ${TumorReplicateId}_RNA_R1.${readsRNAseq_FWD_ext}
        """
    else if (!single_end && !single_end_RNA && have_RNAseq)
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
    else if (!single_end && single_end_RNA && have_RNAseq)
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
    else if (single_end && !single_end_RNA && have_RNAseq)
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
    else if (single_end && !have_RNAseq)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}_R1.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}_R1.${normal_readsFWD_ext}

        fastqc --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext}
        """
    else if (!single_end && !have_RNAseq)
        """
        ln -s $tumor_readsFWD ${TumorReplicateId}_R1.${tumor_readsFWD_ext}
        ln -s $normal_readsFWD ${NormalReplicateId}_R1.${normal_readsFWD_ext}
        ln -s $tumor_readsREV ${TumorReplicateId}_R2.${tumor_readsREV_ext}
        ln -s $normal_readsREV ${NormalReplicateId}_R2.${normal_readsREV_ext}

        $FASTQC --quiet --threads ${task.cpus} \\
            ${TumorReplicateId}_R1.${tumor_readsFWD_ext} ${TumorReplicateId}_R2.${tumor_readsREV_ext} \\
            ${NormalReplicateId}_R1.${normal_readsFWD_ext} ${NormalReplicateId}_R2.${normal_readsREV_ext}
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
            reads_tumor_uBAM_ch,
            reads_tumor_mixcr_DNA_ch,
            fastqc_reads_tumor_trimmed_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*.json")
        ) into ch_fastp_tumor // multiQC


        script:
        trimmedReads_2 = (single_end) ? "NO_FILE_T" : TumorReplicateId + "_trimmed_R2.fastq.gz"

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
            touch NO_FILE_T
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
            reads_normal_uBAM_ch,
            reads_normal_mixcr_DNA_ch,
            fastqc_reads_normal_trimmed_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*.json")
        ) into ch_fastp_normal // multiQC


        script:
        trimmedReads_2 = (single_end) ? "NO_FILE_N" : NormalReplicateId + "_trimmed_R2.fastq.gz"

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
            touch NO_FILE_N
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
    (
        reads_tumor_ch,
        reads_tumor_uBAM_ch,
        reads_tumor_mixcr_DNA_ch,
        ch_fastqc_trimmed
    ) = raw_reads_tumor_ch.into(4)

    (
        reads_normal_ch,
        reads_normal_uBAM_ch,
        reads_normal_mixcr_DNA_ch
    ) = raw_reads_normal_ch.into(3)

    ch_fastqc_trimmed = ch_fastqc_trimmed
                            .map{ it -> tuple(it[0], it[1], "")}

    (ch_fastp_normal, ch_fastp_tumor, ch_fastqc_trimmed) = ch_fastqc_trimmed.into(3)
}

// adapter trimming RNAseq
if (params.trim_adapters_RNAseq && have_RNAseq) {
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
            reads_tumor_optitype_ch,
            reads_tumor_mixcr_RNA_ch,
            fastqc_readsRNAseq_trimmed_ch
        )
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("*.json")
        ) into ch_fastp_RNAseq // multiQC

        script:
        trimmedReads_2 = (single_end_RNA) ? "NO_FILE" : TumorReplicateId + "_RNA_trimmed_R2.fastq.gz"


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

        if(single_end_RNA)
            """
            $FASTP --thread ${task.cpus} \\
                --in1 ${readsRNAseq_FWD} \\
                --out1 ${TumorReplicateId}_RNA_trimmed_R1.fastq.gz \\
                --json ${TumorReplicateId}_RNA_fastp.json \\
                ${fastpAdapter} \\
                ${params.fastpOpts}
            touch NO_FILE
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
            fastqc --quiet --threads ${task.cpus} ${readsRNAseq_FWD}
            """
        else
            """
            $FASTQC --quiet --threads ${task.cpus} \\
                ${readsRNAseq_FWD} ${readsRNAseq_REV}
            """
    }

} else { // no adapter trimming for RNAseq

    (
        reads_tumor_neofuse_ch,
        reads_tumor_optitype_ch,
        reads_tumor_mixcr_RNA_ch,
        ch_fastp_RNAseq
    ) = raw_reads_tumor_neofuse_ch.into(4)

    ch_fastp_RNAseq = ch_fastp_RNAseq
                            .map{ it -> tuple(it[0], it[1], "")}

    (ch_fastqc_trimmed_RNAseq, ch_fastp_RNAseq) = ch_fastp_RNAseq.into(2)

}


// mix tumor normal channels and add sampleType (T/N) so that we can split again
reads_uBAM_ch = Channel.empty()
                        .mix(
                            reads_tumor_uBAM_ch
                                .combine(Channel.of("T")),
                            reads_normal_uBAM_ch
                                .combine(Channel.of("N"))
                        )


/// start processing reads
process 'make_uBAM' {
// make uBAM

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
        sampleType
    ) from reads_uBAM_ch


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(outFileName)
    ) into uBAM_out_ch0

    script:
    outFileName = (sampleType == "T") ? TumorReplicateId + "_unaligned.bam" : NormalReplicateId + "_unaligned.bam"
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    if (single_end)
        """
        mkdir -p ${params.tmpDir}
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} FastqToSam \\
            TMP_DIR=${params.tmpDir} \\
            F1=${readsFWD} \\
            READ_GROUP_NAME=${procSampleName} \\
            SAMPLE_NAME=${procSampleName} \\
            LIBRARY_NAME=${procSampleName} \\
            PLATFORM=ILLUMINA \\
            O=${procSampleName}_unaligned.bam
        """
    else
        """
        mkdir -p ${params.tmpDir}
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} FastqToSam \\
            TMP_DIR=${params.tmpDir} \\
            F1=${readsFWD} \\
            F2=${readsREV} \\
            READ_GROUP_NAME=${procSampleName} \\
            SAMPLE_NAME=${procSampleName} \\
            LIBRARY_NAME=${procSampleName} \\
            PLATFORM=ILLUMINA \\
            O=${procSampleName}_unaligned.bam
        """
}

reads_ch = Channel.empty()
                    .mix(
                        reads_tumor_ch
                            .combine(Channel.of("T")),
                        reads_normal_ch
                            .combine(Channel.of("N"))
                    )


process 'Bwa' {
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
        sampleType
    ) from reads_ch
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
        sampleType,
        file(outFileName)
    ) into Bwa_out_ch0

    script:
    outFileName = (sampleType == "T") ? TumorReplicateId + "_aligned.bam" : NormalReplicateId + "_aligned.bam"
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId

    sort_threads = (task.cpus.compareTo(8) == 1) ? 8 : task.cpus
    if (single_end)
        """
        $BWA mem \\
            -R "@RG\\tID:${procSampleName}\\tLB:${procSampleName}\\tSM:${procSampleName}\\tPL:ILLUMINA" \\
            -M ${RefFasta} \\
            -t ${task.cpus} \\
            -Y \\
            ${readsFWD} | \\
        $SAMTOOLS view -@2 -Shbu - | \\
        $SAMBAMBA sort \\
            --sort-picard \\
            --tmpdir=${params.tmpDir} \\
            -m ${params.SB_sort_mem} \\
            -l 6 \\
            -t ${sort_threads} \\
            -o ${outFileName} \\
            /dev/stdin
        """
    else
        """
        $BWA mem \\
            -R "@RG\\tID:${procSampleName}\\tLB:${procSampleName}\\tSM:${procSampleName}\\tPL:ILLUMINA" \\
            -M ${RefFasta} \\
            -t ${task.cpus} \\
            -Y \\
            ${readsFWD} \\
            ${readsREV}  | \\
        $SAMTOOLS view -@2 -Shbu - | \\
        $SAMBAMBA sort \\
            --sort-picard \\
            --tmpdir=${params.tmpDir} \\
            -m ${params.SB_sort_mem} \\
            -l 6 \\
            -t ${sort_threads} \\
            -o ${outFileName} \\
            /dev/stdin
        """
}


process 'merge_uBAM_BAM' {
// merge alinged BAM and uBAM

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(BAM),
        file(uBAM)
    ) from Bwa_out_ch0
        .combine(uBAM_out_ch0, by: [0,1,2])

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
        sampleType,
        file("${procSampleName}_merged.bam")
    ) into uBAM_BAM_out_ch

    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId

    paired_run = (single_end) ? 'false' : 'true'
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} MergeBamAlignment \\
        TMP_DIR=${params.tmpDir} \\
        VALIDATION_STRINGENCY=SILENT \\
        EXPECTED_ORIENTATIONS=FR \\
        ATTRIBUTES_TO_RETAIN=X0 \\
        REFERENCE_SEQUENCE=${RefFasta} \\
        PAIRED_RUN=${paired_run} \\
        SORT_ORDER="queryname" \\
        IS_BISULFITE_SEQUENCE=false \\
        ALIGNED_READS_ONLY=false \\
        CLIP_ADAPTERS=false \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRamMerge} \\
        ADD_MATE_CIGAR=true \\
        MAX_INSERTIONS_OR_DELETIONS=-1 \\
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \\
        ALIGNER_PROPER_PAIR_FLAGS=true \\
        UNMAP_CONTAMINANT_READS=true \\
        ALIGNED_BAM=${BAM} \\
        UNMAPPED_BAM=${uBAM} \\
        OUTPUT=${procSampleName}_merged.bam
    """
}


process 'MarkDuplicates' {
// Mark duplicates with sambamba

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/01_preprocessing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bam)
    ) from uBAM_BAM_out_ch // BwaTumor_out_ch0

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
        sampleType,
        file("${procSampleName}_aligned_sort_mkdp.bam"),
        file("${procSampleName}_aligned_sort_mkdp.bai")
    ) into (
        MarkDuplicates_out_ch0,
        MarkDuplicates_out_ch1,
        MarkDuplicates_out_ch2,
        MarkDuplicates_out_ch3
    )

    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId

    """
    mkdir -p ${params.tmpDir}
    $SAMBAMBA markdup \\
        -t ${task.cpus} \\
        --tmpdir ${params.tmpDir} \\
        --hash-table-size=${params.SB_hash_table_size } \\
        --overflow-list-size=${params.SB_overflow_list_size} \\
        --io-buffer-size=${params.SB_io_buffer_size} \\
        ${bam} \\
        /dev/stdout | \\
    samtools sort \\
        -@${task.cpus} \\
        -m ${params.STperThreadMem} \\
        -O BAM \\
        -l 0 \\
        /dev/stdin | \\
    $JAVA8 ${params.JAVA_Xmx} -jar $PICARD SetNmMdAndUqTags \\
        TMP_DIR=${params.tmpDir} \\
        R=${RefFasta} \\
        I=/dev/stdin \\
        O=${procSampleName}_aligned_sort_mkdp.bam \\
        CREATE_INDEX=true \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam} \\
        VALIDATION_STRINGENCY=LENIENT

    # samtools index -@${task.cpus} ${procSampleName}_aligned_sort_mkdp.bam
    """
}

// prepare channel for mhc_extract -> hld-hd, optitype
MarkDuplicatesTumor_out_ch0 = MarkDuplicates_out_ch3
                                .filter {
                                    it[2] == "T"
                                }
MarkDuplicatesTumor_out_ch0 = MarkDuplicatesTumor_out_ch0
                                .map{
                                        TumorReplicateId,
                                        NormalReplicateId,
                                        sampleType,
                                        TumorBAM,
                                        TumorBAI  -> tuple(
                                            TumorReplicateId,
                                            NormalReplicateId,
                                            TumorBAM,
                                            TumorBAI
                                        )
                                    }

if(params.WES) {
    process 'alignmentMetrics' {
    // Generate HS metrics using picard

        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/02_QC/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            sampleType,
            file(bam),
            file(bai)
        ) from MarkDuplicates_out_ch0

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
            sampleType,
            file("${procSampleName}.*.txt")
        ) into alignmentMetrics_ch // multiQC

        script:
        procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId

        """
        mkdir -p ${params.tmpDir}
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectHsMetrics \\
            TMP_DIR=${params.tmpDir} \\
            INPUT=${bam} \\
            OUTPUT=${procSampleName}.HS.metrics.txt \\
            R=${RefFasta} \\
            BAIT_INTERVALS=${BaitIntervalsList} \\
            TARGET_INTERVALS=${IntervalsList} \\
            PER_TARGET_COVERAGE=${procSampleName}.perTarget.coverage.txt && \\
        $JAVA8 ${params.JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -jar ${PICARD} CollectAlignmentSummaryMetrics \\
            TMP_DIR=${params.tmpDir} \\
            INPUT=${bam} \\
            OUTPUT=${procSampleName}.AS.metrics.txt \\
            R=${RefFasta} &&
        $SAMTOOLS flagstat -@${task.cpus} ${bam} > ${procSampleName}.flagstat.txt
        """
    }

    alignmentMetricsTumor_ch = Channel.create()
    alignmentMetricsNormal_ch = Channel.create()

    alignmentMetrics_ch
        .choice(
            alignmentMetricsTumor_ch, alignmentMetricsNormal_ch
        ) { it[2] == "T" ? 0 : 1 }

    alignmentMetricsTumor_ch = alignmentMetricsTumor_ch
            .map{
                    TumorReplicateId, NormalReplicateId,
                    sampleTypeN, metricFiles -> tuple(
                        TumorReplicateId,
                        NormalReplicateId,
                        metricFiles
                    )
                }

    alignmentMetricsNormal_ch = alignmentMetricsNormal_ch
            .map{
                    TumorReplicateId, NormalReplicateId,
                    sampleTypeN, metricFiles -> tuple(
                        TumorReplicateId,
                        NormalReplicateId,
                        metricFiles
                    )
                }

} else {
    // bogus channel for multiqc
    alignmentMetricsTumor_ch = MarkDuplicates_out_ch0
                            .map{ it -> tuple(it[0],it[1], "")}
    (alignmentMetricsNormal_ch, alignmentMetricsTumor_ch) = alignmentMetricsTumor_ch.into(2)
}


/*
*********************************************
**             M U T E C T 2               **
*********************************************
*/

process 'scatterBaseRecalGATK4' {
/*
 BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
 Recalibration (BQSR)
 ApplyBQSR (GATK4): apply BQSR table to reads
*/

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bam),
        file(bai),
        file(intervals)
    ) from MarkDuplicates_out_ch1
        .combine(
            SplitIntervals_out_scatterBaseRecalTumorGATK4_ch.flatten()  // TODO: change channel name here and above
        )

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
        sampleType,
        file("${procSampleName}_${intervals}_bqsr.table")
    ) into scatterBaseRecalGATK4_out_ch0


    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    """
    mkdir -p ${params.tmpDir}
    $GATK4 BaseRecalibrator \\
        --java-options '${params.JAVA_Xmx}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_bqsr.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold}
    """
}

process 'gatherGATK4scsatteredBQSRtables' {
// gather scattered bqsr tables

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bqsr_table)
    ) from scatterBaseRecalGATK4_out_ch0
        .groupTuple(by: [0, 1, 2])


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file("${procSampleName}_bqsr.table")
    ) into gatherBQSRtables_out_ch0


    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    """
    mkdir -p ${params.tmpDir}

    $GATK4 GatherBQSRReports \\
        -I ${bqsr_table.join(" -I ")} \\
        -O ${procSampleName}_bqsr.table
    """
}

process 'scatterGATK4applyBQSRS' {
/*
 ApplyBQSR (GATK4): apply BQSR table to reads
*/

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bam),
        file(bai),
        file(bqsr_table),
        file(intervals)
    ) from MarkDuplicates_out_ch2
        .combine(gatherBQSRtables_out_ch0, by: [0, 1, 2])
        .combine(
            SplitIntervals_out_scatterTumorGATK4applyBQSRS_ch.flatten()  // TODO: change channel name here and above
        )

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
        sampleType,
        file("${procSampleName}_${intervals}_recal4.bam"),
        file("${procSampleName}_${intervals}_recal4.bai")
    ) into scatterGATK4applyBQSRS_out_GatherRecalBamFiles_ch0


    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    """
    mkdir -p ${params.tmpDir}
    $GATK4 ApplyBQSR \\
        --java-options '${params.JAVA_Xmx}' \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -R ${RefFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_recal4.bam \\
        --bqsr-recal-file ${bqsr_table}
    """
}

process 'GatherRecalBamFiles' {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_mutect2/processing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bam),
        file(bai)
    ) from scatterGATK4applyBQSRS_out_GatherRecalBamFiles_ch0
        .toSortedList({a, b -> a[3].baseName <=> b[3].baseName})
        .flatten()
        .collate(5)
        .groupTuple(by: [0, 1, 2])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file("${procSampleName}_recal4.bam"),
        file("${procSampleName}_recal4.bam.bai")
    ) into (
        BaseRecalGATK4_out_ch0,
        BaseRecalGATK4_out_ch1,
        BaseRecalGATK4_out_ch2, // into mutect2
        GatherRecalBamFiles_out_IndelRealignerIntervals_ch0
    )

    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 -XX:ParallelGCThreads=${task.cpus} ${params.JAVA_Xmx} -jar $PICARD GatherBamFiles \\
        TMP_DIR=${params.tmpDir} \\
        I=${bam.join(" I=")} \\
        O=/dev/stdout \\
        CREATE_INDEX=false \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam} | \\
    $SAMTOOLS sort \\
        -@${task.cpus} \\
        -m ${params.STperThreadMem} \\
        -o ${procSampleName}_recal4.bam -
    $SAMTOOLS index -@${task.cpus} ${procSampleName}_recal4.bam
    """
}

// make tumor ch BaseRecalTumorGATK4_out_ch3 for CollectSequencingArtifactMetrics
BaseRecalTumorGATK4_out_ch0 = BaseRecalGATK4_out_ch1
                                .filter {
                                    it[2] == "T"
                                }


process 'GetPileup' {
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
        sampleType,
        file(bam),
        file(bai)
    ) from BaseRecalGATK4_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file("${procSampleName}_pileup.table")
    ) into GetPileup_out_ch0

    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    """
    mkdir -p ${params.tmpDir}

    $GATK4 GetPileupSummaries \\
        --tmp-dir ${params.tmpDir} \\
        -I ${bam} \\
        -O ${procSampleName}_pileup.table \\
        -L ${IntervalsList} \\
        --variant ${GnomAD}
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
        sampleType,
        file(bam),
        file(bai)
    ) from BaseRecalTumorGATK4_out_ch0

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

BaseRecalTumor = Channel.create()
BaseRecalNormal = Channel.create()

BaseRecalGATK4_out_ch2
    .choice(
        BaseRecalTumor, BaseRecalNormal
    ) { it[2] == "T" ? 0 : 1 }

(BaseRecalNormal_out_ch0, BaseRecalNormal) = BaseRecalNormal.into(2)
BaseRecalNormal_out_ch0 = BaseRecalNormal_out_ch0
        .map{ TumorReplicateId, NormalReplicateId,
            sampleTypeN, recalNormalBAM, recalNormalBAI -> tuple(
            TumorReplicateId, NormalReplicateId,
            recalNormalBAM, recalNormalBAI
            )}


BaseRecalGATK4_out = BaseRecalTumor.combine(BaseRecalNormal, by: [0,1])
BaseRecalGATK4_out = BaseRecalGATK4_out
    .map{ TumorReplicateId, NormalReplicateId,
          sampleTypeT, recalTumorBAM, recalTumorBAI,
          sampleTypeN, recalNormalBAM, recalNormalBAI -> tuple(
          TumorReplicateId, NormalReplicateId,
          recalTumorBAM, recalTumorBAI,
          recalNormalBAM, recalNormalBAI
        )}


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
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from BaseRecalGATK4_out
        .combine(
            SplitIntervals_out_ch2.flatten()
        )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${intervals}.vcf.gz"),
        file("${TumorReplicateId}_${intervals}.vcf.gz.stats"),
        file("${TumorReplicateId}_${intervals}.vcf.gz.tbi"),
        file("${TumorReplicateId}_${intervals}-f1r2.tar.gz")
    ) into Mutect2_out_ch0


    script:
    if(params.mutect2ponFile != false) {
        panel_of_normals = "--panel-of-normals ${params.mutect2ponFile}"
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
        --f1r2-tar-gz ${TumorReplicateId}_${intervals}-f1r2.tar.gz \\
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
        file(idx),
        file(f1r2_tar_gz)
    ) from Mutect2_out_ch0
        .groupTuple(by: [0, 1])


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz.tbi"),
        file("${TumorReplicateId}_${NormalReplicateId}_mutect2_raw.vcf.gz.stats"),
        file("${TumorReplicateId}_${NormalReplicateId}_read-orientation-model.tar.gz")
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

    $GATK4 LearnReadOrientationModel \\
        --tmp-dir ${params.tmpDir} \\
        -I ${f1r2_tar_gz.join(" -I ")} \\
        -O ${TumorReplicateId}_${NormalReplicateId}_read-orientation-model.tar.gz
    """
}


PileupTumor = Channel.create()
PileupNormal = Channel.create()

GetPileup_out_ch0
    .choice(
        PileupTumor, PileupNormal
    ) { it[2] == "T" ? 0 : 1 }

GetPileup_out = PileupTumor.combine(PileupNormal, by: [0,1])
GetPileup_out = GetPileup_out
    .map{ TumorReplicateId, NormalReplicateId, sampleTypeT, pileupTumor, sampleTypeN, pileupNormal -> tuple(
            TumorReplicateId, NormalReplicateId, pileupTumor, pileupNormal
        )}

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
        file(pileupNormal),
        file(vcf),
        file(vcfIdx),
        file(vcfStats),
        file(f1r2_tar_gz)
    ) from GetPileup_out
        .combine(gatherMutect2VCFs_out_ch0, by :[0,1])

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
            --ob-priors ${f1r2_tar_gz} \\
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
            --ob-priors ${f1r2_tar_gz} \\
            -O ${TumorReplicateId}_${NormalReplicateId}_oncefiltered.vcf.gz && \\
        $GATK4 SelectVariants \\
            --tmp-dir ${params.tmpDir} \\
            --variant ${TumorReplicateId}_${NormalReplicateId}_oncefiltered.vcf.gz \\
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
    ) from BaseRecalNormal_out_ch0
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

    tag "$TumorReplicateId"

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


// END HTC



/*
*********************************************
**             V A R S C A N               **
*********************************************
*/


process 'IndelRealignerIntervals' {
/*
 RealignerTargetCreator (GATK3): define intervals to target for local realignment
 IndelRealigner (GATK3): perform local realignment of reads around indels
*/

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bam),
        file(bai)
    ) from GatherRecalBamFiles_out_IndelRealignerIntervals_ch0

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
        sampleType,
        file("${procSampleName}_recal4_realign_${interval}.bam"),
        file("${procSampleName}_recal4_realign_${interval}.bai")
    ) into IndelRealignerIntervals_out_GatherRealignedBamFiles_ch0

    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
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

process 'GatherRealignedBamFiles' {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/03_varscan/processing/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(bam),
        file(bai)
    ) from IndelRealignerIntervals_out_GatherRealignedBamFiles_ch0
        .toSortedList({a, b -> a[3].baseName <=> b[3].baseName})
        .flatten()
        .collate(5)
        .groupTuple(by: [0,1,2])

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file("${procSampleName}_recal_realign.bam"),
        file("${procSampleName}_recal_realign.bam.bai")
    ) into GatherRealignedBamFiles_out_ch

    script:
    procSampleName = (sampleType == "T") ? TumorReplicateId : NormalReplicateId
    """
    mkdir -p ${params.tmpDir}

    $JAVA8 -XX:ParallelGCThreads=${task.cpus} ${params.JAVA_Xmx} -jar $PICARD GatherBamFiles \\
        TMP_DIR=${params.tmpDir} \\
        I=${bam.join(" I=")} \\
        O=/dev/stdout \\
        CREATE_INDEX=false \\
        MAX_RECORDS_IN_RAM=${params.maxRecordsInRam} | \\
    $SAMTOOLS sort \\
        -@${task.cpus} \\
        -m ${params.STperThreadMem} \\
        -o ${procSampleName}_recal_realign.bam -
    $SAMTOOLS index -@${task.cpus} ${procSampleName}_recal_realign.bam
    """
}

recalRealTumor = Channel.create()
recalRealNormal = Channel.create()

(
    GatherRealignedBamFiles_out_AlleleCounter_ch0,
    GatherRealignedBamFiles_out_Mpileup4ControFREEC_ch0,
    GatherRealignedBamFiles_out_ch
) = GatherRealignedBamFiles_out_ch.into(3)

GatherRealignedBamFiles_out_ch
    .choice(
        recalRealTumor, recalRealNormal
    ) { it[2] == "T" ? 0 : 1 }


 (
    recalRealTumor_tmp,
    recalRealTumor
) = recalRealTumor.into(3)

recalRealTumor_tmp = recalRealTumor_tmp
                        .map{ TumorReplicateId, NormalReplicateId,
                                sampleTypeT, realTumorBAM, realTumorBAI -> tuple(
                                    TumorReplicateId, NormalReplicateId,
                                    realTumorBAM, realTumorBAI
                                )
                            }

(
    GatherRealignedBamFilesTumor_out_FilterVarscan_ch0,
    GatherRealignedBamFilesTumor_out_mkPhasedVCF_ch0
) = recalRealTumor_tmp.into(2)


RealignedBamFiles = recalRealTumor.combine(recalRealNormal, by: [0,1])
RealignedBamFiles = RealignedBamFiles
                    .map{ TumorReplicateId, NormalReplicateId,
                            sampleTypeT, realTumorBAM, realTumorBAI,
                            sampleTypeN, realNormalBAM, realNormalBAI -> tuple(
                                TumorReplicateId, NormalReplicateId,
                                realTumorBAM, realTumorBAI,
                                realNormalBAM, realNormalBAI
                            )
                        }
(
    GatherRealignedBamFiles_out_VarscanSomaticScattered_ch0,
    GatherRealignedBamFiles_out_Mutect1scattered_ch0,
    GatherRealignedBamFiles_out_MantaSomaticIndels_ch0,
    GatherRealignedBamFiles_out_StrelkaSomatic_ch0,
    GatherRealignedBamFiles_out_Sequenza_ch0
) = RealignedBamFiles.into(5)


process 'VarscanSomaticScattered' {

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(Tumorbam),
        file(Tumorbai),
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from GatherRealignedBamFiles_out_VarscanSomaticScattered_ch0
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
        file(snpSomatic),
        file(snpSomaticHc),
        file(snpLOH),
        file(snpLOHhc),
        file(snpGerm),
        file(snpGemHc),
        file(indelSomatic),
        file(indelSomaticHc),
        file(indelLOH),
        file(indelLOHhc),
        file(indelGerm),
        file(indelGemHc)
    ) from GatherRealignedBamFilesTumor_out_FilterVarscan_ch0 // BaseRecalTumorGATK3_out_ch1
        .combine(ProcessVarscanSNP_out_ch0, by :[0,1])
        .combine(ProcessVarscanIndel_out_ch0, by :[0,1])

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
        -q${params.min_map_q} \\
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
        -q${params.min_map_q} \\
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
        file(Normalbam),
        file(Normalbai),
        file(intervals)
    ) from GatherRealignedBamFiles_out_Mutect1scattered_ch0
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
if (! single_end) {
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
            file(Normalbam),
            file(Normalbai),
            file(RegionsBedGz),
            file(RegionsBedGzTbi)
        ) from GatherRealignedBamFiles_out_MantaSomaticIndels_ch0
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
} else {
    GatherRealignedBamFiles_out_MantaSomaticIndels_ch0
        .map {  item -> tuple(item[0],
                              item[1],
                              "NO_FILE",
                              "NO_FILE_IDX") }
        .set { MantaSomaticIndels_out_ch0 }
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
        file(Normalbam),
        file(Normalbai),
        file(manta_indel),
        file(manta_indel_tbi),
        file(RegionsBedGz),
        file(RegionsBedGzTbi)
    ) from GatherRealignedBamFiles_out_StrelkaSomatic_ch0
        .combine(MantaSomaticIndels_out_ch0, by:[0,1])
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
    manta_indel_candidates = single_end ? "" : "--indelCandidates ${manta_indel}"
    """
    configureStrelkaSomaticWorkflow.py --tumorBam ${Tumorbam} --normalBam  ${Normalbam} \\
        --referenceFasta ${RefFasta} \\
        ${manta_indel_candidates} \\
        --runDir strelka_${TumorReplicateId} --callRegions ${RegionsBedGz} --exome
    strelka_${TumorReplicateId}/runWorkflow.py -m local -j ${task.cpus}
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
    conda 'assets/py3_icbi.yml'

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
    ) from GatherRealignedBamFilesTumor_out_mkPhasedVCF_ch0 // BaseRecalTumorGATK3_out_ch2
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
        mkPhasedVCF_out_ch0,
        mkPhasedVCF_out_pVACseq_ch0
    )
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf.gz"),
        file("${TumorReplicateId}_${NormalReplicateId}_tumor_vep.vcf.gz.tbi")
    ) into mkPhasedVCF_out_Clonality_ch0
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

// CNVs: ASCAT + FREEC


// adopted from sarek nfcore
process AlleleCounter {

    tag "$TumorReplicateId - $NormalReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(BAM),
        file(BAI)
    ) from GatherRealignedBamFiles_out_AlleleCounter_ch0

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
        file(acLoci)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.acLoci ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(outFileName)
    ) into AlleleCounter_out_ch0


    script:
    outFileName = (sampleType == "T") ? TumorReplicateId + ".alleleCount" : NormalReplicateId + ".alleleCount"
    """
    $ALLELECOUNT \\
        -l ${acLoci} \\
        -r ${RefFasta} \\
        -b ${BAM} \\
        -o ${outFileName}
    """
}

alleleCountOutNormal = Channel.create()
alleleCountOutTumor = Channel.create()

AlleleCounter_out_ch0
    .choice(
        alleleCountOutTumor, alleleCountOutNormal
    ) { it[2] == "T" ? 0 : 1 }

AlleleCounter_out_ch0 = alleleCountOutTumor.combine(alleleCountOutNormal, by: [0,1])


// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process ConvertAlleleCounts {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/14_ASCAT/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        _,
        file(alleleCountTumor),
        _,
        file(alleleCountNormal),
    ) from AlleleCounter_out_ch0


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}.BAF"),
        file("${TumorReplicateId}.LogR"),
        file("${NormalReplicateId}.BAF"),
        file("${NormalReplicateId}.LogR")
    ) into ConvertAlleleCounts_out_ch

    script:
    gender = genderMap[TumorReplicateId]
    """
    Rscript ${workflow.projectDir}/bin/convertAlleleCounts.r \\
        ${TumorReplicateId} ${alleleCountTumor} ${NormalReplicateId} ${alleleCountNormal} ${gender}
    """
}

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process 'Ascat' {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/14_ASCAT/",
        mode: params.publishDirMode


    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(bafTumor),
        file(logrTumor),
        file(bafNormal),
        file(logrNormal)
    ) from ConvertAlleleCounts_out_ch

    file(acLociGC) from Channel.value(reference.acLociGC)

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}.cnvs.txt"),
        file("${TumorReplicateId}.purityploidy.txt"),
    ) into Ascat_out_Clonality_ch0
    file("${TumorReplicateId}.*.{png,txt}")


    script:
    gender = genderMap[TumorReplicateId]
    """
    # get rid of "chr" string if there is any
    for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
    Rscript ${workflow.projectDir}/bin/run_ascat.r ${bafTumor} ${logrTumor} ${bafNormal} ${logrNormal} ${TumorReplicateId} ${baseDir} ${acLociGC} ${gender}
    """
}

process 'Mpileup4ControFREEC' {

    tag "$TumorReplicateId - $NormalReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(BAM),
        file(BAI)
    ) from GatherRealignedBamFiles_out_Mpileup4ControFREEC_ch0

    each file(interval) from ScatteredIntervalListToBed_out_ch1.flatten()


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
        sampleType,
        file("${TumorReplicateId}_${NormalReplicateId}_${sampleType}_${interval}.pileup.gz")
    ) into Mpileup4ControFREEC_out_ch0

    script:
    """
    $SAMTOOLS mpileup \\
        -q 1 \\
        -f ${RefFasta} \\
        -l ${interval} \\
        ${BAM} | \\
    $BGZIP --threads ${task.cpus} -c > ${TumorReplicateId}_${NormalReplicateId}_${sampleType}_${interval}.pileup.gz
    """


}

 Mpileup4ControFREEC_out_ch0 = Mpileup4ControFREEC_out_ch0.groupTuple(by:[0, 1, 2])

process 'gatherMpileups' {
// Merge scattered Varscan vcfs

    tag "$TumorReplicateId - $NormalReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(mpileup)
    ) from Mpileup4ControFREEC_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        sampleType,
        file(outFileName)
    ) into gatherMpileups_out_ch0

    script:
    outFileName = (sampleType == "T") ? TumorReplicateId + ".pileup.gz" : NormalReplicateId + ".pileup.gz"
    """
    scatters=`ls -1v *.pileup.gz`
    zcat \$scatters | \\
    bgzip --threads ${task.cpus} -c > ${outFileName}
    """
}

mpileupOutNormal = Channel.create()
mpileupOutTumor = Channel.create()

gatherMpileups_out_ch0
    .choice(
        mpileupOutTumor, mpileupOutNormal
    ) { it[2] == "T" ? 0 : 1 }

gatherMpileups_out_ch0 = mpileupOutTumor.combine(mpileupOutNormal, by: [0,1])

// run ControlFREEC : adopted from nfcore sarek
process 'ControlFREEC' {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/15_controlFREEC/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        _,
        file(mpileupTumor),
        _,
        file(mpileupNormal)
    ) from gatherMpileups_out_ch0

    set(
        file(RefChrDir),
        file(RefChrLen),
        file(DBSNP),
        file(DBSNPIdx)
    ) from Channel.value(
        [ reference.RefChrDir,
          reference.RefChrLen,
          database.DBSNP,
          database.DBSNPIdx ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}.pileup.gz_CNVs"),
        file("${TumorReplicateId}.pileup.gz_ratio.txt"),
        file("${TumorReplicateId}.pileup.gz_BAF.txt")
    ) into ControlFREEC_out_ch0

    script:
    config = "${TumorReplicateId}_vs_${NormalReplicateId}.config.txt"
    gender = genderMap[TumorReplicateId]

    read_orientation = (single_end) ? "0" : "FR"
    minimalSubclonePresence = (params.WES) ? 30 : 20
    degree = (params.WES) ? 1 : 4
    noisyData = (params.WES) ? "TRUE" : "FALSE"
    window = (params.WES) ? 0 : 50000
    breakPointType = (params.WES) ? 4 : 2
    breakPointThreshold = (params.WES) ? "1.2" : "0.8"
    printNA = (params.WES) ? "FALSE" :  "TRUE"
    readCountThreshold = (params.WES) ? 50 : 10
    minimalCoveragePerPosition = (params.WES) ? 5 : 0
    captureRegions = (params.WES) ? reference.RegionsBed : ""
    """
    rm -f ${config}
    touch ${config}
    echo "[general]" >> ${config}
    echo "BedGraphOutput = TRUE" >> ${config}
    echo "chrFiles = \${PWD}/${RefChrDir.fileName}" >> ${config}
    echo "chrLenFile = \${PWD}/${RefChrLen.fileName}" >> ${config}
    echo "coefficientOfVariation = 0.05" >> ${config}
    echo "contaminationAdjustment = TRUE" >> ${config}
    echo "forceGCcontentNormalization = 0" >> ${config}
    echo "maxThreads = ${task.cpus}" >> ${config}
    echo "minimalSubclonePresence = ${minimalSubclonePresence}" >> ${config}
    echo "ploidy = 2,3,4" >> ${config}
    echo "degree = ${degree}" >> ${config}
    echo "noisyData = ${noisyData}" >> ${config}
    echo "sex = ${gender}" >> ${config}
    echo "window = ${window}" >> ${config}
    echo "breakPointType = ${breakPointType}" >> ${config}
    echo "breakPointThreshold = ${breakPointThreshold}" >> ${config}
    echo "printNA = ${printNA}" >> ${config}
    echo "readCountThreshold = ${readCountThreshold}" >> ${config}
    echo "" >> ${config}
    echo "[control]" >> ${config}
    echo "inputFormat = pileup" >> ${config}
    echo "mateFile = \${PWD}/${mpileupNormal}" >> ${config}
    echo "mateOrientation = ${read_orientation}" >> ${config}
    echo "" >> ${config}
    echo "[sample]" >> ${config}
    echo "inputFormat = pileup" >> ${config}
    echo "mateFile = \${PWD}/${mpileupTumor}" >> ${config}
    echo "mateOrientation = ${read_orientation}" >> ${config}
    echo "" >> ${config}
    echo "[BAF]" >> ${config}
    echo "SNPfile = ${DBSNP.fileName}" >> ${config}
    echo "minimalCoveragePerPosition = ${minimalCoveragePerPosition}" >> ${config}
    echo "" >> ${config}
    echo "[target]" >> ${config}
    echo "captureRegions = ${captureRegions}" >> ${config}
    $FREEC -conf ${config}
    """
}


process 'ControlFREECviz' {

    tag "$TumorReplicateId"

    // makeGraph.R and assess_significance.R seem to be instable
    errorStrategy 'ignore'

    publishDir "$params.outputDir/$TumorReplicateId/15_controlFREEC/",
        mode: params.publishDirMode


    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(cnvTumor),
        file(ratioTumor),
        file(bafTumor),
    ) from ControlFREEC_out_ch0

    output:
    set(
        file("*.txt"),
        file("*.png"),
        file("*.bed"),
        file("*.circos")
     ) into ControlFREECviz_out_ch0


    script:
    """
    cat ${workflow.projectDir}/bin/assess_significance.R | R --slave --args ${cnvTumor} ${ratioTumor}
    cat ${workflow.projectDir}/bin/makeGraph.R | R --slave --args 2 ${ratioTumor} ${bafTumor}
    $PERL ${workflow.projectDir}/bin/freec2bed.pl -f ${ratioTumor} > ${TumorReplicateId}.bed
    $PERL ${workflow.projectDir}/bin/freec2circos.pl -f ${ratioTumor} > ${TumorReplicateId}.circos
    """
}


Channel
    .fromPath(reference.RefChrLen)
    .splitCsv(sep: "\t")
    .map { row -> row[1] }
    .set { chromosomes_ch }

process 'SequenzaUtils' {

    tag "$TumorReplicateId"

    conda "assets/py3_icbi.yml"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(tumorBAM),
        file(tumorBAI),
        file(normalBAM),
        file(normalBAI)
    ) from GatherRealignedBamFiles_out_Sequenza_ch0
    each chromosome from chromosomes_ch

    set(
        file(RefFasta),
        file(RefIdx),
        file(RefDict),
        file(SequnzaGC)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.SequenzaGC ]
    )

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${chromosome}_${TumorReplicateId}_seqz.gz")
    )  into SequenzaUtils_out_ch0

    script:
    """
    sequenza-utils \\
        bam2seqz \\
        --fasta ${RefFasta} \\
        --tumor ${tumorBAM} \\
        --normal ${normalBAM} \\
        -gc ${SequnzaGC} \\
        --chromosome ${chromosome} \\
        | \\
    sequenza-utils \\
        seqz_binning \\
        -w 50 \\
        -s - \\
        | \\
    $BGZIP \\
        --threads ${task.cpus} -c > ${chromosome}_${TumorReplicateId}_seqz.gz

    """
}

process gatherSequenzaInput {

    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(chromosome_seqz_binned)
    ) from SequenzaUtils_out_ch0
        .groupTuple(by: [0,1])


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_seqz.gz")
    ) into gatherSequenzaInput_out_ch0

    script:
    """
    first=1
    scatters=`ls -1v *_${TumorReplicateId}_seqz.gz`
    for f in \$scatters
    do
        if [ "\$first" ]
        then
            zcat "\$f"
            first=
        else
            zcat "\$f" | tail -n +2
        fi
    done | \\
    bgzip --threads ${task.cpus} -c > ${TumorReplicateId}_seqz.gz
    sync
    """
}

process Sequenza {

    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/16_Sequenza/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(seqz_file)
    ) from gatherSequenzaInput_out_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_segments.txt"),
        file("${TumorReplicateId}_confints_CP.txt")
    ) into Sequenza_out_Clonality_ch0
    file("${TumorReplicateId}_*.{png,pdf,txt}")

    script:
    gender = genderMap[TumorReplicateId]
    """
    $RSCRIPT \\
        ${workflow.projectDir}/bin/SequenzaScript.R \\
        ${seqz_file} \\
        ${TumorReplicateId} \\
        ${gender}
    """
}


clonality_input = Ascat_out_Clonality_ch0.combine(Sequenza_out_Clonality_ch0, by: [0,1])
    .map {

        it ->
        def TumorReplicateId = it[0]
        def NormalReplicateId = it[1]
        def ascat_CNVs = it[2]
        def ascat_purity  = it[3]
        def seqz_CNVs  = it[4]
        def seqz_purity  = it[5]

        def ascatOK = true

        def fileReader = ascat_purity.newReader()

        def line = fileReader.readLine()
        line = fileReader.readLine()
        fileReader.close()
        def (purity, ploidy) = line.split("\t")
        if(purity == "0" || ploidy == "0" ) {
            ascatOK = false
        }

        fileReader = ascat_CNVs.newReader()

        line = fileReader.readLine()
        fileReader.close()
        def fields = line.split("\t")
        if(fields.size() < 5) {
            ascatOK = false
        }
        return [ TumorReplicateId, NormalReplicateId, file(ascat_CNVs), file(ascat_purity), file(seqz_CNVs), file(seqz_purity), ascatOK ]
    }



process 'Clonality' {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/17_Clonality/",
        mode: params.publishDirMode


    conda "assets/py3_icbi.yml"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(hc_vep_vcf),
        file(hc_vep_idx),
        file(ascat_CNVs),
        file(ascat_purity),
        file(seqz_CNVs),
        file(seqz_purity),
        ascatOK
    ) from mkPhasedVCF_out_Clonality_ch0
        .combine(clonality_input, by: [0,1])


    output:
    set(
        TumorReplicateId,
        file("${TumorReplicateId}_CCFest.tsv"),
    ) into Clonality_out_ch0

    script:
    def seg_opt = ""
    if (ascatOK && ! params.use_sequenza_cnvs) {
        seg_opt = "--seg ${ascat_CNVs}"
        purity_opt = "--purity ${ascat_purity}"
    } else {
        if(! params.use_sequenza_cnvs) {
             log.warn "WARNING: changed from ASCAT to Sequenza purity and segments, ASCAT did not produce results"
        }
        seg_opt = "--seg_sequenza ${seqz_CNVs}"
        purity_opt = "--purity_sequenza ${seqz_purity}"
    }
    """
    mkCCF_input.py \\
        --PatientID ${TumorReplicateId} \\
        --vcf ${hc_vep_vcf} \\
        ${seg_opt} \\
        ${purity_opt} \\
        --min_vaf 0.01 \\
        --result_table ${TumorReplicateId}_segments.txt
    ${RSCRIPT} \\
        ${workflow.projectDir}/bin/CCF.R \\
        ${TumorReplicateId}_segments.txt \\
        ${TumorReplicateId}_CCFest.tsv \\
    """
}


// END CNVs


// HLA TYPING

process 'mhc_extract' {
    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file(tumor_BAM_aligned_sort_mkdp),
        file(tumor_BAI_aligned_sort_mkdp)
    ) from MarkDuplicatesTumor_out_ch0


    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${mhcReads_1}"),
        file("${mhcReads_2}")
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
    mhcReads_2 = (single_end) ? "NO_FILE" : TumorReplicateId + "_reads_mhc_R2.fastq.gz"

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
        touch NO_FILE

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
    ) from reads_tumor_hla_ch
    val yaraIdx from Channel.value(reference.YaraIndex)

    output:
    set (
        TumorReplicateId,
        file("mapped_{1,2}.bam")
    ) into fished_reads

    script:
    yara_cpus = ((task.cpus - 2).compareTo(2) == -1) ? 2 : (task.cpus - 2)
    samtools_cpus = ((task.cpus - yara_cpus).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus)
    if (single_end)
        """
        $YARA -e 3 -t $yara_cpus -f bam ${yaraIdx} ${readsFWD} | \\
            $SAMTOOLS view -@ $samtools_cpus -h -F 4 -b1 -o mapped_1.bam
        """
    else
        """
        $YARA -e 3 -t $yara_cpus -f bam ${yaraIdx} ${readsFWD} ${readsREV} | \\
            $SAMTOOLS view -@ $samtools_cpus -h -F 4 -b1 -o mapped_1.bam
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
    ) from fished_reads

    output:
    set (
        TumorReplicateId,
        file("${TumorReplicateId}_optitype_result.tsv")
    ) into optitype_output
    file("${TumorReplicateId}_optitype_coverage_plot.pdf")

    script:
    """
    $PYTHON $OPTITYPE -i ${reads} -e 1 -b 0.009 --dna -o ./tmp && \\
    mv ./tmp/*/*_result.tsv ./${TumorReplicateId}_optitype_result.tsv && \\
    mv ./tmp/*/*_coverage_plot.pdf ./${TumorReplicateId}_optitype_coverage_plot.pdf && \\
    rm -rf ./tmp/
    """
}

if (have_RNAseq) {
    process 'pre_map_hla_RNA' {
        tag "$TumorReplicateId"

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            readRNAFWD,
            readRNAREV
        ) from reads_tumor_optitype_ch
        val yaraIdx from Channel.value(reference.YaraIndex)

        output:
        set (
            TumorReplicateId,
            file("mapped_{1,2}.bam")
        ) into fished_reads_RNA

        script:
        yara_cpus = ((task.cpus - 2).compareTo(2) == -1) ? 2 : (task.cpus - 2)
        samtools_cpus = ((task.cpus - yara_cpus).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus)
        if (single_end_RNA)
            """
            $YARA -e 3 -t $yara_cpus -f bam ${yaraIdx} ${readRNAFWD} | \\
                $SAMTOOLS view -@ $samtools_cpus -h -F 4 -b1 -o mapped_1.bam
            """
        else
            """
            $YARA -e 3 -t $yara_cpus -f bam ${yaraIdx} ${readRNAFWD} ${readRNAREV} | \\
                $SAMTOOLS view -@ $samtools_cpus -h -F 4 -b1 -o mapped_1.bam
            """
    }

    process 'OptiType_RNA' {
        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/08_OptiType/",
            mode: params.publishDirMode

        input:
        set (
            TumorReplicateId,
            file(reads)
        ) from fished_reads_RNA

        output:
        set (
            TumorReplicateId,
            file("${TumorReplicateId}_optitype_RNA_result.tsv")
        ) into optitype_RNA_output
        file("${TumorReplicateId}_optitype_RNA_coverage_plot.pdf")

        script:
        """
        $PYTHON $OPTITYPE -i ${reads} -e 1 -b 0.009 --dna -o ./tmp && \\
        mv ./tmp/*/*_result.tsv ./${TumorReplicateId}_optitype_RNA_result.tsv && \\
        mv ./tmp/*/*_coverage_plot.pdf ./${TumorReplicateId}_optitype_RNA_coverage_plot.pdf && \\
        rm -rf ./tmp/
        """
    }

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
        file(readsREV)
    ) from reads_tumor_hlaHD_ch
    val frData from Channel.value(reference.HLAHDFreqData)
    file gSplit from Channel.value(reference.HLAHDGeneSplit)
    val dict from Channel.value(reference.HLAHDDict)

    output:
    set (
        TumorReplicateId,
        file("**/*_final.result.txt")
    ) into (
        hlahd_output,
        hlahd_mixMHC2_pred_ch0
    )

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
        opti_out_rna,
        hlahd_out,
        custom_hlas
    ) from optitype_output
        .combine(optitype_RNA_output, by: 0)
        .combine(hlahd_output, by: 0)
        .combine(custom_hlas_ch, by: 0)

    output:
    set (
        TumorReplicateId,
        file("${TumorReplicateId}_hlas.txt")
    ) into (hlas, hlas_neoFuse)

    script:
    def user_hlas = custom_hlas.name != 'NO_FILE' ? "--custom $custom_hlas" : ''
    def rna_hlas = have_RNAseq ? "--opti_out_RNA $opti_out_rna" : ''
    """
    HLA_parser.py \\
        --opti_out ${opti_out} \\
        --hlahd_out ${hlahd_out} \\
        ${rna_hlas} \\
        ${user_hlas} \\
        --ref_hlas ${baseDir}/assets/pVACseqAlleles.txt \\
        > ./${TumorReplicateId}_hlas.txt
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

if (have_RNAseq) {
    process Neofuse {

        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/10_NeoFuse/",
            mode: params.publishDirMode

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            readRNAFWD,
            readRNAREV,
            hla_types
        ) from reads_tumor_neofuse_ch
            .combine(hlas_neoFuse, by: 0)

        file STARidx from file(reference.STARidx)
        file RefFasta from file(reference.RefFasta)
        file AnnoFile from file(reference.AnnoFile)

        output:
        set (
            TumorReplicateId,
            file("**/${TumorReplicateId}.tpm.txt") // TODO: consider changing wildcards to expected dirnames and filenames
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
                -g ${RefFasta} \\
                -a ${AnnoFile} \\
                -N ${params.netMHCpan} \\
                -C ${hla_types} \\
                -k true
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
                -g ${RefFasta} \\
                -a ${AnnoFile} \\
                -N ${params.netMHCpan} \\
                -C ${hla_types} \\
                -k true
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
            file(vep_phased_vcf_gz),
            file(vep_phased_vcf_gz_tbi),
            file(vep_somatic_vcf_gz),
            file(vep_somatic_vcf_gz_tbi),
            file(final_tpm),
            file(RNA_bam),
            file(RNA_bai),
        ) from mkPhasedVCF_out_ch0
            .combine(final_file, by: 0)
            .combine(star_bam_file, by: 0)
        set(
        file(RefFasta),
        file(RefIdx),
        ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx ])

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("${TumorReplicateId}_vep_somatic_gx.vcf.gz"),
            file("${TumorReplicateId}_vep_somatic_gx.vcf.gz.tbi")
        ) into vcf_vep_ex_gz
        set(
            TumorReplicateId,
            NormalReplicateId,
            file("${TumorReplicateId}_${NormalReplicateId}_vep_phased.vcf.gz"),
            file("${TumorReplicateId}_${NormalReplicateId}_vep_phased.vcf.gz.tbi"),
            file("${TumorReplicateId}_vep_somatic_gx.vcf.gz"),
            file("${TumorReplicateId}_vep_somatic_gx.vcf.gz.tbi")
        ) into (
            gene_annotator_out_getseq_ch0,
            gene_annotator_out_mixMHC2pred_ch0
        )

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
} else { // no RNAseq data

    process no_gene_annotator {

        tag "$TumorReplicateId"

        input:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(vep_phased_vcf_gz),
            file(vep_phased_vcf_gz_tbi),
            file(vep_somatic_vcf_gz),
            file(vep_somatic_vcf_gz_tbi)
        ) from mkPhasedVCF_out_ch0

        output:
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(vep_somatic_vcf_gz),
            file(vep_somatic_vcf_gz_tbi)
        ) into vcf_vep_ex_gz
        set(
            TumorReplicateId,
            NormalReplicateId,
            file(vep_phased_vcf_gz),
            file(vep_phased_vcf_gz_tbi),
            file(vep_somatic_vcf_gz),
            file(vep_somatic_vcf_gz_tbi)
        ) into (
            gene_annotator_out_getseq_ch0,
            gene_annotator_out_mixMHC2pred_ch0
        )

        script:
        // do nothing
        """
        """
    }
}

/*
Run pVACseq
*/

process 'pVACseq' {
    tag "$TumorReplicateId"

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
    ) from mkPhasedVCF_out_pVACseq_ch0
        .combine(vcf_vep_ex_gz, by: 0)
        .combine(hlas.splitText(), by: 0)

    output:
    set(
        TumorReplicateId,
        file("**/MHC_Class_I/*filtered.tsv"),
        file("**/MHC_Class_I/*all_epitopes.tsv")
    ) optional true into mhcI_out_f

    set(
        TumorReplicateId,
        file("**/MHC_Class_II/*filtered.tsv"),
        file("**/MHC_Class_II/*all_epitopes.tsv")
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

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/MHC_Class_I/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        '*filtered.tsv',
        '*all_epitopes.tsv'
    ) from mhcI_out_f.groupTuple(by: 0)

    output:
    set(
        TumorReplicateId,
        file("${TumorReplicateId}_MHCI_filtered.tsv")
    ) into (
        MHCI_final_ranked,
        MHCI_final_immunogenicity,
        concat_mhcI_filtered_files_out_addCCF_ch0
    )
    set(
        TumorReplicateId,
        file("${TumorReplicateId}_MHCI_all_epitopes.tsv")
    ) into (
        MHCI_all_epitopes,
        concat_mhcI_all_files_out_addCCF_ch0
    )

    script:
    """
    sed -e '2,\${/^Chromosome/d' -e '}' *filtered.tsv > ${TumorReplicateId}_MHCI_filtered.tsv
    sed -e '2,\${/^Chromosome/d' -e '}' *all_epitopes.tsv > ${TumorReplicateId}_MHCI_all_epitopes.tsv
    """
}

process concat_mhcII_files {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/MHC_Class_II/",
        mode: params.publishDirMode

    input:
    set(
        TumorReplicateId,
        '*filtered.tsv',
        '*all_epitopes.tsv'
    ) from mhcII_out_f.groupTuple(by: 0)

    output:
    set(
        TumorReplicateId,
        file("${TumorReplicateId}_MHCII_filtered.tsv")
    ) into (
        MHCII_final_ranked,
        concat_mhcII_filtered_files_out_addCCF_ch0
    )
    set(
        TumorReplicateId,
        file("${TumorReplicateId}_MHCII_all_epitopes.tsv")
    ) into (
        MHCII_all_epitopes,
        concat_mhcII_all_files_out_addCCF_ch0
    )

    script:
    """
    sed -e '2,\${/^Chromosome/d' -e '}' *filtered.tsv > ${TumorReplicateId}_MHCII_filtered.tsv
    sed -e '2,\${/^Chromosome/d' -e '}' *all_epitopes.tsv > ${TumorReplicateId}_MHCII_all_epitopes.tsv
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
    mkdir ./MHC_Class_I/
    pvacseq generate_condensed_ranked_report \\
        -m lowest \\
        $pvacseq_mhcI_file \\
        ./MHC_Class_I/${TumorReplicateId}_MHCI_filtered.condensed.ranked.tsv
    mkdir ./MHC_Class_II/
    pvacseq generate_condensed_ranked_report \\
        -m lowest $pvacseq_mhcII_file \\
        ./MHC_Class_II/${TumorReplicateId}_MHCII_filtered.condensed.ranked.tsv
    """
}

process 'pVACtools_generate_protein_seq' {
    tag "$TumorReplicateId"

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        vep_phased_vcf_gz,
        vep_phased_vcf_idx,
        vep_tumor_vcf_gz,
        vep_tumor_vcf_idx
    ) from gene_annotator_out_getseq_ch0

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_protSeq.fasta")
    ) optional true into pVACtools_generate_protein_seq

    script:
    """
    pvacseq generate_protein_fasta \\
        -p ${vep_phased_vcf_gz} \\
        -s ${TumorReplicateId} \\
        ${vep_tumor_vcf_gz} \\
        31 \\
        ${TumorReplicateId}_protSeq.fasta
    """
}

process 'pepare_mixMHC2_seq' {
    tag "$TumorReplicateId"

    // publishDir "$params.outputDir/$TumorReplicateId/12_mixMHC2pred/",
    //     mode: params.publishDirMode

    conda 'assets/py3_icbi.yml'

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        protSeq_fasta,
        hlahd_allel_file
    ) from pVACtools_generate_protein_seq
        .combine(hlahd_mixMHC2_pred_ch0, by:0)

    output:
    set(
        TumorReplicateId,
        NormalReplicateId,
        file("${TumorReplicateId}_peptides.fasta")
    ) optional true into pepare_mixMHC2_seq_out_ch0
    file("${TumorReplicateId}_alleles_translated.txt") optional true into pepare_mixMHC2_seq_out_ch1

    script:
    """
    pepChopper.py \\
        --pep_len ${params.mhcii_epitope_len.split(",").join(" ")} \\
        --fasta_in ${protSeq_fasta} \\
        --fasta_out ${TumorReplicateId}_peptides.fasta
    translate_HLAHD2mixMHC2pred.py \\
        --HLAHD_file ${hlahd_allel_file} \\
        --Allele_file ${baseDir}/assets/Alleles_list.txt \\
        --translated_file ${TumorReplicateId}_alleles_translated.txt
    """
}


process mixMHC2pred {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/12_mixMHC2pred",
        mode: params.publishDirMode

    conda 'assets/py3_icbi.yml'

    input:
    set(
        TumorReplicateId,
        NormalReplicateId,
        mut_peps,
        _,
        _,
        vep_somatic_gx_vcf_gz,
        vep_somatic_gx_vcf_gz_tbi
    ) from pepare_mixMHC2_seq_out_ch0
        .combine(gene_annotator_out_mixMHC2pred_ch0, by: [0, 1])
    val allelesFile from pepare_mixMHC2_seq_out_ch1

    output:
    file("${TumorReplicateId}_mixMHC2pred_all.tsv")
    file("${TumorReplicateId}_mixMHC2pred_filtered.tsv")

    script:
    alleles = file(allelesFile).readLines().join(" ")
    """
    ${MiXMHC2PRED} \\
        -i ${mut_peps} \\
        -o ${TumorReplicateId}_mixMHC2pred.tsv \\
        -a ${alleles}
    parse_mixMHC2pred.py \\
        --vep_vcf ${vep_somatic_gx_vcf_gz} \\
        --pep_fasta ${mut_peps} \\
        --mixMHC2pred_result ${TumorReplicateId}_mixMHC2pred.tsv \\
        --out ${TumorReplicateId}_mixMHC2pred_all.tsv \\
        --sample_name ${TumorReplicateId} \\
        --normal_name ${NormalReplicateId}
    awk \\
        '{
            if (\$0 ~ /\\#/) { print }
            else { if (\$18 <= 2) { print } }
        }' ${TumorReplicateId}_mixMHC2pred_all.tsv > ${TumorReplicateId}_mixMHC2pred_filtered.tsv
    """
}

// add CCF clonality to neoepitopes result files
process addCCF {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/17_Clonality",
        mode: params.publishDirMode

    conda 'assets/py3_icbi.yml'

    input:
    set(
        TumorReplicateId,
        file(epitopes),
        file(CCF)
    ) from concat_mhcI_filtered_files_out_addCCF_ch0
        .concat(
                concat_mhcI_all_files_out_addCCF_ch0,
                concat_mhcII_filtered_files_out_addCCF_ch0,
                concat_mhcII_all_files_out_addCCF_ch0
        )
        .combine(Clonality_out_ch0, by: 0)

    output:
    file(outfile)

    script:
    outfile = epitopes.baseName + "_ccf.tsv"
    """
    add_CCF.py \\
        --neoepitopes ${epitopes} \\
        --ccf ${CCF} \\
        --outfile ${outfile}
    """
}


/*
  Immunogenicity scoring
*/

process csin {
    tag "$TumorReplicateId"

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/",
        mode: params.publishDirMode

    input:
    set (
        TumorReplicateId,
        "*_MHCI_all_epitopes.tsv",
        "*_MHCII_all_epitopes.tsv"
    ) from MHCI_all_epitopes
    .combine(MHCII_all_epitopes, by:0)

    output:
    file("${TumorReplicateId}_CSiN.tsv")

    script:
    """
    CSiN.py --MHCI_tsv *_MHCI_all_epitopes.tsv \\
        --MHCII_tsv *_MHCII_all_epitopes.tsv \\
        --rank $params.csin_rank \\
        --ic50 $params.csin_ic50 \\
        --gene_exp $params.csin_gene_exp \\
        --output ./${TumorReplicateId}_CSiN.tsv
    """
}


process immunogenicity_scoring {
    tag "$TumorReplicateId"

    // TODO: check why sometimes this fails: workaround ignore errors
    errorStrategy 'ignore'

    publishDir "$params.outputDir/$TumorReplicateId/11_pVACseq/MHC_Class_I/",
        mode: params.publishDirMode

    input:
    set (
        TumorReplicateId,
        pvacseq_file
    ) from MHCI_final_immunogenicity
    // val(TumorReplicateId) from mhCI_tag_immunogenicity
    // file pvacseq_file from MHCI_final_immunogenicity

    output:
    file("${TumorReplicateId}_immunogenicity.tsv") optional true

    script:
    """
    get_epitopes.py \\
        --pvacseq_out $pvacseq_file \\
        --sample_id $TumorReplicateId \\
        --output ./${TumorReplicateId}_epitopes.tsv
    NR_EPI=`wc -l ./${TumorReplicateId}_epitopes.tsv | cut -d" " -f 1`
    if [ \$NR_EPI -gt 1 ]; then
        NeoAg_immunogenicity_predicition_GBM.R \\
            ./${TumorReplicateId}_epitopes.tsv ./${TumorReplicateId}_temp_immunogenicity.tsv
        immuno_score.py \\
            --pvacseq_tsv $pvacseq_file \\
            --score_tsv ${TumorReplicateId}_temp_immunogenicity.tsv \\
            --output ${TumorReplicateId}_immunogenicity.tsv
    fi
    """
}

if(params.TCR) {
    process mixcr_DNA_tumor {
        tag "$TumorReplicateId"

        publishDir "$params.outputDir/$TumorReplicateId/13_MiXCR",
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
            file("${TumorReplicateId}_mixcr_DNA.clonotypes.ALL.txt"),
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

        publishDir "$params.outputDir/$TumorReplicateId/13_MiXCR",
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
            file("${NormalReplicateId}_mixcr_DNA.clonotypes.ALL.txt"),
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

    if (have_RNAseq) {
        process mixcr_RNA {
            tag "$TumorReplicateId"

            publishDir "$params.outputDir/$TumorReplicateId/13_MiXCR",
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
                file("${TumorReplicateId}_mixcr_RNA.clonotypes.ALL.txt"),
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
        file("*")
    )   from ch_fastqc
            .combine(ch_fastp_tumor, by: [0,1])
            .combine(ch_fastp_normal, by: [0,1])
            .combine(ch_fastqc_trimmed, by: [0,1])
            .combine(ch_fastp_RNAseq, by: [0,1])
            .combine(ch_fastqc_trimmed_RNAseq, by: [0,1])
            .combine(alignmentMetricsTumor_ch, by: [0,1])
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
        if (params.references.size() != 18) exit 1, """
        ERROR: Not all References needed found in configuration
        Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
        """
        return [
            'RefFasta'          : checkParamReturnFileReferences("RefFasta"),
            'RefIdx'            : checkParamReturnFileReferences("RefIdx"),
            'RefDict'           : checkParamReturnFileReferences("RefDict"),
            'RefChrLen'         : checkParamReturnFileReferences("RefChrLen"),
            'RefChrDir'         : checkParamReturnFileReferences("RefChrDir"),
            'BwaRef'            : checkParamReturnFileReferences("BwaRef"),
            'VepFasta'          : checkParamReturnFileReferences("VepFasta"),
            'BaitsBed'          : checkParamReturnFileReferences("BaitsBed"),
            'RegionsBed'        : checkParamReturnFileReferences("RegionsBed"),
            'YaraIndex'         : checkParamReturnFileReferences("YaraIndex"),
            'HLAHDFreqData'     : checkParamReturnFileReferences("HLAHDFreqData"),
            'HLAHDGeneSplit'    : checkParamReturnFileReferences("HLAHDGeneSplit"),
            'HLAHDDict'         : checkParamReturnFileReferences("HLAHDDict"),
            'STARidx'           : checkParamReturnFileReferences("STARidx"),
            'AnnoFile'          : checkParamReturnFileReferences("AnnoFile"),
            'acLoci'            : checkParamReturnFileReferences("acLoci"),
            'acLociGC'          : checkParamReturnFileReferences("acLociGC"),
            'SequenzaGC'        : checkParamReturnFileReferences("SequenzaGC")
        ]
    } else {
        if (params.references.size() != 17) exit 1, """
        ERROR: Not all References needed found in configuration
        Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
        """
        return [
            'RefFasta'          : checkParamReturnFileReferences("RefFasta"),
            'RefIdx'            : checkParamReturnFileReferences("RefIdx"),
            'RefDict'           : checkParamReturnFileReferences("RefDict"),
            'RefChrLen'         : checkParamReturnFileReferences("RefChrLen"),
            'RefChrDir'         : checkParamReturnFileReferences("RefChrDir"),
            'BwaRef'            : checkParamReturnFileReferences("BwaRef"),
            'VepFasta'          : checkParamReturnFileReferences("VepFasta"),
            'YaraIndex'         : checkParamReturnFileReferences("YaraIndex"),
            'HLAHDFreqData'     : checkParamReturnFileReferences("HLAHDFreqData"),
            'HLAHDGeneSplit'    : checkParamReturnFileReferences("HLAHDGeneSplit"),
            'HLAHDDict'         : checkParamReturnFileReferences("HLAHDDict"),
            'STARidx'           : checkParamReturnFileReferences("STARidx"),
            'AnnoFile'          : checkParamReturnFileReferences("AnnoFile"),
            'acLoci'            : checkParamReturnFileReferences("acLoci"),
            'acLociGC'          : checkParamReturnFileReferences("acLociGC"),
            'SequenzaGC'        : checkParamReturnFileReferences("SequenzaGC")
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
    log.info "--        U S A G E       NOT UP TO DATE PLEASE FIX --"
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
    log.info "tumorSampleName,readsTumorFWD,readsTumorREV,normalSampleName,readsNormalFWD,readsNormalREV,HLAfile,group"
    log.info "tumor1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,normal1,Normal1_reads_1.fastq,Normal1_reads_2.fastq,None,group1"
    log.info "tumor2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,normal2,Normal2_reads_1.fastq,Normal2_reads_2.fastq,None,group1"
    log.info "..."
    log.info "sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,NormalN_reads_1.fastq,NormalN_reads_2.fastq,None,groupX"
    log.info ""
    log.info "CSV-file, single-end T only reads:"
    log.info "tumorSampleName,readsTumorFWD,readsTumorREV,readsNormalFWD,readsNormalREV,HLAfile,group"
    log.info "tumor1,Tumor1_reads_1.fastq,Tumor1_reads_2.fastq,noname,NO_FILE,NO_FILE,None,group1"
    log.info "tumor2,Tumor2_reads_1.fastq,Tumor2_reads_2.fastq,noname,NO_FILE,NO_FILE,None,group1"
    log.info "..."
    log.info "sampleN,TumorN_reads_1.fastq,TumorN_reads_2.fastq,None,,None,groupX"
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
