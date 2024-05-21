#!/usr/bin/env nextflow

// enable DSL 1
nextflow.enable.dsl = 1

import org.yaml.snakeyaml.Yaml
import java.nio.file.Files

log.info ""
log.info " NEXTFLOW ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
log.info "-------------------------------------------------------------------------"
log.info "          Nextfow NeoAntigen Prediction Pipeline - nextNEOpi    "
log.info "-------------------------------------------------------------------------"
log.info ""
log.info " Features: "
log.info " - somatic variants from tumor + matched normal samples"
log.info " - CNV analysis"
log.info " - tumor muational burden"
log.info " - class I and class II HLA typing"
log.info " - gene fusion peptide prediction using RNAseq data"
log.info " - peptide MHC binding perdiction"
log.info " - clonality of neoantigens"
log.info " - expression of neoantigens"
log.info ""
log.info "-------------------------------------------------------------------------"
log.info "C O N F I G U R A T I O N"
log.info ""
log.info "Command Line: \t\t " + workflow.commandLine
log.info "Working Directory: \t " + workflow.workDir
log.info "Output Directory: \t " + params.outputDir
log.info ""
log.info "I N P U T"
log.info ""
log.info "batch file: \t\t " + params.batchFile
log.info ""
log.info "Please check --help for further instruction"
log.info "-------------------------------------------------------------------------"

// Check if License(s) were accepted
params.accept_license = false

if (params.accept_license) {
    acceptLicense()
} else {
    checkLicense()
}

/*
________________________________________________________________________________

                            C O N F I G U R A T I O N
________________________________________________________________________________
*/
if (params.help) exit 0, helpMessage()

// switch for enable/disable processes (debug/devel only: use if(params.RUNTHIS) { .... })
params.RUNTHIS = false

// default is not to get bams as input data
bamInput = false

// initialize RNA tag seq
have_RNA_tag_seq = params.RNA_tag_seq

// set and initialize the Exome capture kit
setExomeCaptureKit(params.exomeCaptureKit)

// check conda channels
if (params.enable_conda) {
    checkCondaChannels()
}

// check IEDB dir
check_iedb_dir(params.databases.IEDB_dir)

// check MHCflurry dir
check_mhcflurry_dir(params.databases.MHCFLURRY_dir)

// set and check references and databases
reference = defineResources('references', params.WES, params.HLAHD_DIR)
database = defineResources('databases', params.WES, params.HLAHD_DIR)

// create tmp dir and make sure we have the realpath for it
tmpDir = mkTmpDir(params.tmpDir)

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
summary['Pipeline Name']                 = 'icbi/nextNEOpi'
summary['Pipeline Version']              = workflow.manifest.version
summary['Batch file']                    = params.batchFile
summary['Read length']                   = params.readLength
summary['Exome capture kit']             = params.exomeCaptureKit
summary['Fasta Ref']                     = params.references.RefFasta
summary['MillsGold']                     = params.databases.MillsGold
summary['hcSNPS1000G']                   = params.databases.hcSNPS1000G
summary['HapMap']                        = params.databases.HapMap
summary['Cosmic']                        = params.databases.Cosmic
summary['DBSNP']                         = params.databases.DBSNP
summary['GnomAD']                        = params.databases.GnomAD
summary['GnomADfull']                    = params.databases.GnomADfull
summary['KnownIndels']                   = params.databases.KnownIndels
summary['BlastDB']                       = params.references.ProteinBlastDBdir
summary['priority variant Caller']       = params.primaryCaller
summary['Mutect 1 and 2 minAD']          = params.minAD
summary['VarScan min_cov']               = params.min_cov
summary['VarScan min_cov_tumor']         = params.min_cov_tumor
summary['VarScan min_cov_normal']        = params.min_cov_normal
summary['VarScan min_freq_for_hom']      = params.min_freq_for_hom
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
summary['Output dir']                    = params.outputDir
summary['Working dir']                   = workflow.workDir
summary['TMP dir']                       = tmpDir
summary['Current home']                  = "$HOME"
summary['Current user']                  = "$USER"
summary['Current path']                  = "$PWD"
summary['Picard maxRecordsInRam']        = params.maxRecordsInRam
summary['Script dir']                    = workflow.projectDir
summary['Config Profile']                = workflow.profile


if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(30)}: $v" }.join("\n")
log.info "-------------------------------------------------------------------------"

// End Summary

// determine the publishDirMode
def publishDirMode = get_publishMode(params.outputDir, params.publishDirMode)

// Check if we got a batch file
params.batchFile = null
if (params.batchFile == null) {
    log.error "No sample sheet specified, please use --batchFile to pass the sample sheet"
    exit(1)
}

def batchCSV = file(params.batchFile).splitCsv(header:true)

def validFQfields = ["sampleName",
                     "reads1",
                     "reads2",
                     "sampleType",
                    "HLAfile",
                    "sex"]

def validBAMfields = ["sampleName",
                      "bam",
                      "sampleType",
                      "HLAfile",
                      "sex"]

def validSampleTypes = ["tumor_DNA", "normal_DNA", "tumor_RNA"]

if (batchCSV.size() > 0) {

    if (batchCSV[0].keySet().sort() == validFQfields.sort()) {
        bamInput = false
    } else if (batchCSV[0].keySet().sort() == validBAMfields.sort()) {
        bamInput = true
    } else {
        exit 1, "Error: Incorrect fields in batch file, please check your batchFile"
    }
} else {
    exit 1, "Error: No samples found, please check your batchFile"
}

def raw_data = []
def custom_HLA_data = []

def t_map = [:]
def n_map = [:]
def r_map = [:]
def lt_map = [:]

for ( row in batchCSV ) {
    def meta  = [:]
    def reads = []

    meta.sampleName = row.sampleName.toString()
    meta.sampleType = row.sampleType.toString()

    t_map[meta.sampleName] = (t_map[meta.sampleName] == null) ? 0 : t_map[meta.sampleName]
    n_map[meta.sampleName] = (n_map[meta.sampleName] == null) ? 0 : n_map[meta.sampleName]

    if(row.sex) {
        if (row.sex.toLowerCase() in ["xx", "female"]) {
           meta.sex = "XX"
        } else if (row.sex.toLowerCase() in ["xy", "male"]) {
            meta.sex = "XY"
        } else if (row.sex.toLowerCase() in ["none", "na"]) {
            meta.sex = "None"
            println("WARNING: " + row.sampleName + " sex not specified will infer from data")
        } else {
            exit 1, "sex should be one of: XX, xx, XY, xy, Female, female, Male, male, None, none, NA, got: " + row.sex
        }
        meta.maleRef = (meta.sex ==  "XY") ? true : false
    } else {
        println("WARNING: " + row.sampleName + " sex not specified will infer from data")
        meta.sex = "None"
        meta.maleRef = true
    }

    if (row.sampleType in validSampleTypes) {
        t_map[meta.sampleName] += (row.sampleType == "tumor_DNA") ? 1 : 0
        n_map[meta.sampleName] += (row.sampleType == "normal_DNA") ? 1 : 0
    } else {
        exit 1, "Error: Incorrect sampleType [got: " + row.sampleType + " - allowed: " + validSampleTypes + "]: please check your batchFile"
    }

    // remove stuff that is not required in the HLAfile channel
    def meta_min = meta.clone()
    meta_min.keySet().removeAll(['sampleType', 'libType'])
    if (row.HLAfile) {
        custom_HLA_data.add([meta_min, file(row.HLAfile, checkIfExists: true)])
    } else {
        custom_HLA_data.add([meta_min, []])
    }

    if (row.sampleType == "tumor_RNA") {
        meta.have_RNA = true
        r_map[row.sampleName] = true
    }


    if (! row.bam) {
        meta.libType = "SE"
        if (row.reads1) { reads.add(file(row.reads1, checkIfExists: true)) }
        if (row.reads2) {
            reads.add(file(row.reads2, checkIfExists: true))
            meta.libType = "PE"
        }
        if (meta.sampleType != "tumor_RNA") {
            if (lt_map[meta.sampleName] == null) {
                lt_map[meta.sampleName] = meta.libType
            } else {
                if (lt_map[meta.sampleName] != meta.libType) {
                    exit 1, "Please do not mix pe and se for tumor/normal pairs: " + meta.sampleName + " - Not supported"
                }
            }
        }
        raw_data.add([meta, reads])
    } else {
        raw_data.add([meta, file(row.bam, checkIfExists: true)])
    }

}

// check if we have T/N DNA pairs for all patients
raw_data.each {
    record ->
    if (t_map[record[0].sampleName] == 0) {
        exit 1, "NO tumor DNA sample specified for: " + record[0].sampleName
    }
    if (n_map[record[0].sampleName] == 0) {
        exit 1, "NO normal DNA sample specified for: " + record[0].sampleName
    }
}

// update meta for RNA
raw_data.each {
    record ->
    record[0].have_RNA = (r_map[record[0].sampleName]) ? true : false
}
custom_HLA_data = custom_HLA_data.each {
    record ->
    record[0].have_RNA = (r_map[record[0].sampleName]) ? true : false
}.unique()

use_custom_hlas = (custom_HLA_data.size() > 0) ? true : false

batch_raw_data_ch = Channel.fromList(raw_data)
batch_custom_HLA_data_ch = Channel.fromList(custom_HLA_data)

if (bamInput == false) {
    batch_raw_data_ch.map {
            meta, fastq ->
            [ meta, fastq ]
        }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.unique().size() == 1
                    return [ meta, fastq.flatten().unique() ]
                multiple: fastq.unique().size() > 1
                    return [ meta, fastq.flatten().unique() ]
        }
        .set { fastq_ch }
} else {
    fastq_ch = Channel.of().branch{single: []; multiple: []}
}

// optional panel of normals file
pon_file = file(params.mutect2ponFile)

scatter_count = Channel.from(params.scatter_count)
padding = params.readLength + 100

MiXMHC2PRED   = ( params.MiXMHC2PRED != "" ) ? file(params.MiXMHC2PRED) : ""

// check HLAHD & OptiType
have_HLAHD = false
run_OptiType = (params.disable_OptiType) ? false : true

if (params.HLAHD_DIR != "") {
    HLAHD = file(params.HLAHD_DIR + "/bin/hlahd.sh")
    if (checkToolAvailable(HLAHD, "exists", "warn")) {
        HLAHD_DIR  = file(params.HLAHD_DIR)
        HLAHD_PATH = HLAHD_DIR + "/bin"
        if(params.HLAHD_module != "") {
            if (checkToolAvailable("bowtie2", "inPath", "warn", module=params.HLAHD_module)) {
                have_HLAHD = true
            }
        } else {
            if (checkToolAvailable("bowtie2", "inPath", "warn")) {
                have_HLAHD = true
            }
        }
    }
}
if (! have_HLAHD && run_OptiType) {
    log.warn "WARNING: HLAHD not available - can not predict Class II neoepitopes"
} else if (! have_HLAHD && ! run_OptiType && use_custom_hlas) {
    log.warn "WARNING: HLAHD not available and OptiType disabled - using only user supplied HLA types"
} else if (! have_HLAHD && ! run_OptiType && ! use_custom_hlas) {
    exit 1, "ERROR: HLAHD not available and OptiType disabled - can not predict HLA types"
}

// check if all tools are installed when not running conda or singularity
have_vep = false
if (! workflow.profile.contains('conda') && ! workflow.profile.contains('singularity')) {
    def execTools = ["fastqc", "fastp", "bwa", "samtools", "sambamba", "gatk", "vep", "bam-readcount",
                     "perl", "bgzip", "tabix", "bcftools", "yara_mapper", "python", "cnvkit.py",
                     "OptiTypePipeline.py", "alleleCounter", "freec", "Rscript", "java", "multiqc",
                     "sequenza-utils"]

    for (tool in execTools) {
        checkToolAvailable(tool, "inPath", "error")
    }

    VARSCAN = "java -jar " + file(params.VARSCAN)
    have_vep = true
} else {
    VARSCAN = "varscan "
}

// check if we have mutect1 installed
have_Mutect1 = false
if (params.MUTECT1 != "" && file(params.MUTECT1) && params.JAVA7 != "" && file(params.JAVA7)) {
    if(checkToolAvailable(params.JAVA7, "exists", "warn") && checkToolAvailable(params.MUTECT1, "exists", "warn")) {
        JAVA7 = file(params.JAVA7)
        MUTECT1 = file(params.MUTECT1)
        have_Mutect1 = true
    }
}

// check if we have GATK3 installed
have_GATK3 = false
if (params.GATK3 != "" && file(params.GATK3) && params.JAVA8 != "" && file(params.JAVA8) && ! workflow.profile.contains('conda') && ! workflow.profile.contains('singularity')) {
    if(checkToolAvailable(params.JAVA8, "inPath", "warn") && checkToolAvailable(params.GATK3, "exists", "warn")) {
        JAVA8 = file(params.JAVA8)
        GATK3 = file(params.GATK3)
        have_GATK3 = true
    }
} else if (workflow.profile.contains('singularity')) {
    JAVA8 = "java"
    GATK3 = "/usr/local/opt/gatk-3.8/GenomeAnalysisTK.jar"
    have_GATK3 = true
} else if (workflow.profile.contains('conda')) {
    JAVA8 = "java"
    GATK3 = "\$CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar"
    have_GATK3 = true
}


// check MIXCR licence
if (params.TCR && params.MIXCR_lic != "") {
    checkToolAvailable(params.MIXCR_lic, "exists", "error")
} else if (params.TCR && params.MIXCR_lic == "") {
    exit 1, "ERROR: no MiXCR license file specified, please provide a MiXCR license file in params.config or by using the --MIXCR_lic option.\nIf you do not have a MiXCR license you may:\n\ta) run nextNEOpi with --TCR false\n\tb) request one at https://licensing.milaboratories.com"
}

/*
________________________________________________________________________________

                                P R O C E S S E S
________________________________________________________________________________
*/

gatk4_chck_file = file(baseDir + "/assets/.gatk4_install_ok.chck")
gatk4_chck_file.append()
nextNEOpiENV_setup_ch0 = Channel.value("OK")

/*
*********************************************
**       P R E P R O C E S S I N G         **
*********************************************
*/

// Handle BAM input files. We need to convert BAMs to Fastq
if(bamInput) {
    process check_PE {
        label 'nextNEOpiENV'

        tag "$meta.sampleName : $meta.sampleType"

        input:
        tuple(
            val(meta),
            path(bam)
        ) from batch_raw_data_ch

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path(bam),
            stdout
        ) into bam_ch


        script:
        """
        check_pe.py $bam
        """
    }
    (bam_ch, check_seqLib_ch) = bam_ch.into(2)
    seqLibTypes_ok = Channel.value(check_seqLibTypes_ok(check_seqLib_ch))


    process bam2fastq {
        label 'nextNEOpiENV'

        tag "$meta.sampleName : $meta.sampleType"

        publishDir "${params.outputDir}/analyses/${meta.sampleName}/01_preprocessing",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(bam),
            val(libType),
            val(libOK)
        ) from bam_ch
            .combine(seqLibTypes_ok)

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path("${prefix}_R*.fastq.gz"),
        ) into (
            bam_fastq_ch
        )

        script:
        def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 4) / task.cpus)), 1)
        prefix = meta.sampleName + "_" + meta.sampleType
        meta.libType = libType
        if (libType == "PE")
            """
            samtools sort -@ ${task.cpus} -m ${STperThreadMem}G -l 0 -n ${bam} | \\
            samtools fastq \\
                -@ ${task.cpus} \\
                -c 5 \\
                -1 ${prefix}_R1.fastq.gz \\
                -2 ${prefix}_R2.fastq.gz \\
                -0 /dev/null -s /dev/null \\
                -n \\
                /dev/stdin
            """
        else if (libType == "SE")
            """
            samtools fastq \\
                -@ ${task.cpus} \\
                -n \\
                ${bam} | \\
                bgzip -@ ${task.cpus} -c /dev/stdin > ${prefix}_R1.fastq.gz
            """
    }
} else {
    bam_fastq_ch = channel.empty()
}
// END BAM input handling

// merge multi run or multi lane sample fastq
process merge_fastq {

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "${params.outputDir}/${meta.sampleName}/01_preprocessing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple meta, path("R1/R1.*"), path("R2/R2.*") from fastq_ch.multiple.map{
        meta, r -> {
            def r1 = []
            def r2 = []
            r.eachWithIndex{ v, ix -> ( ix & 1 ? r2 : r1 ) << v }
            tuple( meta, r1.sort(), r2.sort())
        }
    }

    output:
    tuple val(meta), path("*_R{1,2}.merged.fastq.gz") into merged_fastq_ch

    script:
    def prefix = meta.sampleName + "_" + meta.sampleType
    """
    cat R1/* > ${prefix}_R1.merged.fastq.gz
    cat R2/* > ${prefix}_R2.merged.fastq.gz
    """
}


// here we have our final raw fastq files:
// merged if multi lane runs are used
// bam2fastq if bam files are used
(raw_reads_ch, fastqc_reads_ch) = merged_fastq_ch.mix(fastq_ch.single, bam_fastq_ch).into(2)



// Common region files preparation for faster processing
if (params.WES) {
    process 'RegionsBedToIntervalList' {

        label 'nextNEOpiENV'

        tag 'RegionsBedToIntervalList'

        publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
            mode: publishDirMode

        input:
        tuple(
            path(RefDict),
            path(RegionsBed)
        ) from Channel.value(
            [ reference.RefDict,
            reference.RegionsBed ]
        )

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        path(
            "${RegionsBed.baseName}.interval_list"
        ) into (
            RegionsBedToIntervalList_out_ch0,
            RegionsBedToIntervalList_out_ch1
        )

        script:
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
        """
        gatk --java-options ${JAVA_Xmx} BedToIntervalList \\
            -I ${RegionsBed} \\
            -O ${RegionsBed.baseName}.interval_list \\
            -SD $RefDict
        """
    }

    process 'BaitsBedToIntervalList' {

        label 'nextNEOpiENV'

        tag 'BaitsBedToIntervalList'

        publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
            mode: publishDirMode

        input:
        tuple(
            path(RefDict),
            path(BaitsBed)
        ) from Channel.value(
            [ reference.RefDict,
            reference.BaitsBed ]
        )

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        path(
            "${BaitsBed.baseName}.interval_list"
        ) into BaitsBedToIntervalList_out_ch0

        script:
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
        """
        gatk --java-options ${JAVA_Xmx} BedToIntervalList \\
            -I ${BaitsBed} \\
            -O ${BaitsBed.baseName}.interval_list \\
            -SD $RefDict
        """
    }
} else {
    RegionsBedToIntervalList_out_ch0 = Channel.fromPath('NO_FILE')
    RegionsBedToIntervalList_out_ch1 = Channel.empty()
    BaitsBedToIntervalList_out_ch0 = Channel.empty()
}

process 'preprocessIntervalList' {

    label 'nextNEOpiENV'

    tag 'preprocessIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
        mode: publishDirMode

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )
    path(interval_list) from RegionsBedToIntervalList_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    path(
        outFileName
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
    outFileName = (params.WES) ? interval_list.baseName + "_merged_padded.interval_list" : "wgs_ScatterIntervalsByNs.interval_list"
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    if(params.WES)
        """
        gatk PreprocessIntervals \\
            -R $RefFasta \\
            -L ${interval_list} \\
            --bin-length 0 \\
            --padding ${padding} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            -O ${outFileName}
        """
    else
        """
        gatk --java-options ${JAVA_Xmx} ScatterIntervalsByNs \\
            --REFERENCE $RefFasta \\
            --OUTPUT_TYPE ACGT \\
            --OUTPUT ${outFileName}
        """
}

// Splitting interval file in 20(default) files for scattering Mutect2
process 'SplitIntervals' {

    label 'nextNEOpiENV'

    tag "SplitIntervals"

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/SplitIntervals/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    path(IntervalsList) from preprocessIntervalList_out_ch0

    val x from scatter_count
    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    path(
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
    mkdir -p ${tmpDir}

    gatk SplitIntervals \\
        --tmp-dir ${tmpDir} \\
        -R ${RefFasta}  \\
        -scatter ${x} \\
        --interval-merging-rule ALL \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
        -L ${IntervalsList} \\
        -O ${IntervalName}

    """
}


// convert padded interval list to Bed file (used by varscan)
// generate a padded tabix indexed region BED file for strelka
process 'IntervalListToBed' {

    label 'nextNEOpiENV'

    tag 'BedFromIntervalList'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/",
        mode: publishDirMode

    input:
    path(paddedIntervalList) from preprocessIntervalList_out_ch1
    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    path("${paddedIntervalList.baseName}.{bed.gz,bed.gz.tbi}") into (
        RegionsBedToTabix_out_ch0,
        RegionsBedToTabix_out_ch1
    )

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${JAVA_Xmx} IntervalListToBed \\
        -I ${paddedIntervalList} \\
        -O ${paddedIntervalList.baseName}.bed

    bgzip -c ${paddedIntervalList.baseName}.bed > ${paddedIntervalList.baseName}.bed.gz &&
    tabix -p bed ${paddedIntervalList.baseName}.bed.gz
    """
}

// convert scattered padded interval list to Bed file (used by varscan)
process 'ScatteredIntervalListToBed' {

    label 'nextNEOpiENV'

    tag 'ScatteredIntervalListToBed'

    publishDir "$params.outputDir/supplemental/00_prepare_Intervals/SplitIntervals/${IntervalName}",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        val(IntervalName),
        file(IntervalsList)
    ) from SplitIntervals_out_ch0_Name
        .combine(
            SplitIntervals_out_ch0.flatten()
        )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0


    output:
    path(
        "*.bed"
    ) into (
        ScatteredIntervalListToBed_out_ch0,
        ScatteredIntervalListToBed_out_ch1
    )

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    gatk --java-options ${JAVA_Xmx} IntervalListToBed \\
        -I ${IntervalsList} \\
        -O ${IntervalsList.baseName}.bed
    """
}

// FastQC
process FastQC {

    label 'fastqc'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "${params.outputDir}/analyses/${meta.sampleName}/QC/fastqc",
        mode: publishDirMode,
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}


    input:
    tuple val(meta), path(reads) from fastqc_reads_ch


    output:
    tuple val(meta), path("*_fastqc.zip") into ch_fastqc // multiQC

    script:
    def reads_R1_ext = (reads[0].getExtension() == "gz") ? "fastq.gz" : reads[0].getExtension()
    def reads_R1     = meta.sampleName + "_" + meta.sampleType + "_R1." + reads_R1_ext

    // do we have PE reads?
    def reads_R2 = "_missing_"
    if(meta.libType == "PE") {
        def reads_R2_ext = (reads[1].getExtension() == "gz") ? "fastq.gz" : reads[1].getExtension()
        reads_R2     = meta.sampleName + "_" + meta.sampleType + "_R2." + reads_R2_ext
    }
    """
    if [ ! -e ${reads_R1} ]; then
        ln -s ${reads[0]} ${reads_R1}
    fi

    if [ "${reads_R2}" != "_missing_" ] && [ ! -e ${reads_R2} ]; then
        ln -s ${reads[1]} ${reads_R2}
    fi

    fastqc --quiet --threads ${task.cpus} \\
        ${reads_R1} ${reads_R2}
    """
}


// adapter trimming
//
// We first check if we need to trim DNA and RNA or only one of them.
// If it is only one of them we need to combine the trimmed and raw
// reads again

def trim_adapters = false

if (params.trim_adapters || params.trim_adapters_RNAseq) {
    trim_adapters = true
    reads_to_trim = raw_reads_ch

    if (params.trim_adapters && params.trim_adapters_RNAseq) {
        reads_to_trim_ch = raw_reads_ch
        reads_to_keep_ch = Channel.empty()
    }

    if (params.trim_adapters && ! params.trim_adapters_RNAseq) {
        raw_reads_ch.branch {
            DNA: it[0].sampleType != "tumor_RNA"
            RNA: it[0].sampleType == "tumor_RNA"
        }
        .set{ raw_reads_ch }

        reads_to_trim_ch = raw_reads_ch.DNA
        reads_to_keep_ch = raw_reads_ch.RNA
    }

    if (! params.trim_adapters && params.trim_adapters_RNAseq) {
        raw_reads_ch.branch {
            DNA: it[0].sampleType != "tumor_RNA"
            RNA: it[0].sampleType == "tumor_RNA"
        }
        .set{ raw_reads_ch }

        reads_to_trim_ch = raw_reads_ch.RNA
        reads_to_keep_ch = raw_reads_ch.DNA
    }
}

if (trim_adapters) {
    process fastp {

        label 'fastp'

        tag "$meta.sampleName : $meta.sampleType"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/",
            mode: publishDirMode,
            saveAs: {
                filename ->
                    if(filename.indexOf(".json") > 0) {
                        return "QC/fastp/$filename"
                    } else if(filename.indexOf("NO_FILE") >= 0) {
                        return null
                    } else {
                        return  "01_preprocessing/$filename"
                    }
            }

        input:
        tuple val(meta), path(reads) from reads_to_trim_ch

        output:
        tuple val(meta), path("*_trimmed_R{1,2}.fastq.gz") into reads_trimmed_ch, fastqc_reads_trimmed_ch
        tuple val(meta), path("*.json") into ch_fastp // multiQC


        script:
        def reads_R1         = "--in1 " + reads[0]
        def trimmed_reads_R1 = "--out1 " + meta.sampleName + "_" + meta.sampleType + "_trimmed_R1.fastq.gz"

        // do we have PE reads?
        def reads_R2         = ""
        def trimmed_reads_R2 = ""
        if(meta.libType == "PE") {
            reads_R2          = "--in2 " + reads[1]
            trimmed_reads_R2  = "--out2 " + meta.sampleName + "_" + meta.sampleType + "_trimmed_R2.fastq.gz"
        }

        def fastpAdapter = ''
        def adapterSeqFile
        def aseq = false
        def aseqR2 = false
        def afile = false

        if (meta.sampleType.indexOf("DNA") > 0) {
            afile = params.adapterSeqFile
            aseq = params.adapterSeq
            aseqR2 = params.adapterSeqR2
        } else {
             afile = params.adapterSeqFileRNAseq
             aseq = params.adapterSeqRNAseq
             aseqR2 = params.adapterSeqR2RNAseq
        }
        if(afile != false) {
            adapterSeqFile = Channel.fromPath(afile)
            fastpAdapter = "--adapter_fasta " + adapterSeqFile
        } else {
            if(aseq != false) {
                adapterSeq   = Channel.value(aseq)
                fastpAdapter = "--adapter_sequence " + aseq.getVal()

                if(aseqR2 != false && meta.libType == "PE") {
                    adapterSeqR2   = Channel.value(aseqR2)
                    fastpAdapter += " --adapter_sequence_r2 " + adapterSeqR2.getVal()
                }
            }
        }

        """
        fastp --thread ${task.cpus} \\
            ${reads_R1} \\
            ${reads_R2} \\
            ${trimmed_reads_R1} \\
            ${trimmed_reads_R2} \\
            --json ${meta.sampleName}_${meta.sampleType}_fastp.json \\
            ${fastpAdapter} \\
            ${params.fastpOpts}
        """
    }


    // FastQC after adapter trimming
    process FastQC_trimmed {

        label 'fastqc'

        tag "$meta.sampleName : $meta.sampleType"

        publishDir "${params.outputDir}/analyses/${meta.sampleName}/QC/fastqc",
            mode: publishDirMode,
            saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}


        input:
        tuple val(meta), path(reads) from fastqc_reads_trimmed_ch


        output:
        tuple val(meta), path("*_fastqc.zip") into ch_fastqc_trimmed // multiQC

        script:
        def reads_R1 = reads[0]
        def reads_R2 = (meta.libType == "PE") ? reads[1] : ""
        """
        fastqc --quiet --threads ${task.cpus} \\
            ${reads_R1} ${reads_R2}
        """
    }

    // combine trimmed reads ch with reads channel of reads that
    // did not need trimming
    reads_ch = reads_trimmed_ch.mix(reads_to_keep_ch)

} else { // no adapter trimming
    ch_fastqc_trimmed = Channel.empty()
    reads_ch = raw_reads_ch
    ch_fastp = Channel.empty()
}

// get DNA/RNA reads
reads_ch.branch {
    DNA: it[0].sampleType != "tumor_RNA"
    RNA: it[0].sampleType == "tumor_RNA"
}
.set{ reads_ch }

(reads_BAM_ch, reads_uBAM_ch, reads_mixcr_DNA_ch, dummy_ch) = reads_ch.DNA.into(4)
(reads_tumor_optitype_ch, reads_tumor_hlahd_RNA_ch, reads_tumor_neofuse_ch, reads_tumor_mixcr_RNA_ch) = reads_ch.RNA.into(4)

reads_mixcr_ch = reads_mixcr_DNA_ch.mix(reads_tumor_mixcr_RNA_ch)

// setup dummy RNA channels
dummy_ch.filter{
    it[0].have_RNA == false
}.map{
    it ->
        return [it[0], []]
}
.into{ no_RNA_0; no_RNA_1 }

/////// start processing reads ///////

// make uBAM
process 'make_uBAM' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/01_preprocessing/",
        mode: publishDirMode

    input:
    tuple val(meta), path(reads) from reads_uBAM_ch

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple val(meta), path(ubam) into uBAM_out_ch0

    script:
    ubam = meta.sampleName + "_" + meta.sampleType + "_unaligned.bam"
    def read_group = meta.sampleName + "_" + meta.sampleType.replaceAll("_DNA", "")
    def reads_in = "-F1 " + reads[0]
    reads_in += (meta.libType == "PE") ? " -F2 " + reads[1] : ""
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'

    """
    mkdir -p ${tmpDir}
    gatk --java-options ${java_opts} FastqToSam \\
        --TMP_DIR ${tmpDir} \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} \\
        ${reads_in} \\
        --READ_GROUP_NAME ${read_group} \\
        --SAMPLE_NAME ${read_group} \\
        --LIBRARY_NAME ${read_group} \\
        --PLATFORM ILLUMINA \\
        -O ${ubam}
    """
}

// Aligning reads to reference, sort and index; create BAMs
process 'Bwa' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/02_alignments/",
        mode: publishDirMode

    input:
    tuple val(meta), path(reads) from reads_BAM_ch

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
        path(BwaRef)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.BwaRef ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple val(meta), path(bam) into BWA_out_ch0

    script:
    bam = meta.sampleName + "_" + meta.sampleType + "_aligned.bam"
    def read_group = meta.sampleName + "_" + meta.sampleType.replaceAll("_DNA", "")

    def sort_threads = (task.cpus.compareTo(8) == 1) ? 8 : task.cpus
    def SB_sort_mem =  Math.max((task.memory.toGiga() - 4), 1) + "G"
    """
    bwa mem \\
        -R "@RG\\tID:${read_group}\\tLB:${read_group}\\tSM:${read_group}\\tPL:ILLUMINA" \\
        -M ${RefFasta} \\
        -t ${task.cpus} \\
        -Y \\
        ${reads.join(" ")} | \\
    samtools view -@2 -Shbu - | \\
    sambamba sort \\
        --sort-picard \\
        --tmpdir=${tmpDir} \\
        -m ${SB_sort_mem} \\
        -l 6 \\
        -t ${sort_threads} \\
        -o ${bam} \\
        /dev/stdin
    """
}

// merge alinged BAM and uBAM
process 'merge_uBAM_BAM' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/02_alignments/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple val(meta), path(bam), path(ubam) from BWA_out_ch0.join(uBAM_out_ch0, by: [0])

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple val(meta), path("${procSampleName}_aligned_uBAM_merged.bam") into uBAM_BAM_out_ch

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType

    def paired_run = (meta.libType == 'SE') ? 'false' : 'true'
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${tmpDir}

    gatk --java-options ${java_opts} MergeBamAlignment \\
        --TMP_DIR ${tmpDir} \\
        --VALIDATION_STRINGENCY SILENT \\
        --EXPECTED_ORIENTATIONS FR \\
        --ATTRIBUTES_TO_RETAIN X0 \\
        --REFERENCE_SEQUENCE ${RefFasta} \\
        --PAIRED_RUN ${paired_run} \\
        --SORT_ORDER "queryname" \\
        --IS_BISULFITE_SEQUENCE false \\
        --ALIGNED_READS_ONLY false \\
        --CLIP_ADAPTERS false \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRamMerge} \\
        --ADD_MATE_CIGAR true \\
        --MAX_INSERTIONS_OR_DELETIONS -1 \\
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \\
        --ALIGNER_PROPER_PAIR_FLAGS true \\
        --UNMAP_CONTAMINANT_READS true \\
        --ALIGNED_BAM ${bam} \\
        --UNMAPPED_BAM ${ubam} \\
        --OUTPUT ${procSampleName}_aligned_uBAM_merged.bam
    """
}


// Mark duplicates with sambamba
process 'MarkDuplicates' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/02_alignments/",
        mode: publishDirMode

    input:
    tuple val(meta), path(bam) from uBAM_BAM_out_ch

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path(bam_out)
    ) into (
        MarkDuplicates_out_ch0,
        MarkDuplicates_out_ch1,
        MarkDuplicates_out_ch2,
        MarkDuplicates_out_ch3,
        MarkDuplicates_out_ch4
    )

    script:
    def procSampleName = meta.sampleName + "_" + meta.sampleType
    def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 8) / task.cpus)), 1)
    def JAVA_Xmx = '-Xmx4G'
    bam_out = [procSampleName + "_aligned_sort_mkdp.bam", procSampleName + "_aligned_sort_mkdp.bai"]
    """
    mkdir -p ${tmpDir}
    sambamba markdup \\
        -t ${task.cpus} \\
        --tmpdir ${tmpDir} \\
        --hash-table-size=${params.SB_hash_table_size } \\
        --overflow-list-size=${params.SB_overflow_list_size} \\
        --io-buffer-size=${params.SB_io_buffer_size} \\
        ${bam} \\
        /dev/stdout | \\
    samtools sort \\
        -@${task.cpus} \\
        -m ${STperThreadMem}G \\
        -O BAM \\
        -l 0 \\
        /dev/stdin | \\
    gatk --java-options ${JAVA_Xmx} SetNmMdAndUqTags \\
        --TMP_DIR ${tmpDir} \\
        -R ${RefFasta} \\
        -I /dev/stdin \\
        -O ${procSampleName}_aligned_sort_mkdp.bam \\
        --CREATE_INDEX true \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} \\
        --VALIDATION_STRINGENCY LENIENT

    """
}

// prepare channel for mhc_extract -> hlad-hd, optitype
 MarkDuplicates_out_ch3.filter {
                                    it[0].sampleType == "tumor_DNA"
                                }
                                .set { MarkDuplicatesTumor_out_ch0 }

// spilt T/N and remove differing/unused info from meta for joining
// this prepares T/N channel for CNVkit

MarkDuplicates_out_ch4.branch {
        meta_ori, bam ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor_DNA"
                meta.remove('sampleType')
                return [meta, bam]

            normal: meta.sampleType == "normal_DNA"
                meta.remove('sampleType')
                return [meta, bam]
}.set{ MarkDuplicates_out_ch4 }

MarkDuplicates_out_CNVkit_ch0 = MarkDuplicates_out_ch4.tumor.join(MarkDuplicates_out_ch4.normal, by:[0])

if(params.WES) {
    // Generate HS metrics using picard
    process 'alignmentMetrics' {

        label 'nextNEOpiENV'

        tag "$meta.sampleName : $meta.sampleType"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/QC/alignments/",
            mode: publishDirMode

        input:
        tuple val(meta), path(bam) from MarkDuplicates_out_ch0

        tuple(
            path(RefFasta),
            path(RefIdx)
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx ]
        )

        path(BaitIntervalsList) from BaitsBedToIntervalList_out_ch0
        path(IntervalsList) from RegionsBedToIntervalList_out_ch1

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple val(meta), path("${procSampleName}.*.txt") into alignmentMetrics_ch // multiQC

        script:
        procSampleName = meta.sampleName + "_" + meta.sampleType
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
        def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
        """
        mkdir -p ${tmpDir}
        gatk --java-options ${java_opts} CollectHsMetrics \\
            --TMP_DIR ${tmpDir} \\
            --INPUT ${bam[0]} \\
            --OUTPUT ${procSampleName}.HS.metrics.txt \\
            -R ${RefFasta} \\
            --BAIT_INTERVALS ${BaitIntervalsList} \\
            --TARGET_INTERVALS ${IntervalsList} \\
            --PER_TARGET_COVERAGE ${procSampleName}.perTarget.coverage.txt && \\
        gatk --java-options ${java_opts} CollectAlignmentSummaryMetrics \\
            --TMP_DIR ${tmpDir} \\
            --INPUT ${bam[0]} \\
            --OUTPUT ${procSampleName}.AS.metrics.txt \\
            -R ${RefFasta} &&
        samtools flagstat -@${task.cpus} ${bam[0]} > ${procSampleName}.flagstat.txt
        """
    }
} else {
    // bogus channel for multiqc
    alignmentMetrics_ch = Channel.empty()
}


/*
 BaseRecalibrator (GATK4): generates recalibration table for Base Quality Score
 Recalibration (BQSR)
 ApplyBQSR (GATK4): apply BQSR table to reads
*/
process 'scatterBaseRecalGATK4' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    input:
    tuple(
        val(meta),
        path(bam),
        path(intervals)
    ) from MarkDuplicates_out_ch1
        .combine(
            SplitIntervals_out_scatterBaseRecalTumorGATK4_ch.flatten()
        )

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    tuple(
        path(MillsGold),
        path(MillsGoldIdx),
        path(DBSNP),
        path(DBSNPIdx),
        path(KnownIndels),
        path(KnownIndelsIdx)
    ) from Channel.value(
        [ database.MillsGold,
          database.MillsGoldIdx,
          database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple val(meta), path("${procSampleName}_${intervals}_bqsr.table") into scatterBaseRecalGATK4_out_ch0


    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${tmpDir}
    gatk  --java-options ${JAVA_Xmx} BaseRecalibrator \\
        --tmp-dir ${tmpDir} \\
        -I ${bam[0]} \\
        -R ${RefFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_bqsr.table \\
        --known-sites ${DBSNP} \\
        --known-sites ${KnownIndels} \\
        --known-sites ${MillsGold}
    """
}

// gather scattered bqsr tables
process 'gatherGATK4scsatteredBQSRtables' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/03_baserecalibration/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple val(meta), path(bqsr_table) from scatterBaseRecalGATK4_out_ch0.groupTuple(by: [0])
    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple val(meta), path("${procSampleName}_bqsr.table") into gatherBQSRtables_out_ch0

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType

    """
    mkdir -p ${tmpDir}

    gatk GatherBQSRReports \\
        -I ${bqsr_table.join(" -I ")} \\
        -O ${procSampleName}_bqsr.table
    """
}


// ApplyBQSR (GATK4): apply BQSR table to reads
process 'scatterGATK4applyBQSRS' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    input:
    tuple(
        meta,
        path(bam),
        path(bqsr_table),
        path(intervals)
    ) from MarkDuplicates_out_ch2
        .join(gatherBQSRtables_out_ch0, by: [0])
        .combine(
            SplitIntervals_out_scatterTumorGATK4applyBQSRS_ch.flatten()
        )

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    tuple(
        path(MillsGold),
        path(MillsGoldIdx),
        path(DBSNP),
        path(DBSNPIdx),
        path(KnownIndels),
        path(KnownIndelsIdx)
    ) from Channel.value(
        [ database.MillsGold,
          database.MillsGoldIdx,
          database.DBSNP,
          database.DBSNPIdx,
          database.KnownIndels,
          database.KnownIndelsIdx ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path(bam_out)
    ) into scatterGATK4applyBQSRS_out_GatherRecalBamFiles_ch0


    script:
    def procSampleName = meta.sampleName + "_" + meta.sampleType
    bam_out = [ procSampleName + "_" + intervals + "_recal4.bam",
                procSampleName + "_" + intervals + "_recal4.bai"]
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${tmpDir}
    gatk ApplyBQSR \\
        --java-options ${JAVA_Xmx} \\
        --tmp-dir ${tmpDir} \\
        -I ${bam[0]} \\
        -R ${RefFasta} \\
        -L ${intervals} \\
        -O ${procSampleName}_${intervals}_recal4.bam \\
        --bqsr-recal-file ${bqsr_table}
    """
}

process 'GatherRecalBamFiles' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/03_baserecalibration/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(bam),
        path(bai)
    ) from scatterGATK4applyBQSRS_out_GatherRecalBamFiles_ch0
        .toSortedList({a, b -> a[1][0].baseName <=> b[1][0].baseName})
        .flatten()
        .collate(3)
        .groupTuple(by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${procSampleName}_recalibrated.{bam,bam.bai}")
    ) into (
        BaseRecalGATK4_out_ch0,
        BaseRecalGATK4_out_ch1,
        GatherRecalBamFiles_out_IndelRealignerIntervals_ch0
    )

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType
    def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 4) / task.cpus)), 1)
    def JAVA_Xmx = "4G"
    def java_opts = '"-Xmx' + JAVA_Xmx + ' -XX:ParallelGCThreads=2"'
    """
    mkdir -p ${tmpDir}

    rm -f ${procSampleName}_gather.fifo
    mkfifo ${procSampleName}_gather.fifo
    gatk --java-options ${java_opts} GatherBamFiles \\
        --TMP_DIR ${tmpDir} \\
        -I ${bam.join(" -I ")} \\
        -O ${procSampleName}_gather.fifo \\
        --CREATE_INDEX false \\
        --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} &
    samtools sort \\
        -@${task.cpus} \\
        -m ${STperThreadMem}G \\
        -o ${procSampleName}_recalibrated.bam ${procSampleName}_gather.fifo
    samtools index -@${task.cpus} ${procSampleName}_recalibrated.bam
    rm -f ${procSampleName}_gather.fifo
    """
}


// GetPileupSummaries (GATK4): tabulates pileup metrics for inferring contamination
process 'GetPileup' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName : $meta.sampleType"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/mutect2/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(GnomAD),
        path(GnomADIdx)
    ) from Channel.value(
        [ database.GnomAD,
          database.GnomADIdx ]
    )

    tuple(
        val(meta),
        path(bam),
        path(IntervalsList)
    ) from BaseRecalGATK4_out_ch0
        .combine(
            preprocessIntervalList_out_ch3
        )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple val(meta), path("${procSampleName}_pileup.table") into GetPileup_out_ch0

    script:
    procSampleName = meta.sampleName + "_" + meta.sampleType
    """
    mkdir -p ${tmpDir}

    gatk GetPileupSummaries \\
        --tmp-dir ${tmpDir} \\
        -I ${bam[0]} \\
        -O ${procSampleName}_pileup.table \\
        -L ${IntervalsList} \\
        --variant ${GnomAD}
    """
}

(BaseRecalNormal_out_ch0, BaseRecalGATK4_out) = BaseRecalGATK4_out_ch1.into(2)

BaseRecalNormal_out_ch0.filter {
                                it[0].sampleType == "normal_DNA"
                            }
                            .map {
                                meta_ori, bam ->
                                    def meta = meta_ori.clone()
                                    meta.remove('sampleType')
                                    return [meta, bam]
                            }
                            .set { BaseRecalNormal_out_ch0 }

BaseRecalGATK4_out.branch {
        meta_ori, bam ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor_DNA"
                meta.remove('sampleType')
                return [meta, bam]

            normal: meta.sampleType == "normal_DNA"
                meta.remove('sampleType')
                return [meta, bam]
}.set{ BaseRecalGATK4_out }

BaseRecalGATK4_out = BaseRecalGATK4_out.tumor.join(BaseRecalGATK4_out.normal, by: [0])

if (have_GATK3) {
    (
        BaseRecalGATK4_out_Mutect2_ch0,
        BaseRecalGATK4_out_MantaSomaticIndels_ch0,
        BaseRecalGATK4_out_StrelkaSomatic_ch0,
        BaseRecalGATK4_out_MutationalBurden_ch0,
        BaseRecalGATK4_out_MutationalBurden_ch1
    ) = BaseRecalGATK4_out.into(5)
} else {
    (
        BaseRecalGATK4_out_Mutect2_ch0,
        BaseRecalGATK4_out_MantaSomaticIndels_ch0,
        BaseRecalGATK4_out_StrelkaSomatic_ch0,
        BaseRecalGATK4_out_MutationalBurden_ch0,
        BaseRecalGATK4_out_MutationalBurden_ch1,
        BaseRecalGATK4_out
    ) = BaseRecalGATK4_out.into(6)
}

/*
    MUTECT2
    Call somatic SNPs and indels via local re-assembly of haplotypes; tumor sample
    and matched normal sample
*/
process 'Mutect2' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/mutect2/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
        path(gnomADfull),
        path(gnomADfullIdx)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          database.GnomADfull,
          database.GnomADfullIdx ]
    )
    path pon from pon_file

    tuple(
        val(meta),
        path(Tumorbam),
        path(Normalbam),
        path(intervals)
    ) from BaseRecalGATK4_out_Mutect2_ch0
        .combine(
            SplitIntervals_out_ch2.flatten()
        )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_${intervals}.vcf.gz"),
        path("${meta.sampleName}_${intervals}.vcf.gz.tbi"),
        path("${meta.sampleName}_${intervals}.vcf.gz.stats"),
        path("${meta.sampleName}_${intervals}-f1r2.tar.gz")
    ) into Mutect2_out_ch0


    script:
    def tumorName  = meta.sampleName + "_tumor"
    def normalName = meta.sampleName + "_normal"

    def panel_of_normals = (pon.name != 'NO_FILE') ? "--panel-of-normals $pon" : ""
    def mk_pon_idx = (pon.name != 'NO_FILE') ? "tabix -f $pon" : ""
    """
    mkdir -p ${tmpDir}

    ${mk_pon_idx}

    gatk Mutect2 \\
        --tmp-dir ${tmpDir} \\
        -R ${RefFasta} \\
        -I ${Tumorbam[0]} -tumor ${tumorName} \\
        -I ${Normalbam[0]} -normal ${normalName} \\
        --germline-resource ${gnomADfull} \\
        ${panel_of_normals} \\
        -L ${intervals} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --f1r2-tar-gz ${meta.sampleName}_${intervals}-f1r2.tar.gz \\
        -O ${meta.sampleName}_${intervals}.vcf.gz
    """
}

// Merge scattered Mutect2 vcfs
process 'gatherMutect2VCFs' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/mutect2/",
        mode: publishDirMode,
        saveAs: {
            filename ->
                if(filename.indexOf("_read-orientation-model.tar.gz") > 0 && params.fullOutput) {
                    return "processing/$filename"
                } else if(filename.indexOf("_read-orientation-model.tar.gz") > 0 && ! params.fullOutput) {
                    return null
                } else {
                    return "raw/$filename"
                }
        }


    input:
    tuple(
        val(meta),
        path(vcf),
        path(tbi),
        path(stats),
        path(f1r2_tar_gz)
    ) from Mutect2_out_ch0
        .groupTuple(by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_mutect2_raw.{vcf.gz,vcf.gz.tbi}"),
        path("${meta.sampleName}_mutect2_raw.vcf.gz.stats"),
        path("${meta.sampleName}_read-orientation-model.tar.gz")
    ) into gatherMutect2VCFs_out_ch0

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${tmpDir}

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${vcf.join(" -I ")} \\
        -O ${meta.sampleName}_mutect2_raw.vcf.gz

    gatk MergeMutectStats \\
        --tmp-dir ${tmpDir} \\
        --stats ${stats.join(" --stats ")} \\
        -O ${meta.sampleName}_mutect2_raw.vcf.gz.stats

    gatk LearnReadOrientationModel \\
        --tmp-dir ${tmpDir} \\
        -I ${f1r2_tar_gz.join(" -I ")} \\
        -O ${meta.sampleName}_read-orientation-model.tar.gz
    """
}

GetPileup_out_ch0.branch {
        meta_ori, pileup ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor_DNA"
                meta.remove('sampleType')
                return [meta, pileup]

            normal: meta.sampleType == "normal_DNA"
                meta.remove('sampleType')
                return [meta, pileup]
}.set{ GetPileup_out }

GetPileup_out = GetPileup_out.tumor.join(GetPileup_out.normal, by: [0])

/*
CalculateContamination (GATK4): calculate fraction of reads coming from
cross-sample contamination
FilterMutectCalls (GATK4): filter somatic SNVs and indels
FilterByOrientationBias (GATK4): filter variant calls using orientation bias
SelectVariants (GATK4): select subset of variants from a larger callset
VariantFiltration (GATK4): filter calls based on INFO and FORMAT annotations
*/
process 'FilterMutect2' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/mutect2/",
        mode: publishDirMode

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    tuple(
        val(meta),
        path(pileupTumor),
        path(pileupNormal),
        path(vcf),
        path(vcfStats),
        path(f1r2_tar_gz)
    ) from GetPileup_out
        .join(gatherMutect2VCFs_out_ch0, by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        val("mutect2"),
        path("${meta.sampleName}_mutect2_final.{vcf.gz,vcf.gz.tbi}")
    ) into (
        FilterMutect2_out_ch0,
        FilterMutect2_out_ch1
    )

    script:
    """
    mkdir -p ${tmpDir}

    gatk CalculateContamination \\
        --tmp-dir ${tmpDir} \\
        -I ${pileupTumor} \\
        --matched-normal ${pileupNormal} \\
        -O ${meta.sampleName}_cont.table && \\
    gatk FilterMutectCalls \\
        --tmp-dir ${tmpDir} \\
        -R ${RefFasta} \\
        -V ${vcf[0]} \\
        --contamination-table ${meta.sampleName}_cont.table \\
        --ob-priors ${f1r2_tar_gz} \\
        -O ${meta.sampleName}_oncefiltered.vcf.gz && \\
    gatk SelectVariants \\
        --tmp-dir ${tmpDir} \\
        --variant ${meta.sampleName}_oncefiltered.vcf.gz \\
        -R ${RefFasta} \\
        --exclude-filtered true \\
        --select 'vc.getGenotype(\"${meta.sampleName}_tumor\").getAD().1 >= ${params.minAD}' \\
        --output ${meta.sampleName}_mutect2_final.vcf.gz
    """
}


// HaploTypeCaller
/*
    Call germline SNPs and indels via local re-assembly of haplotypes; normal sample
    germline variants are needed for generating phased vcfs for pVACtools
*/
process 'HaploTypeCaller' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/haplotypecaller/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
        path(DBSNP),
        path(DBSNPIdx)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          database.DBSNP,
          database.DBSNPIdx ]
    )

    tuple(
        val(meta),
        path(Normalbam),
        path(intervals)
    ) from BaseRecalNormal_out_ch0
        .combine(
            SplitIntervals_out_ch5.flatten()
        )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_germline_${intervals}.{vcf.gz,vcf.gz.tbi}"),
        path(Normalbam)
    ) into (
        HaploTypeCaller_out_ch0
    )


    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    mkdir -p ${tmpDir}

    gatk --java-options ${JAVA_Xmx} HaplotypeCaller \\
        --tmp-dir ${tmpDir} \\
        -R ${RefFasta} \\
        -I ${Normalbam[0]} \\
        -L ${intervals} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --dbsnp ${DBSNP} \\
        -O ${meta.sampleName}_germline_${intervals}.vcf.gz
    """
}


/*
    Run a Convolutional Neural Net to filter annotated germline variants; normal sample
    germline variants are needed for generating phased vcfs for pVACtools
*/
process 'CNNScoreVariants' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/haplotypecaller/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    tuple(
        val(meta),
        path(raw_germline_vcf),
        path(Normalbam)
    ) from HaploTypeCaller_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${raw_germline_vcf[0].baseName}_CNNScored.vcf.gz"),
        path("${raw_germline_vcf[0].baseName}_CNNScored.vcf.gz.tbi")
    ) into CNNScoreVariants_out_ch0


    script:
    """
    mkdir -p ${tmpDir}

    gatk CNNScoreVariants \\
        --tmp-dir ${tmpDir} \\
        -R ${RefFasta} \\
        -I ${Normalbam[0]} \\
        -V ${raw_germline_vcf[0]} \\
        -tensor-type read_tensor \\
        --inter-op-threads ${task.cpus} \\
        --intra-op-threads ${task.cpus} \\
        --transfer-batch-size ${params.transferBatchSize} \\
        --inference-batch-size ${params.inferenceBatchSize} \\
        -O ${raw_germline_vcf[0].baseName}_CNNScored.vcf.gz
    """
}


// Merge scattered filtered germline vcfs
process 'MergeHaploTypeCallerGermlineVCF' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/haplotypecaller/raw/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(filtered_germline_vcf),
        path(filtered_germline_vcf_tbi)
    ) from CNNScoreVariants_out_ch0
        .groupTuple(by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_germline_CNNscored.{vcf.gz,vcf.gz.tbi}")
    ) into MergeHaploTypeCallerGermlineVCF_out_ch0

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${tmpDir}

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${filtered_germline_vcf.join(" -I ")} \\
        -O ${meta.sampleName}_germline_CNNscored.vcf.gz
    """
}

/*
    Apply a Convolutional Neural Net to filter annotated germline variants; normal sample
    germline variants are needed for generating phased vcfs for pVACtools
*/
process 'FilterGermlineVariantTranches' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/haplotypecaller/",
        mode: publishDirMode

    input:
    tuple(
        path(MillsGold),
        path(MillsGoldIdx),
        path(HapMap),
        path(HapMapIdx),
        path(hcSNPS1000G),
        path(hcSNPS1000GIdx)
    ) from Channel.value(
        [ database.MillsGold,
          database.MillsGoldIdx,
          database.HapMap,
          database.HapMapIdx,
          database.hcSNPS1000G,
          database.hcSNPS1000GIdx ]
    )

    tuple(
        val(meta),
        path(scored_germline_vcf)
    ) from MergeHaploTypeCallerGermlineVCF_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${scored_germline_vcf[0].simpleName}_Filtered.{vcf.gz,vcf.gz.tbi}")
    ) into FilterGermlineVariantTranches_out_ch0


    script:
    """
    mkdir -p ${tmpDir}

    gatk FilterVariantTranches \\
        --tmp-dir ${tmpDir} \\
        -V ${scored_germline_vcf[0]} \\
        --resource ${hcSNPS1000G} \\
        --resource ${HapMap} \\
        --resource ${MillsGold} \\
        --info-key CNN_2D \\
        --snp-tranche 99.95 \\
        --indel-tranche 99.4 \\
        --invalidate-previous-filters \\
        -O ${scored_germline_vcf[0].simpleName}_Filtered.vcf.gz
    """
}


// END HTC


/*
 RealignerTargetCreator (GATK3): define intervals to target for local realignment
 IndelRealigner (GATK3): perform local realignment of reads around indels
*/
if (have_GATK3) {
    process 'IndelRealignerIntervals' {

        label 'GATK3'

        tag "$meta.sampleName : $meta.sampleType"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/03_realignment/processing/",
            mode: publishDirMode,
            enabled: params.fullOutput

        input:
        tuple(
            val(meta),
            path(bam)
        ) from GatherRecalBamFiles_out_IndelRealignerIntervals_ch0

        tuple(
            path(RefFasta),
            path(RefIdx),
            path(RefDict)
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx,
            reference.RefDict ]
        )

        tuple(
            path(KnownIndels),
            path(KnownIndelsIdx),
            path(MillsGold),
            path(MillsGoldIdx)
        ) from Channel.value(
            [ database.KnownIndels,
            database.KnownIndelsIdx,
            database.MillsGold,
            database.MillsGoldIdx ]
        )

        each path(interval) from SplitIntervals_out_ch3.flatten()

        output:
        tuple(
            val(meta),
            path(bam_out)
        ) into IndelRealignerIntervals_out_GatherRealignedBamFiles_ch0

        script:
        def procSampleName = meta.sampleName + "_" + meta.sampleType
        bam_out = [ procSampleName + "_recalibrated_realign_" + interval + ".bam",
                    procSampleName + "_recalibrated_realign_" + interval + ".bai"]
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"

        """
        mkdir -p ${tmpDir}

        $JAVA8 ${JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${tmpDir} -jar $GATK3 \\
            -T RealignerTargetCreator \\
            --known ${MillsGold} \\
            --known ${KnownIndels} \\
            -R ${RefFasta} \\
            -L ${interval} \\
            -I ${bam[0]} \\
            -o ${interval}_target.list \\
            -nt ${task.cpus} && \\
        $JAVA8 -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${tmpDir} -jar $GATK3 \\
            -T IndelRealigner \\
            -R ${RefFasta} \\
            -L ${interval} \\
            -I ${bam[0]} \\
            -targetIntervals ${interval}_target.list \\
            -known ${KnownIndels} \\
            -known ${MillsGold} \\
            -nWayOut _realign_${interval}.bam && \\
        rm ${interval}_target.list
        """
    }

    process 'GatherRealignedBamFiles' {

        label 'nextNEOpiENV'

        tag "$meta.sampleName : $meta.sampleType"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/03_realignment/",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(bam),
            path(bai)
        ) from IndelRealignerIntervals_out_GatherRealignedBamFiles_ch0
            .toSortedList({a, b -> a[1][0].baseName <=> b[1][0].baseName})
            .flatten()
            .collate(3)
            .groupTuple(by: [0])

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path("${procSampleName}_recalibrated_realign.{bam,bam.bai}")
        ) into GatherRealignedBamFiles_out_ch

        script:
        procSampleName = meta.sampleName + "_" + meta.sampleType
        def STperThreadMem = (int) Math.max(((int) Math.floor((task.memory.toGiga() - 4) / task.cpus)), 1)
        def JAVA_Xmx = "4G"
        def java_opts = '"-Xmx' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
        """
        mkdir -p ${tmpDir}

        rm -f ${procSampleName}_gather.fifo
        mkfifo ${procSampleName}_gather.fifo
        gatk --java-options ${java_opts} GatherBamFiles \\
            --TMP_DIR ${tmpDir} \\
            -I ${bam.join(" -I ")} \\
            -O ${procSampleName}_gather.fifo \\
            --CREATE_INDEX false \\
            --MAX_RECORDS_IN_RAM ${params.maxRecordsInRam} &
        samtools sort \\
            -@${task.cpus} \\
            -m ${STperThreadMem}G \\
            -o ${procSampleName}_recalibrated_realign.bam ${procSampleName}_gather.fifo
        samtools index -@${task.cpus}  ${procSampleName}_recalibrated_realign.bam
        rm -f ${procSampleName}_gather.fifo
        """
    }

    (
        GatherRealignedBamFiles_out_AlleleCounter_ch0,
        GatherRealignedBamFiles_out_Mpileup4ControFREEC_ch0,
        GatherRealignedBamFiles_out_ch
    ) = GatherRealignedBamFiles_out_ch.into(3)

    GatherRealignedBamFiles_out_ch.branch {
            meta_ori, bam ->
                def meta = meta_ori.clone()
                tumor : meta.sampleType == "tumor_DNA"
                    meta.remove('sampleType')
                    return [meta, bam]

                normal: meta.sampleType == "normal_DNA"
                    meta.remove('sampleType')
                    return [meta, bam]
    }.set{ GatherRealignedBamFiles_out_ch }

    (
        GatherRealignedBamFilesTumor_out_FilterVarscan_ch0,
        GatherRealignedBamFilesTumor_out_mkPhasedVCF_ch0,
        GatherRealignedBamFilesTumor
    ) = GatherRealignedBamFiles_out_ch.tumor.into(3)


    RealignedBamFiles = GatherRealignedBamFilesTumor.join(GatherRealignedBamFiles_out_ch.normal, by: [0])
    (
        VarscanBAMfiles_ch,
        GatherRealignedBamFiles_out_Mutect1scattered_ch0,
        GatherRealignedBamFiles_out_Sequenza_ch0
    ) = RealignedBamFiles.into(3)

} else {

    log.info "INFO: GATK3 not installed! Can not generate indel realigned BAMs for varscan and mutect1\n"

    (
        VarscanBAMfiles_ch,
        GatherRealignedBamFiles_out_Mutect1scattered_ch0,
        GatherRealignedBamFiles_out_Sequenza_ch0,
        BaseRecalGATK4_out
    ) = BaseRecalGATK4_out.into(4)

    GatherRealignedBamFilesTumor_out_FilterVarscan_ch0 = BaseRecalGATK4_out
                                                            .map{ meta,
                                                                  recalTumorBAM, recalNormalBAM -> tuple(
                                                                      meta, recalTumorBAM
                                                                  )
                                                            }
    GatherRealignedBamFiles_out_AlleleCounter_ch0 = GatherRecalBamFiles_out_IndelRealignerIntervals_ch0

} // END if have GATK3

process 'VarscanSomaticScattered' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/varscan/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(Tumorbam),
        path(Normalbam),
        path(intervals)
    ) from VarscanBAMfiles_ch // GatherRealignedBamFiles_out_VarscanSomaticScattered_ch0
        .combine(
            ScatteredIntervalListToBed_out_ch0.flatten()
        )

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_${intervals}_varscan.snp.vcf"),
        path("${meta.sampleName}_${intervals}_varscan.indel.vcf")
    ) into VarscanSomaticScattered_out_ch0

    script:
    // awk filters at the end needed, found at least one occurence of "W" in Ref field of
    // varscan vcf (? wtf). Non ACGT seems to cause MergeVCF (picard) crashing
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    rm -f ${meta.sampleName}_${intervals}_mpileup.fifo
    mkfifo ${meta.sampleName}_${intervals}_mpileup.fifo
    samtools mpileup \\
        -q 1 \\
        -f ${RefFasta} \\
        -l ${intervals} \\
        ${Normalbam[0]} ${Tumorbam[0]} > ${meta.sampleName}_${intervals}_mpileup.fifo &
    varscan ${JAVA_Xmx} somatic \\
        ${meta.sampleName}_${intervals}_mpileup.fifo \\
        ${meta.sampleName}_${intervals}_varscan_tmp \\
        --output-vcf 1 \\
        --mpileup 1 \\
        --min-coverage ${params.min_cov} \\
        --min-coverage-normal ${params.min_cov_normal} \\
        --min-coverage-tumor ${params.min_cov_tumor} \\
        --min-freq-for-hom ${params.min_freq_for_hom} \\
        --tumor-purity 0.95 \\
        --p-value ${params.somatic_pvalue} \\
        --somatic-p-value ${params.somatic_somaticpvalue} \\
        --strand-filter ${params.strand_filter} && \\
    rm -f ${meta.sampleName}_${intervals}_mpileup.fifo && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]/) { print } } else { print } }' \\
        ${meta.sampleName}_${intervals}_varscan_tmp.snp.vcf \\
        > ${meta.sampleName}_${intervals}_varscan.snp.vcf && \\
    awk '{OFS=FS="\t"} { if(\$0 !~ /^#/) { if (\$4 ~ /[ACGT]+/) { print } } else { print } }' \\
        ${meta.sampleName}_${intervals}_varscan_tmp.indel.vcf \\
        > ${meta.sampleName}_${intervals}_varscan.indel.vcf

    rm -f ${meta.sampleName}_${intervals}_varscan_tmp.*
    """
}

// Merge scattered Varscan vcfs
process 'gatherVarscanVCFs' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/varscan/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    tuple(
        val(meta),
        path(snp_vcf),
        path(indel_vcf)
    ) from VarscanSomaticScattered_out_ch0
        .toSortedList({a, b -> a[1].baseName <=> b[1].baseName})
        .flatten()
        .collate(3)
        .groupTuple(by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_varscan.snp.vcf"),
        path("${meta.sampleName}_varscan.indel.vcf")
    ) into gatherVarscanVCFs_out_ch0

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${tmpDir}

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${snp_vcf.join(" -I ")} \\
        -O ${meta.sampleName}_varscan.snp.vcf \\
        --SEQUENCE_DICTIONARY ${RefDict}

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${indel_vcf.join(" -I ")} \\
        -O ${meta.sampleName}_varscan.indel.vcf \\
        --SEQUENCE_DICTIONARY ${RefDict}

    """
}

// Filter variants by somatic status and confidences
process 'ProcessVarscan' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/varscan/raw/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(snp),
        path(indel)
    ) from gatherVarscanVCFs_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_varscan.snp.Somatic.vcf"),
        path("${meta.sampleName}_varscan.snp.Somatic.hc.vcf"),
        path("${meta.sampleName}_varscan.snp.LOH.vcf"),
        path("${meta.sampleName}_varscan.snp.LOH.hc.vcf"),
        path("${meta.sampleName}_varscan.snp.Germline.vcf"),
        path("${meta.sampleName}_varscan.snp.Germline.hc.vcf")
    ) into ProcessVarscanSNP_out_ch0

    tuple(
        val(meta),
        path("${meta.sampleName}_varscan.indel.Somatic.vcf"),
        path("${meta.sampleName}_varscan.indel.Somatic.hc.vcf"),
        path("${meta.sampleName}_varscan.indel.LOH.vcf"),
        path("${meta.sampleName}_varscan.indel.LOH.hc.vcf"),
        path("${meta.sampleName}_varscan.indel.Germline.vcf"),
        path("${meta.sampleName}_varscan.indel.Germline.hc.vcf")
    ) into ProcessVarscanIndel_out_ch0

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    varscan ${JAVA_Xmx} processSomatic \\
        ${snp} \\
        --min-tumor-freq ${params.min_tumor_freq} \\
        --max-normal-freq ${params.max_normal_freq} \\
        --p-value ${params.processSomatic_pvalue} && \\
    varscan ${JAVA_Xmx} processSomatic \\
        ${indel} \\
        --min-tumor-freq ${params.min_tumor_freq} \\
        --max-normal-freq ${params.max_normal_freq} \\
        --p-value ${params.processSomatic_pvalue}
    """
}

/*
    AWK-script: calcualtes start-end position of variant
    Bamreadcount: generate metrics at single nucleotide positions for filtering
    fpfilter (Varscan): apply false-positive filter to variants
*/
process 'FilterVarscan' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/varscan/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(bam),
        path(snpSomatic),
        path(snpSomaticHc),
        path(snpLOH),
        path(snpLOHhc),
        path(snpGerm),
        path(snpGemHc),
        path(indelSomatic),
        path(indelSomaticHc),
        path(indelLOH),
        path(indelLOHhc),
        path(indelGerm),
        path(indelGemHc)
    ) from GatherRealignedBamFilesTumor_out_FilterVarscan_ch0
        .join(ProcessVarscanSNP_out_ch0, by :[0])
        .join(ProcessVarscanIndel_out_ch0, by :[0])

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_varscan.snp.Somatic.hc.filtered.vcf"),
        path("${meta.sampleName}_varscan.indel.Somatic.hc.filtered.vcf")
    ) into FilterVarscan_out_ch0

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    """
    cat ${snpSomaticHc} | \\
    awk '{if (!/^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' | \\
    bam-readcount \\
        -q${params.min_map_q} \\
        -b${params.min_base_q} \\
        -w1 \\
        -l /dev/stdin \\
        -f ${RefFasta} \\
        ${bam[0]} | \\
    varscan ${JAVA_Xmx} fpfilter \\
        ${snpSomaticHc} \\
        /dev/stdin \\
        --output-file ${meta.sampleName}_varscan.snp.Somatic.hc.filtered.vcf && \\
    cat ${indelSomaticHc} | \\
    awk '{if (! /^#/) { x = length(\$5) - 1; print \$1,\$2,(\$2+x); }}' | \\
    bam-readcount \\
        -q${params.min_map_q} \\
        -b${params.min_base_q} \\
        -w1 \\
        -l /dev/stdin \\
        -f ${RefFasta} ${bam[0]} | \\
    varscan ${JAVA_Xmx} fpfilter \\
        ${indelSomaticHc} \\
        /dev/stdin \\
        --output-file ${meta.sampleName}_varscan.indel.Somatic.hc.filtered.vcf
    """
}


/*
    1. Merge filtered SNPS and INDELs from VarScan
    2. Rename the sample names (TUMOR/NORMAL) from varscan vcfs to the real samplenames
*/
process 'MergeAndRenameSamplesInVarscanVCF' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/varscan/",
        mode: publishDirMode

    input:
    path(RefDict) from Channel.value(reference.RefDict)

    tuple(
        val(meta),
        path(VarScanSNP_VCF),
        path(VarScanINDEL_VCF)
    ) from FilterVarscan_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        val("varscan"),
        path("${meta.sampleName}_varscan.Somatic.hc.filtered.{vcf.gz,vcf.gz.tbi}")
    ) into (
        MergeAndRenameSamplesInVarscanVCF_out_ch0,
        MergeAndRenameSamplesInVarscanVCF_out_ch1
    )

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${tmpDir}

    bgzip -c ${VarScanSNP_VCF} > ${VarScanSNP_VCF}.gz
    tabix -p vcf ${VarScanSNP_VCF}.gz
    bgzip -c ${VarScanINDEL_VCF} > ${VarScanINDEL_VCF}.gz
    tabix -p vcf ${VarScanINDEL_VCF}.gz

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${VarScanSNP_VCF}.gz \\
        -I ${VarScanINDEL_VCF}.gz \\
        -O ${meta.sampleName}_varscan_combined.vcf.gz \\
        --SEQUENCE_DICTIONARY ${RefDict}

    gatk --java-options ${java_opts} SortVcf \\
        --TMP_DIR ${tmpDir} \\
        -I ${meta.sampleName}_varscan_combined.vcf.gz \\
        -O ${meta.sampleName}_varscan_combined_sorted.vcf.gz \\
        --SEQUENCE_DICTIONARY ${RefDict}

    # rename samples in varscan vcf
    printf "TUMOR ${meta.sampleName}_tumor\nNORMAL ${meta.sampleName}_normal\n" > vcf_rename_${meta.sampleName}_tmp

    bcftools reheader \\
        -s vcf_rename_${meta.sampleName}_tmp \\
        ${meta.sampleName}_varscan_combined_sorted.vcf.gz \\
        > ${meta.sampleName}_varscan.Somatic.hc.filtered.vcf.gz

    tabix -p vcf ${meta.sampleName}_varscan.Somatic.hc.filtered.vcf.gz
    rm -f vcf_rename_${meta.sampleName}_tmp

    """

}

if(have_Mutect1) {
    // Mutect1: calls SNPS from tumor and matched normal sample
    process 'Mutect1scattered' {

        tag "$meta.sampleName"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/mutect1/processing/",
            mode: publishDirMode,
            enabled: params.fullOutput

        input:
        tuple(
            val(meta),
            path(Tumorbam),
            path(Normalbam),
            path(intervals)
        ) from GatherRealignedBamFiles_out_Mutect1scattered_ch0
            .combine(
                SplitIntervals_out_ch6.flatten()
            )

        tuple(
            path(RefFasta),
            path(RefIdx),
            path(RefDict),
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx,
            reference.RefDict ]
        )


        tuple(
            path(DBSNP),
            path(DBSNPIdx)
        ) from Channel.value(
            [ database.DBSNP,
              database.DBSNPIdx ]
        )

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}_${intervals}.raw.{vcf.gz,vcf.gz.idx}"),
            path("${meta.sampleName}_${intervals}.raw.stats.txt")
        ) into Mutect1scattered_out_ch0

        script:
        def cosmic = ( file(params.databases.Cosmic).exists() &&
                       file(params.databases.CosmicIdx).exists()
                     ) ? "--cosmic " + file(params.databases.Cosmic)
                       : ""
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
        """
        mkdir -p ${tmpDir}

        $JAVA7 ${JAVA_Xmx} -Djava.io.tmpdir=${tmpDir} -jar $MUTECT1 \\
            --analysis_type MuTect \\
            --reference_sequence ${RefFasta} \\
            ${cosmic} \\
            --dbsnp ${DBSNP} \\
            -L ${intervals} \\
            --input_file:normal ${Normalbam[0]} \\
            --input_file:tumor ${Tumorbam[0]} \\
            --out ${meta.sampleName}_${intervals}.raw.stats.txt \\
            --vcf ${meta.sampleName}_${intervals}.raw.vcf.gz
        """
    }

    // Merge scattered Mutect1 vcfs
    process 'gatherMutect1VCFs' {

        label 'nextNEOpiENV'

        tag "$meta.sampleName"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/mutect1/",
            saveAs: {
                fileName ->
                    if(fileName.indexOf("_mutect1_raw") >= 0) {
                        targetFile = "raw/" + fileName
                    } else {
                        targetFile = fileName
                    }
                    return targetFile
            },
            mode: publishDirMode

        input:
        tuple(
            path(RefFasta),
            path(RefIdx),
            path(RefDict)
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx,
            reference.RefDict ]
        )

        tuple(
            val(meta),
            path(vcf),
            path(idx),
            path(stats)
        ) from Mutect1scattered_out_ch0
            .toSortedList({a, b -> a[1][0].baseName <=> b[1][0].baseName})
            .flatten()
            .collate(4)
            .groupTuple(by: [0])

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        path("${meta.sampleName}_mutect1_raw.{vcf.gz,vcg.gz.tbi}")
        path("${meta.sampleName}_mutect1_raw.stats.txt")

        tuple(
            val(meta),
            val("mutect1"),
            path("${meta.sampleName}_mutect1_final.{vcf.gz,vcf.gz.tbi}")
        ) into (
            Mutect1_out_ch0,
            Mutect1_out_ch1
        )


        script:
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
        def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
        """
        mkdir -p ${tmpDir}

        gatk --java-options ${java_opts} MergeVcfs \\
            --TMP_DIR ${tmpDir} \\
            -I ${vcf[0].join(" -I ")} \\
            -O ${meta.sampleName}_mutect1_raw.vcf.gz

        gatk SelectVariants \\
            --tmp-dir ${tmpDir} \\
            --variant ${meta.sampleName}_mutect1_raw.vcf.gz \\
            -R ${RefFasta} \\
            --exclude-filtered true \\
            --select 'vc.getGenotype(\"${meta.sampleName}_tumor\").getAD().1 >= ${params.minAD}' \\
            --output ${meta.sampleName}_mutect1_final.vcf.gz


        head -2 ${stats[0]} > ${meta.sampleName}_mutect1_raw.stats.txt
        tail -q -n +3 ${stats.join(" ")} >> ${meta.sampleName}_mutect1_raw.stats.txt
        """
    }
} else {

    log.info "INFO: Mutect1 not available, skipping...."

    Channel.empty().set { Mutect1_out_ch0 }

    GatherRealignedBamFiles_out_Mutect1scattered_ch0
        .map {  item -> tuple(item[0], "mutect1", []) } // pass meta and emptpy [] would be [vcf, idx]
        .set { Mutect1_out_ch1 }

} // END if have MUTECT1


// Strelka2 and Manta
// Will only run with PE libraries

MantaSomaticIndels_out_ch0 = Channel.create()
MantaSomaticIndels_out_NeoFuse_in_ch0 = Channel.create()

process 'MantaSomaticIndels' {

    label 'Manta'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/manta/",
        saveAs: {
            filename ->
                if((filename.indexOf("_diploidSV.vcf") > 0 ||
                    filename.indexOf("_svCandidateGenerationStats.tsv") > 0 ||
                    filename.indexOf("_candidateSV.vcf") > 0) && params.fullOutput) {
                    return "allSV/$filename"
                } else if((filename.indexOf("_diploidSV.vcf") > 0 ||
                            filename.indexOf("_svCandidateGenerationStats.tsv") > 0 ||
                            filename.indexOf("_candidateSV.vcf") > 0) && ! params.fullOutput) {
                    return null
                } else {
                    return "$filename"
                }
        },
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(Tumorbam),
        path(Normalbam),
        path(RegionsBedGz),
        path(RegionsBedGzTbi)
    ) from BaseRecalGATK4_out_MantaSomaticIndels_ch0
        .combine(RegionsBedToTabix_out_ch1)

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
        reference.RefIdx,
        reference.RefDict ]
    )

    output:
    tuple(
        val(meta),
        path(indel_vcf)
    ) into MantaSomaticIndels_out_ch0

    tuple(
        val(meta),
        path(sv_vcf)
    ) into MantaSomaticIndels_out_NeoFuse_in_ch0

    path("${meta.sampleName}_*.{vcf.gz,vcf.gz.tbi}") optional true

    script:
    indel_vcf = (meta.libType == "PE") ? meta.sampleName + "_candidateSmallIndels.{vcf.gz,vcf.gz.tbi}" : []
    sv_vcf = (meta.libType == "PE") ? meta.sampleName + "_somaticSV.{vcf.gz,vcf.gz.tbi}" : []

    def exome_options = params.WES ? "--callRegions ${RegionsBedGz} --exome" : ""

    // only run with PE samples, Manta is not working with SE samples
    if(meta.libType == "PE")
        """
        configManta.py --tumorBam ${Tumorbam[0]} --normalBam  ${Normalbam[0]} \\
            --referenceFasta ${RefFasta} \\
            --runDir manta_${meta.sampleName} ${exome_options}
        manta_${meta.sampleName}/runWorkflow.py -m local -j ${task.cpus}
        cp manta_${meta.sampleName}/results/variants/diploidSV.vcf.gz ${meta.sampleName}_diploidSV.vcf.gz
        cp manta_${meta.sampleName}/results/variants/diploidSV.vcf.gz.tbi ${meta.sampleName}_diploidSV.vcf.gz.tbi
        cp manta_${meta.sampleName}/results/variants/candidateSV.vcf.gz ${meta.sampleName}_candidateSV.vcf.gz
        cp manta_${meta.sampleName}/results/variants/candidateSV.vcf.gz.tbi ${meta.sampleName}_candidateSV.vcf.gz.tbi
        cp manta_${meta.sampleName}/results/variants/candidateSmallIndels.vcf.gz ${meta.sampleName}_candidateSmallIndels.vcf.gz
        cp manta_${meta.sampleName}/results/variants/candidateSmallIndels.vcf.gz.tbi ${meta.sampleName}_candidateSmallIndels.vcf.gz.tbi
        cp manta_${meta.sampleName}/results/variants/somaticSV.vcf.gz ${meta.sampleName}_somaticSV.vcf.gz
        cp manta_${meta.sampleName}/results/variants/somaticSV.vcf.gz.tbi ${meta.sampleName}_somaticSV.vcf.gz.tbi
        cp manta_${meta.sampleName}/results/stats/svCandidateGenerationStats.tsv ${meta.sampleName}_svCandidateGenerationStats.tsv
        """
    else
        """
        """
}

// combine bam/bai with manta indel vcf/idx
// return null if not manta indels in case of single end library
StrelkaSomatic_in_ch = BaseRecalGATK4_out_StrelkaSomatic_ch0
                        .join(MantaSomaticIndels_out_ch0, by:0, remainder: true)
                        .map {
                            it ->
                            meta       = it[0]
                            tumor_bam  = it[1]
                            normal_bam = it[2]

                            if (it[3] != null) {
                                manta_indel_vfc = it[3]
                                return [meta, tumor_bam, normal_bam, manta_indel_vfc]
                            } else {
                                return [meta, tumor_bam, normal_bam, []]
                            }
                        }

process StrelkaSomatic {

    label 'Strelka'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/strelka/",
        saveAs: { filename -> filename.indexOf("_runStats") > 0 ? "stats/$filename" : "raw/$filename"},
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(tumor_bam),
        path(normal_bam),
        path(manta_indel),
        path(RegionsBedGz),
        path(RegionsBedGzTbi)
    ) from StrelkaSomatic_in_ch
        .combine(RegionsBedToTabix_out_ch0)

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_somatic.snvs.{vcf.gz,vcf.gz.tbi}"),
        path("${meta.sampleName}_somatic.indels.{vcf.gz,vcf.tbi}")
    ) into (
        StrelkaSomatic_out_ch0
    )
    path("${meta.sampleName}_runStats.tsv")
    path("${meta.sampleName}_runStats.xml")

    script:
    def manta_indel_candidates = (manta_indel[0] == null) ? "" : "--indelCandidates " + manta_indel[0]
    def exome_options = params.WES ? "--callRegions ${RegionsBedGz} --exome" : ""

    """
    configureStrelkaSomaticWorkflow.py --tumorBam ${tumor_bam[0]} --normalBam  ${normal_bam[0]} \\
        --referenceFasta ${RefFasta} \\
        ${manta_indel_candidates} \\
        --runDir strelka_${meta.sampleName} ${exome_options}
    strelka_${meta.sampleName}/runWorkflow.py -m local -j ${task.cpus}
    cp strelka_${meta.sampleName}/results/variants/somatic.indels.vcf.gz ${meta.sampleName}_somatic.indels.vcf.gz
    cp strelka_${meta.sampleName}/results/variants/somatic.indels.vcf.gz.tbi ${meta.sampleName}_somatic.indels.vcf.gz.tbi
    cp strelka_${meta.sampleName}/results/variants/somatic.snvs.vcf.gz ${meta.sampleName}_somatic.snvs.vcf.gz
    cp strelka_${meta.sampleName}/results/variants/somatic.snvs.vcf.gz.tbi ${meta.sampleName}_somatic.snvs.vcf.gz.tbi
    cp strelka_${meta.sampleName}/results/stats/runStats.tsv ${meta.sampleName}_runStats.tsv
    cp strelka_${meta.sampleName}/results/stats/runStats.xml ${meta.sampleName}_runStats.xml

    """
}

process 'finalizeStrelkaVCF' {
    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/strelka/",
        saveAs: { filename -> filename.indexOf("_strelka_combined_somatic.vcf.gz") > 0 ? "raw/$filename" : "$filename"},
        mode: publishDirMode


    input:
    tuple(
        val(meta),
        path(somatic_snvs),
        path(somatic_indels),
    ) from StrelkaSomatic_out_ch0

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        val("strelka"),
        path("${meta.sampleName}_strelka_somatic_final.{vcf.gz,vcf.gz.tbi}")
    ) into (
        StrelkaSomaticFinal_out_ch0,
        StrelkaSomaticFinal_out_ch1
    )
    path("${meta.sampleName}_strelka_combined_somatic.vcf.gz")

    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${somatic_snvs[0]} \\
        -I ${somatic_indels[0]} \\
        -O ${meta.sampleName}_strelka_combined.vcf.gz \\
        --SEQUENCE_DICTIONARY ${RefDict}

    gatk --java-options ${java_opts} SortVcf \\
        --TMP_DIR ${tmpDir} \\
        -I ${meta.sampleName}_strelka_combined.vcf.gz \\
        -O ${meta.sampleName}_strelka_combined_sorted.vcf.gz \\
        --SEQUENCE_DICTIONARY ${RefDict}

    # rename samples in strelka vcf
    printf "TUMOR ${meta.sampleName}_tumor\nNORMAL ${meta.sampleName}_normal\n" > vcf_rename_${meta.sampleName}_tmp

    bcftools reheader \\
        -s vcf_rename_${meta.sampleName}_tmp \\
        ${meta.sampleName}_strelka_combined_sorted.vcf.gz \\
        > ${meta.sampleName}_strelka_combined_somatic.vcf.gz

    tabix -p vcf ${meta.sampleName}_strelka_combined_somatic.vcf.gz
    rm -f vcf_rename_${meta.sampleName}_tmp

    gatk SelectVariants \\
        --tmp-dir ${tmpDir} \\
        --variant ${meta.sampleName}_strelka_combined_somatic.vcf.gz \\
        -R ${RefFasta} \\
        --exclude-filtered true \\
        --output ${meta.sampleName}_strelka_somatic_final.vcf.gz

    """
}

// END Strelka2 and Manta

/*
    Creates a VCF that is based on the primary caller (e.g. mutect2) vcf but that contains only variants
    that are confirmed by any of the confirming callers (e..g. mutect1, varscan)
*/
process 'mkHCsomaticVCF' {

    label 'nextNEOpiENV'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/high_confidence/",
        mode: publishDirMode

    input:
    path(RefDict) from Channel.value(reference.RefDict)

    tuple(
        val(meta),
        _,
        path(Mutect2_VCF),
        _,
        path(Mutect1_VCF),
        _,
        path(VarScan_VCF),
        _,
        path(Strelka_VCF)
    ) from FilterMutect2_out_ch1
        .join(Mutect1_out_ch1, by: [0])
        .join(MergeAndRenameSamplesInVarscanVCF_out_ch1, by: [0])
        .join(StrelkaSomaticFinal_out_ch0, by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        val("hc"),
        path("${meta.sampleName}_Somatic.hc.{vcf.gz,vcf.gz.tbi}")
    ) into (
        mkHCsomaticVCF_out_ch0,
        mkHCsomaticVCF_out_ch1,
        mkHCsomaticVCF_out_ch2
    )
    tuple(
        val(meta),
        env(VARCOUNT)
    ) into hc_var_count

    script:
    def callerMap = [:]
    callerMap.M2 = (Mutect2_VCF.size() > 0) ? Mutect2_VCF[0] : ""
    callerMap.M1 = (Mutect1_VCF.size() > 0) ? Mutect1_VCF[0] : ""
    callerMap.VS = (VarScan_VCF.size() > 0) ? VarScan_VCF[0] : ""
    callerMap.ST = (Strelka_VCF.size() > 0) ? Strelka_VCF[0] : ""

    if(!have_Mutect1) {
        callerMap.remove("M1")
    }

    def primary_caller_file = callerMap[params.primaryCaller]

    callerMap.remove(params.primaryCaller)
    def confirming_caller_files = callerMap.values().join(" ")
    def confirming_caller_names = callerMap.keySet().join(" ")
    """
    make_hc_vcf.py \\
        --primary ${primary_caller_file} \\
        --primary_name ${params.primaryCaller} \\
        --confirming ${confirming_caller_files} \\
        --confirming_names ${confirming_caller_names} \\
        --out_vcf ${meta.sampleName}_Somatic.hc.vcf \\
        --out_single_vcf ${meta.sampleName}_Somatic.single.vcf \\
    > var_count.txt

    VARCOUNT=\$(cut -f2 var_count.txt)

    bgzip -c ${meta.sampleName}_Somatic.hc.vcf > ${meta.sampleName}_Somatic.hc.vcf.gz
    tabix -p vcf ${meta.sampleName}_Somatic.hc.vcf.gz
    """

}

vep_cache_chck_file_name = "." + params.vep_species + "_" + params.vep_assembly + "_" + params.vep_cache_version + "_cache_ok.chck"
vep_cache_chck_file = file(params.databases.vep_cache + "/" + vep_cache_chck_file_name)
if(!vep_cache_chck_file.exists() || vep_cache_chck_file.isEmpty()) {

    log.warn "WARNING: VEP cache not installed, starting installation. This may take a while."

    process 'installVEPcache' {

        label 'VEP'

        tag 'installVEPcache'

        // do not cache
        cache false

        output:
        path("${vep_cache_chck_file_name}") into (
            vep_cache_ch0,
            vep_cache_ch1,
            vep_cache_ch2
        )

        script:
        if(!have_vep)
            """
            mkdir -p ${params.databases.vep_cache}
            vep_install \\
                -a cf \\
                -s ${params.vep_species} \\
                -y ${params.vep_assembly} \\
                -c ${params.databases.vep_cache} \\
                --CACHE_VERSION ${params.vep_cache_version} \\
                --CONVERT 2> vep_errors.txt && \\
            echo "OK" > ${vep_cache_chck_file_name} && \\
            cp -f  ${vep_cache_chck_file_name} ${vep_cache_chck_file}
            """
        else
            """
            echo "OK" > ${vep_cache_chck_file_name} && \\
            cp -f  ${vep_cache_chck_file_name} ${vep_cache_chck_file}
            """
    }

} else {

    vep_cache_ch = Channel.fromPath(vep_cache_chck_file)
    (vep_cache_ch0, vep_cache_ch1, vep_cache_ch2) = vep_cache_ch.into(3)

}

vep_plugins_chck_file_name = "." + params.vep_cache_version + "_plugins_ok.chck"
vep_plugins_chck_file = file(params.databases.vep_cache + "/" + vep_plugins_chck_file_name)
if(!vep_plugins_chck_file.exists() || vep_plugins_chck_file.isEmpty()) {

    log.warn "WARNING: VEP plugins not installed, starting installation. This may take a while."

    process 'installVEPplugins' {

        label 'VEP'

        tag 'installVEPplugins'

        // do not cache
        cache false

        input:
        path(vep_cache_chck_file) from vep_cache_ch2

        output:
        path("${vep_plugins_chck_file_name}") into (
            vep_plugins_ch0,
            vep_plugins_ch1
        )

        script:
        if(!have_vep)
            """
            mkdir -p ${params.databases.vep_cache}
            vep_install \\
                -a p \\
                -c ${params.databases.vep_cache} \\
                --PLUGINS all 2> vep_errors.txt && \\
            cp -f ${baseDir}/assets/Wildtype.pm ${params.databases.vep_cache}/Plugins && \\
            cp -f ${baseDir}/assets/Frameshift.pm ${params.databases.vep_cache}/Plugins && \\
            echo "OK" > ${vep_plugins_chck_file_name} && \\
            cp -f  ${vep_plugins_chck_file_name} ${vep_plugins_chck_file}
            """
        else
            """
            echo "OK" > ${vep_plugins_chck_file_name} && \\
            cp -f  ${vep_plugins_chck_file_name} ${vep_plugins_chck_file}
            """
    }

} else {

    vep_plugins_ch = Channel.fromPath(vep_plugins_chck_file)
    (vep_plugins_ch0, vep_plugins_ch1) = vep_plugins_ch.into(2)

}

// Variant Effect Prediction: using ensembl vep
process 'VEPtab' {

    label 'VEP'

    tag "$meta.sampleName"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/05_vep/tables/",
        saveAs: {
            filename ->
                (filename.indexOf("$CallerName") > 0 && CallerName != "hc")
                ? "$CallerName/$filename"
                : "high_confidence/$filename"
        },
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        val(CallerName),
        path(Vcf),
        path(tbi),
        path(vep_cache_chck_file),
        path(vep_plugin_chck_file)
    ) from FilterMutect2_out_ch0
        .concat(MergeAndRenameSamplesInVarscanVCF_out_ch0)
        .concat(Mutect1_out_ch0)
        .concat(StrelkaSomaticFinal_out_ch1)
        .concat(mkHCsomaticVCF_out_ch0)
        .flatten()
        .collate(4)
        .combine(vep_cache_ch0)
        .combine(vep_plugins_ch0)

    output:
    path("${meta.sampleName}_${CallerName}_vep.txt")
    path("${meta.sampleName}_${CallerName}_vep_summary.html")

    script:
    """
    vep -i ${Vcf[0]} \\
        -o ${meta.sampleName}_${CallerName}_vep.txt \\
        --fork ${task.cpus} \\
        --stats_file ${meta.sampleName}_${CallerName}_vep_summary.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --dir ${params.databases.vep_cache} \\
        --cache \\
        --cache_version ${params.vep_cache_version} \\
        --dir_cache ${params.databases.vep_cache} \\
        --fasta ${params.references.VepFasta} \\
        --format "vcf" \\
        ${params.vep_options} \\
        --tab 2> vep_errors.txt
    """
}

// CREATE phased VCF
/*
    make phased vcf for pVACseq using tumor and germline variants:
    based on https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/proximal_vcf.html
*/

// combined germline and somatic variants
process 'mkCombinedVCF' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/high_confidence_readbacked_phased/processing/",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict ]
    )

    tuple(
        val(meta),
        path(germlineVCF),
        _,
        path(tumorVCF)
    ) from FilterGermlineVariantTranches_out_ch0
        .join(mkHCsomaticVCF_out_ch1, by: [0])          // uses confirmed mutect2 variants

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_tumor_germline_combined_sorted.{vcf.gz,vcf.gz.tbi}")
    ) into mkCombinedVCF_out_ch


    script:
    def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
    def java_opts = '"' + JAVA_Xmx + ' -XX:ParallelGCThreads=' + task.cpus + '"'
    """
    mkdir -p ${tmpDir}

    gatk --java-options ${JAVA_Xmx} SelectVariants \\
        --tmp-dir ${tmpDir} \\
        -R ${RefFasta} \\
        -V ${tumorVCF[0]} \\
        --sample-name ${meta.sampleName}_tumor \\
        -O ${meta.sampleName}_tumor.vcf.gz

    gatk --java-options ${java_opts} RenameSampleInVcf \\
        --TMP_DIR ${tmpDir} \\
        -I ${germlineVCF[0]} \\
        --NEW_SAMPLE_NAME ${meta.sampleName}_tumor \\
        -O ${meta.sampleName}_germline_rename2tumorID.vcf.gz

    gatk --java-options ${java_opts} MergeVcfs \\
        --TMP_DIR ${tmpDir} \\
        -I ${meta.sampleName}_tumor.vcf.gz \\
        -I ${meta.sampleName}_germline_rename2tumorID.vcf.gz \\
        -O ${meta.sampleName}_tumor_germline_combined.vcf.gz

    gatk --java-options ${java_opts} SortVcf \\
        --TMP_DIR ${tmpDir} \\
        -I ${meta.sampleName}_tumor_germline_combined.vcf.gz \\
        -O ${meta.sampleName}_tumor_germline_combined_sorted.vcf.gz \\
        --SEQUENCE_DICTIONARY ${RefDict}
    """
}

process 'VEPvcf' {

    label 'VEP'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/05_vep/vcf/high_confidence/",
        saveAs: {
            filename ->
                if (filename.indexOf("_vep_pick.vcf") > 0 && params.fullOutput) {
                    return "combined/$filename"
                } else if (filename.indexOf("_vep_pick.vcf") > 0 && ! params.fullOutput) {
                    return null
                } else if (filename.endsWith(".fa")) {
                    return "../../../06_proteinseq/$filename"
                } else {
                    return "$filename"
                }
        },
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(combinedVCF),
        _,
        path(tumorVCF),
        val(var_count),
        path(vep_cache_chck_file),
        path(vep_plugin_chck_file)
    ) from mkCombinedVCF_out_ch
        .join(mkHCsomaticVCF_out_ch2, by: [0])
        .join(hc_var_count)
        .combine(vep_cache_ch1)
        .combine(vep_plugins_ch1)

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_tumor_germline_combined_sorted_vep_pick.{vcf.gz,vcf.gz.tbi}")
    ) into (
        VEPvcf_out_ch0
    )

    tuple(
        val(meta),
        path("${meta.sampleName}_hc_vep_pick.{vcf.gz,vcf.gz.tbi}")
    ) into (
        VEPvcf_out_ch2,
    )

    tuple(
        val(meta),
        path("${meta.sampleName}_hc_vep.{vcf.gz,vcf.gz.tbi}")
    ) into (
        VEPvcf_out_ch1, // mkPhasedVCF_out_Clonality_ch0
        VEPvcf_out_ch3,
        VEPvcf_out_ch4
    )
    path("${meta.sampleName}_hc_reference.fa")
    path("${meta.sampleName}_hc_mutated.fa")


    when:
    var_count.toInteger() != 0

    script:
    """
    mkdir -p ${tmpDir}

    # pVACSeq
    vep -i ${combinedVCF[0]} \\
        -o ${meta.sampleName}_tumor_germline_combined_sorted_vep_pick.vcf \\
        --fork ${task.cpus} \\
        --stats_file ${meta.sampleName}_tumor_germline_combined_sorted_vep_summary_pick.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --cache \\
        --cache_version ${params.vep_cache_version} \\
        --dir ${params.databases.vep_cache} \\
        --dir_cache ${params.databases.vep_cache} \\
        --hgvs \\
        --fasta ${params.references.VepFasta} \\
        --pick --plugin Frameshift --plugin Wildtype \\
        --symbol --terms SO --transcript_version --tsl \\
        --format vcf \\
        --vcf 2> vep_errors_0.txt

    # pVACSeq
    vep -i ${tumorVCF[0]} \\
        -o ${meta.sampleName}_hc_vep_pick.vcf \\
        --fork ${task.cpus} \\
        --stats_file ${meta.sampleName}_hc_vep_summary_pick.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --cache \\
        --cache_version ${params.vep_cache_version} \\
        --dir ${params.databases.vep_cache} \\
        --dir_cache ${params.databases.vep_cache} \\
        --hgvs \\
        --fasta ${params.references.VepFasta} \\
        --pick --plugin Frameshift --plugin Wildtype \\
        --symbol --terms SO --transcript_version --tsl \\
        --vcf 2>> vep_errors_1.txt

    # All variants
    vep -i ${tumorVCF[0]} \\
        -o ${meta.sampleName}_hc_vep.vcf \\
        --fork ${task.cpus} \\
        --stats_file ${meta.sampleName}_hc_vep_summary.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --offline \\
        --cache \\
        --cache_version ${params.vep_cache_version} \\
        --dir ${params.databases.vep_cache} \\
        --dir_cache ${params.databases.vep_cache} \\
        --hgvs \\
        --fasta ${params.references.VepFasta} \\
        --plugin ProteinSeqs,${meta.sampleName}_hc_reference.fa,${meta.sampleName}_hc_mutated.fa \\
        --symbol --terms SO --transcript_version --tsl \\
        --vcf 2>> vep_errors_1.txt


    bgzip -c ${meta.sampleName}_tumor_germline_combined_sorted_vep_pick.vcf \\
        > ${meta.sampleName}_tumor_germline_combined_sorted_vep_pick.vcf.gz

    tabix -p vcf ${meta.sampleName}_tumor_germline_combined_sorted_vep_pick.vcf.gz && \\
        sleep 2

    bgzip -c ${meta.sampleName}_hc_vep_pick.vcf \\
        > ${meta.sampleName}_hc_vep_pick.vcf.gz

    tabix -p vcf ${meta.sampleName}_hc_vep_pick.vcf.gz && \\
        sleep 2

    bgzip -c ${meta.sampleName}_hc_vep.vcf \\
        > ${meta.sampleName}_hc_vep.vcf.gz

    tabix -p vcf ${meta.sampleName}_hc_vep.vcf.gz && \\
        sleep 2

    sync
    """
}

if(have_GATK3) {
    process 'ReadBackedphasing' {

        label 'GATK3'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/04_variations/high_confidence_readbacked_phased/",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(tumorBAM),
            path(combinedVCF)
        ) from GatherRealignedBamFilesTumor_out_mkPhasedVCF_ch0
            .join(VEPvcf_out_ch0, by: [0])

        tuple(
            path(RefFasta),
            path(RefIdx),
            path(RefDict)
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx,
            reference.RefDict ]
        )

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}_vep_phased.{vcf.gz,vcf.gz.tbi}"),
        ) into (
            mkPhasedVCF_out_ch0,
            mkPhasedVCF_out_pVACseq_ch0,
            generate_protein_fasta_phased_vcf_ch0
        )

        script:
        def JAVA_Xmx = '-Xmx' + task.memory.toGiga() + "G"
        """
        $JAVA8 ${JAVA_Xmx} -XX:ParallelGCThreads=${task.cpus} -Djava.io.tmpdir=${tmpDir} -jar $GATK3 \\
            -T ReadBackedPhasing \\
            -R ${RefFasta} \\
            -I ${tumorBAM[0]} \\
            -V ${combinedVCF[0]} \\
            -L ${combinedVCF[0]} \\
            -o ${meta.sampleName}_vep_phased.vcf.gz
        """
    }
} else {

    log.warn "WARNING: GATK3 not installed! Can not generate readbacked phased VCF:\n" +
        "You should manually review the sequence data for all candidates (e.g. in IGV) for proximal variants and\n" +
        " either account for these manually, or eliminate these candidates. Failure to do so may lead to inclusion\n" +
        " of incorrect peptide sequences."

    (mkPhasedVCF_out_ch0, mkPhasedVCF_out_pVACseq_ch0, generate_protein_fasta_phased_vcf_ch0) = VEPvcf_out_ch0.into(3)
}
// END CREATE phased VCF


// CNVs: ASCAT + FREEC


// adopted from sarek nfcore
process AlleleCounter {

    label 'AlleleCounter'

    tag "$meta.sampleName : $meta.sampleType"

    input:
    tuple(
        val(meta),
        path(BAM)
    ) from GatherRealignedBamFiles_out_AlleleCounter_ch0

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
        path(acLoci)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.acLoci ]
    )

    output:
    tuple(
        val(meta),
        path(outFileName)
    ) into AlleleCounter_out_ch0


    script:
    outFileName = (meta.sampleType == "tumor_DNA") ? meta.sampleName + "_tumor.alleleCount" : meta.sampleName + "_normal.alleleCount"
    def single_end = (meta.libType == "SE") ? "-f 0" : ""
    """
    alleleCounter \\
        -l ${acLoci} \\
        -d \\
        -r ${RefFasta} \\
        -b ${BAM[0]} \\
        ${single_end} \\
        -o ${outFileName}
    """

}

AlleleCounter_out_ch0.branch {
        meta_ori, countFile ->
            def meta = meta_ori.clone()
            tumor : meta.sampleType == "tumor_DNA"
                meta.remove('sampleType')
                return [meta, countFile]

            normal: meta.sampleType == "normal_DNA"
                meta.remove('sampleType')
                return [meta, countFile]
}.set{ AlleleCounter_out_ch0 }

AlleleCounter_out_ch = AlleleCounter_out_ch0.tumor.join(AlleleCounter_out_ch0.normal, by: [0])

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process ConvertAlleleCounts {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/ASCAT/processing",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(alleleCountTumor),
        path(alleleCountNormal),
    ) from AlleleCounter_out_ch

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}.BAF"),
        path("${meta.sampleName}.LogR"),
        path("${meta.sampleName}_normal.BAF"),
        path("${meta.sampleName}_normal.LogR")
    ) into ConvertAlleleCounts_out_ch

    script:
    def sex = (meta.sex == "None") ? "XY" : meta.sex
    """
    Rscript ${baseDir}/bin/convertAlleleCounts.r \\
        ${meta.sampleName} ${alleleCountTumor} ${meta.sampleName}_normal ${alleleCountNormal} ${sex}
    """
}

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process 'Ascat' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/ASCAT/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(bafTumor),
        path(logrTumor),
        path(bafNormal),
        path(logrNormal)
    ) from ConvertAlleleCounts_out_ch

    path(acLociGC) from Channel.value(reference.acLociGC)

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        meta,
        path("${meta.sampleName}.cnvs.txt"),
        path("${meta.sampleName}.purityploidy.txt"),
    ) into Ascat_out_Clonality_ch0
    path("${meta.sampleName}.*.{png,txt}")


    script:
    def sex = (meta.sex == "None") ? "XY" : meta.sex
    """
    # get rid of "chr" string if there is any
    for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
    Rscript ${baseDir}/bin/run_ascat.r ${bafTumor} ${logrTumor} ${bafNormal} ${logrNormal} ${meta.sampleName} ${baseDir} ${acLociGC} ${sex}
    """
}

(Ascat_out_Clonality_ch1, Ascat_out_Clonality_ch0) = Ascat_out_Clonality_ch0.into(2)

if (params.controlFREEC) {
    process 'Mpileup4ControFREEC' {

        label 'nextNEOpiENV'

        tag "$meta.sampleName - $meta.sampleType"

        input:
        tuple(
            val(meta),
            path(BAM)
        ) from GatherRealignedBamFiles_out_Mpileup4ControFREEC_ch0

        each path(interval) from ScatteredIntervalListToBed_out_ch1.flatten()


        tuple(
            path(RefFasta),
            path(RefIdx),
            path(RefDict),
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx,
            reference.RefDict ]
        )

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}_${meta.sampleType}_${interval}.pileup.gz")
        ) into Mpileup4ControFREEC_out_ch0

        script:
        """
        samtools mpileup \\
            -q 1 \\
            -f ${RefFasta} \\
            -l ${interval} \\
            ${BAM[0]} | \\
        bgzip --threads ${task.cpus} -c > ${meta.sampleName}_${meta.sampleType}_${interval}.pileup.gz
        """


    }

    Mpileup4ControFREEC_out_ch0 = Mpileup4ControFREEC_out_ch0.groupTuple(by:[0])

    // Merge scattered pileups
    process 'gatherMpileups' {

        label 'nextNEOpiENV'

        tag "${meta.sampleName} : ${meta.sampleType}"

        input:
        tuple(
            val(meta),
            path(mpileup)
        ) from Mpileup4ControFREEC_out_ch0

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path(outFileName)
        ) into gatherMpileups_out_ch0

        script:
        outFileName = (meta.sampleType == "tumor_DNA") ? meta.sampleName + "_tumor.pileup.gz" : meta.sampleName + "_normal.pileup.gz"
        """
        scatters=`ls -1v *.pileup.gz`
        zcat \$scatters | \\
        bgzip --threads ${task.cpus} -c > ${outFileName}
        """
    }

    gatherMpileups_out_ch0.branch {
            meta_ori, pileupFile ->
                def meta = meta_ori.clone()
                tumor : meta.sampleType == "tumor_DNA"
                    meta.remove('sampleType')
                    return [meta, pileupFile]

                normal: meta.sampleType == "normal_DNA"
                    meta.remove('sampleType')
                    return [meta, pileupFile]
    }.set{ gatherMpileups_out_ch0 }

    gatherMpileups_out_ch0 = gatherMpileups_out_ch0.tumor.join(gatherMpileups_out_ch0.normal, by: [0])

    // run ControlFREEC : adopted from nfcore sarek
    process 'ControlFREEC' {

        label 'Freec'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/controlFREEC/",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(mpileupTumor),
            path(mpileupNormal)
        ) from gatherMpileups_out_ch0

        tuple(
            path(RefChrDir),
            path(RefChrLen),
            path(DBSNP),
            path(DBSNPIdx)
        ) from Channel.value(
            [ reference.RefChrDir,
            reference.RefChrLen,
            database.DBSNP,
            database.DBSNPIdx ]
        )

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}.pileup.gz_CNVs"),
            path("${meta.sampleName}.pileup.gz_ratio.txt"),
            path("${meta.sampleName}.pileup.gz_BAF.txt")
        ) into ControlFREEC_out_ch0

        script:
        def config = meta.sampleName + "_tumor_vs_normal.config.txt"
        def sex = (meta.sex == "None") ? "XY" : meta.sex

        def read_orientation = (meta.libType == "SE") ? "0" : "FR"
        def minimalSubclonePresence = (params.WES) ? 30 : 20
        def degree = (params.WES) ? 1 : 4
        def noisyData = (params.WES) ? "TRUE" : "FALSE"
        def window = (params.WES) ? 0 : 50000
        def breakPointType = (params.WES) ? 4 : 2
        def breakPointThreshold = (params.WES) ? "1.2" : "0.8"
        def printNA = (params.WES) ? "FALSE" :  "TRUE"
        def readCountThreshold = (params.WES) ? 50 : 10
        def minimalCoveragePerPosition = (params.WES) ? 5 : 0
        def captureRegions = (params.WES) ? "captureRegions = ${reference.RegionsBed}" : ""
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
        echo "sex = ${sex}" >> ${config}
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
        echo "${captureRegions}" >> ${config}
        freec -conf ${config}
        """
    }


    process 'ControlFREECviz' {

        tag "${meta.sampleName}"

        // makeGraph.R and assess_significance.R seem to be instable
        errorStrategy 'ignore'

        publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/controlFREEC/",
            mode: publishDirMode


        input:
        tuple(
            val(meta),
            path(cnvTumor),
            path(ratioTumor),
            path(bafTumor),
        ) from ControlFREEC_out_ch0

        output:
        tuple(
            path("*.txt"),
            path("*.png"),
            path("*.bed"),
            path("*.circos")
        ) into ControlFREECviz_out_ch0


        script:
        """
        cat ${baseDir}/bin/assess_significance.R | R --slave --args ${cnvTumor} ${ratioTumor}
        cat ${baseDir}/bin/makeGraph.R | R --slave --args 2 ${ratioTumor} ${bafTumor}
        perl ${baseDir}/bin/freec2bed.pl -f ${ratioTumor} > ${meta.sampleName}.bed
        perl ${baseDir}/bin/freec2circos.pl -f ${ratioTumor} > ${meta.sampleName}.circos
        """
    }
}

Channel
    .fromPath(reference.RefChrLen)
    .splitCsv(sep: "\t")
    .map { row -> row[1] }
    .set { chromosomes_ch }

process 'SequenzaUtils' {

    label 'sequenzaUtils'

    tag "${meta.sampleName}"

    input:
    tuple(
        val(meta),
        path(tumorBAM),
        path(normalBAM)
    ) from GatherRealignedBamFiles_out_Sequenza_ch0
    each chromosome from chromosomes_ch

    tuple(
        path(RefFasta),
        path(RefIdx),
        path(RefDict),
        path(SequnzaGC)
    ) from Channel.value(
        [ reference.RefFasta,
          reference.RefIdx,
          reference.RefDict,
          reference.SequenzaGC ]
    )

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${chromosome}_${meta.sampleName}_seqz.gz")
    )  into SequenzaUtils_out_ch0

    script:
    """
    sequenza-utils \\
        bam2seqz \\
        --fasta ${RefFasta} \\
        --tumor ${tumorBAM[0]} \\
        --normal ${normalBAM[0]} \\
        -gc ${SequnzaGC} \\
        --chromosome ${chromosome} \\
        | \\
    sequenza-utils \\
        seqz_binning \\
        -w 50 \\
        -s - \\
        | \\
    bgzip \\
        --threads ${task.cpus} -c > ${chromosome}_${meta.sampleName}_seqz.gz

    """
}

process gatherSequenzaInput {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/Sequenza/processing",
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(chromosome_seqz_binned)
    ) from SequenzaUtils_out_ch0
        .groupTuple(by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_seqz.gz")
    ) into gatherSequenzaInput_out_ch0

    script:
    """
    first=1
    scatters=`ls -1v *_${meta.sampleName}_seqz.gz`
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
    bgzip --threads ${task.cpus} -c > ${meta.sampleName}_seqz.gz
    sync
    """
}

process Sequenza {

    label 'sequenzaR'

    tag "${meta.sampleName}"

    errorStrategy "ignore"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/Sequenza/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(seqz_file),
    ) from gatherSequenzaInput_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_segments.txt"),
        path("${meta.sampleName}_confints_CP.txt")
    ) into Sequenza_out_Clonality_ch0
    path("${meta.sampleName}_*.{png,pdf,txt}")

    script:
    def sex = (meta.sex == "None") ? "XY" : meta.sex
    """
    Rscript \\
        ${baseDir}/bin/SequenzaScript.R \\
        ${seqz_file} \\
        ${meta.sampleName} \\
        ${sex} || \\
        touch ${meta.sampleName}_segments.txt && \\
        touch ${meta.sampleName}_confints_CP.txt
    """
}

(Sequenza_out_Clonality_ch1, Sequenza_out_Clonality_ch0) = Sequenza_out_Clonality_ch0.into(2)

// get purity for CNVKit
purity_estimate = Ascat_out_Clonality_ch0.join(Sequenza_out_Clonality_ch0, by: [0])
    .map {

        it ->
        def meta = it[0]
        def ascat_CNVs = it[1]
        def ascat_purity  = it[2]
        def seqz_CNVs  = it[3]
        def seqz_purity  = it[4]

        def ascatOK = true
        def sequenzaOK = true

        def purity = 0.95 // default
        def ploidy = 2.0 // default
        def sample_purity = 0.95// default

        def fileReader = ascat_purity.newReader()

        def line = fileReader.readLine()
        line = fileReader.readLine()
        fileReader.close()
        if(line) {
            (purity, ploidy) = line.split("\t")
            if(purity == "0" || ploidy == "0" ) {
                ascatOK = false
            }
        } else {
            ascatOK = false
        }


        if(ascatOK && ! params.use_sequenza_cnvs) {
            sample_purity = purity
        } else {
            fileReader = seqz_purity.newReader()

            line = fileReader.readLine()
            def fields
            if(line) {
                fields = line.split("\t")
                if(fields.size() < 3) {
                    sequenzaOK = false
                } else {
                    line = fileReader.readLine()
                    line = fileReader.readLine()
                    (purity, ploidy, _) = line.split("\t")
                }
            } else {
                sequenzaOK = false
            }
            fileReader.close()

            if(sequenzaOK) {
                sample_purity = purity
                log.warn "WARNING (" + meta.sampleName + "): changed from ASCAT to Sequenza purity and segments, ASCAT did not produce results"
            } else {
                log.warn "WARNING (" + meta.sampleName + "): neither ASCAT nor Sequenza produced results, using purity of 0.95"
            }
        }

        return [ meta, sample_purity ]
    }

(purity_estimate_ch0, purity_estimate) = purity_estimate.into(2)

purity_estimate.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}.into{purity_estimate_ch1; purity_estimate_ch2}

// CNVkit
if (params.CNVkit) {
    process make_CNVkit_access_file {

        label 'CNVkit'

        tag 'mkCNVkitaccess'

        publishDir "$params.outputDir/supplemental/01_prepare_CNVkit/",
            mode: publishDirMode

        input:
        tuple(
            path(RefFasta),
            path(RefIdx),
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx ]
        )

        output:
        path(
            "access-5kb.${RefFasta.simpleName}.bed"
        ) into make_CNVkit_access_file_out_ch0

        script:
        """
        cnvkit.py \\
            access \\
            ${RefFasta} \\
            -s 5000 \\
            -o access-5kb.${RefFasta.simpleName}.bed
        """
    }

    process CNVkit {

        label 'CNVkit'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/08_CNVs/CNVkit/",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(tumorBAM),
            path(normalBAM),
            val(sample_purity)
        ) from MarkDuplicates_out_CNVkit_ch0
            .join(purity_estimate_ch0, by: [0])

        path(CNVkit_accessFile) from make_CNVkit_access_file_out_ch0

        tuple(
            path(RefFasta),
            path(RefIdx),
            path(AnnoFile),
        ) from Channel.value(
            [ reference.RefFasta,
            reference.RefIdx,
            reference.AnnoFile ]
        )

        path(BaitsBed) from Channel.fromPath(reference.BaitsBed)

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}*")
        ) into CNVkit_out_ch0

        script:
        def gender = (meta.sex == "None") ? "" : "--sample-sex " + ((meta.sex == "XX") ? "female" : "male")
        def maleRef = (meta.maleRef == "true") ? "-y" : ""
        def method = (params.WES) ? "--method hybrid" : "--method wgs"
        def targets = (params.WES) ? "--targets ${BaitsBed}" : ""

        """
        # set Agg as backend for matplotlib
        export MATPLOTLIBRC="./matplotlibrc"
        echo "backend : Agg" > \$MATPLOTLIBRC

        cnvkit.py \\
            batch \\
            ${tumorBAM[0]} \\
            --normal ${normalBAM[0]} \\
            ${method} \\
            ${targets} \\
            --fasta ${RefFasta} \\
            --annotate ${AnnoFile} \\
            --access ${CNVkit_accessFile} \\
            ${maleRef} \\
            -p ${task.cpus} \\
            --output-reference output_reference.cnn \\
            --output-dir ./

        cnvkit.py segmetrics \\
            -s ${tumorBAM[0].baseName}.cn{s,r} \\
            --ci \\
            --pi

        cnvkit.py call \\
            ${tumorBAM[0].baseName}.cns \\
            --filter ci \\
            -m clonal \\
            --purity ${sample_purity} \\
            ${gender} \\
            ${maleRef} \\
            -o ${tumorBAM[0].baseName}.call.cns

        cnvkit.py \\
            scatter \\
            ${tumorBAM[0].baseName}.cnr \\
            -s ${tumorBAM[0].baseName}.cns \\
            -o ${tumorBAM[0].baseName}_scatter.png

        cnvkit.py \\
            diagram \\
            ${tumorBAM[0].baseName}.cnr \\
            -s ${tumorBAM[0].baseName}.cns \\
            ${gender} \\
            ${maleRef} \\
            -o ${tumorBAM[0].baseName}_diagram.pdf

        cnvkit.py \\
            breaks \\
            ${tumorBAM[0].baseName}.cnr ${tumorBAM[0].baseName}.cns \\
            -o ${tumorBAM[0].baseName}_breaks.tsv

        cnvkit.py \\
            genemetrics \\
            ${tumorBAM[0].baseName}.cnr \\
            -s ${tumorBAM[0].baseName}.cns \\
            ${gender} \\
            ${maleRef} \\
            -t 0.2 -m 5 \\
            -o ${tumorBAM[0].baseName}_gainloss.tsv

        # run PDF to PNG conversion if mogrify and gs is installed
        mogrify -version > /dev/null 2>&1 && \\
        gs -v > /dev/null 2>&1 && \\
            mogrify -density 600 -resize 2000 -format png *.pdf

        # clean up
        rm -f \$MATPLOTLIBRC
        """
    }
}

clonality_input = Ascat_out_Clonality_ch1.join(Sequenza_out_Clonality_ch1, by: [0])
    .map {

        it ->
        def meta = it[0]
        def ascat_CNVs = it[1]
        def ascat_purity  = it[2]
        def seqz_CNVs  = it[3]
        def seqz_purity  = it[4]

        def ascatOK = true
        def sequenzaOK = true

        def fileReader = ascat_purity.newReader()

        def line = fileReader.readLine()
        line = fileReader.readLine()
        fileReader.close()
        if(line) {
            def (purity, ploidy) = line.split("\t")
            if(purity == "0" || ploidy == "0" ) {
                ascatOK = false
            }
        } else {
            ascatOK = false
        }

        fileReader = ascat_CNVs.newReader()

        def fields
        line = fileReader.readLine()
        fileReader.close()
        if(line) {
            fields = line.split("\t")
            if(fields.size() < 5) {
                ascatOK = false
            }
        } else {
            ascatOK = false
        }


        fileReader = seqz_CNVs.newReader()

        line = fileReader.readLine()
        fileReader.close()
        if(line) {
            fields = line.split("\t")
            if(fields.size() < 13) {
                sequenzaOK = false
            }
        } else {
            sequenzaOK = false
        }

        fileReader = seqz_purity.newReader()

        line = fileReader.readLine()
        fileReader.close()
        if(line) {
            fields = line.split("\t")
            if(fields.size() < 3) {
                sequenzaOK = false
            }
        } else {
            sequenzaOK = false
        }

        return [ meta, file(ascat_CNVs), file(ascat_purity), file(seqz_CNVs), file(seqz_purity), ascatOK, sequenzaOK ]
    }


process 'Clonality' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/09_CCF/",
        mode: publishDirMode

    cache 'lenient'

    input:
    tuple(
        val(meta),
        path(hc_vep_vcf),
        path(ascat_CNVs),
        path(ascat_purity),
        path(seqz_CNVs),
        path(seqz_purity),
        val(ascatOK),
        val(sequenzaOK)
    ) from VEPvcf_out_ch1 // mkPhasedVCF_out_Clonality_ch0
        .join(clonality_input, by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_CCFest.tsv"),
        val(ascatOK),
        val(sequenzaOK)
    ) into (
        Clonality_out_ch0,
        ccf_ch0,
        ccf_ch1
    )

    script:
    def seg_opt = ""
    def purity_opt = ""
    if (ascatOK && ! params.use_sequenza_cnvs) {
        seg_opt = "--seg ${ascat_CNVs}"
        purity_opt = "--purity ${ascat_purity}"
    } else if (sequenzaOK) {
        if(! params.use_sequenza_cnvs) {
            log.warn "WARNING: changed from ASCAT to Sequenza purity and segments, ASCAT did not produce results"
        }
        seg_opt = "--seg_sequenza ${seqz_CNVs}"
        purity_opt = "--purity_sequenza ${seqz_purity}"
    } else {
        log.warn "WARNING: neither ASCAT nor Sequenza did produce results"
    }

    if (ascatOK || sequenzaOK)
        """
        mkCCF_input.py \\
            --PatientID ${meta.sampleName}_tumor \\
            --vcf ${hc_vep_vcf[0]} \\
            ${seg_opt} \\
            ${purity_opt} \\
            --min_vaf 0.01 \\
            --result_table ${meta.sampleName}_segments_CCF_input.txt && \\
        Rscript \\
            ${baseDir}/bin/CCF.R \\
            ${meta.sampleName}_segments_CCF_input.txt \\
            ${meta.sampleName}_CCFest.tsv \\
        """
    else
        """
        echo "Not avaliable" > ${meta.sampleName}_CCFest.tsv
        """
}

// mutational burden all variants all covered positions
process 'MutationalBurden' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/07_MutationalBurden/",
        mode: publishDirMode

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    cache 'lenient'

    input:
    tuple(
        val(meta),
        path(Tumorbam),
        path(Normalbam),
        path(vep_somatic_vcf_gz),
        path(ccf_file),
        val(ascatOK),
        val(sequenzaOK)
    ) from BaseRecalGATK4_out_MutationalBurden_ch0
        .join(VEPvcf_out_ch3, by: [0])
        .join(ccf_ch0, by: [0])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_mutational_burden.txt")
    ) into sample_info_tmb


    script:
    def ccf_opts = ""

    if (ascatOK || sequenzaOK) {
        ccf_opts =  "--ccf ${ccf_file} --ccf_clonal_thresh ${params.CCFthreshold} --p_clonal_thresh ${params.pClonal}"
    }
    """
    mutationalLoad.py \\
        --normal_bam ${Normalbam[0]} \\
        --tumor_bam ${Tumorbam[0]} \\
        --vcf ${vep_somatic_vcf_gz[0]} \\
        --min_coverage 5 \\
        --min_BQ 20 \\
        ${ccf_opts} \\
        --cpus ${task.cpus} \\
        --output_file ${meta.sampleName}_mutational_burden.txt
    """
}

// mutational burden coding variants coding (exons) covered positions
process 'MutationalBurdenCoding' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/07_MutationalBurden/",
        mode: publishDirMode

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3

    cache 'lenient'

    input:
    tuple(
        val(meta),
        path(Tumorbam),
        path(Normalbam),
        path(vep_somatic_vcf_gz),
        path(ccf_file),
        val(ascatOK),
        val(sequenzaOK)
    ) from BaseRecalGATK4_out_MutationalBurden_ch1
        .join(VEPvcf_out_ch4, by: [0])
        .join(ccf_ch1, by: [0])
    path (exons) from Channel.value(reference.ExonsBED)

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_mutational_burden_coding.txt")
    ) into sample_info_tmb_coding


    script:
    def ccf_opts = ""

    if (ascatOK || sequenzaOK) {
        ccf_opts =  "--ccf ${ccf_file} --ccf_clonal_thresh ${params.CCFthreshold} --p_clonal_thresh ${params.pClonal}"
    }
    """
    mutationalLoad.py \\
        --normal_bam ${Normalbam[0]} \\
        --tumor_bam ${Tumorbam[0]} \\
        --vcf ${vep_somatic_vcf_gz[0]} \\
        --min_coverage 5 \\
        --min_BQ 20 \\
        --bed ${exons} \\
        --variant_type coding \\
        ${ccf_opts} \\
        --cpus ${task.cpus} \\
        --output_file ${meta.sampleName}_mutational_burden_coding.txt
    """
}


// END CNVs


// HLA TYPING

process 'mhc_extract' {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/mhc_extract",
        mode: publishDirMode,
        saveAs: {
            filename ->
                if(filename.indexOf("NO_FILE") >= 0) {
                    return null
                } else {
                    return "$filename"
                }
        },
        enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(tumor_BAM_aligned_sort_mkdp)
    ) from MarkDuplicatesTumor_out_ch0

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${prefix}_R*.fastq.gz")
    ) into (
        reads_tumor_hla_ch,
        reads_tumor_hlaHD_ch
    )


    script:
    def mhc_region = params.HLA_HD_genome_version ? params.MHC_genomic_region[ params.HLA_HD_genome_version ].region ?: false : false

    if (!mhc_region) {
        exit 1, "MHC region not found for genome version: ${params.HLA_HD_genome_version}"
    }

    prefix = meta.sampleName + "_reads_mhc"

    if(meta.libType == "SE")
        """
        rm -f unmapped_bam mhc_mapped_bam R.fastq
        mkfifo unmapped_bam
        mkfifo mhc_mapped_bam
        mkfifo R.fastq

        samtools  view -@4 -h -b -u -f 4 ${tumor_BAM_aligned_sort_mkdp[0]} > unmapped_bam &
        samtools  view -@4 -h -b -u ${tumor_BAM_aligned_sort_mkdp[0]} ${mhc_region} > mhc_mapped_bam &

        samtools merge -@4 -u - mhc_mapped_bam unmapped_bam | \\
            samtools sort -@4 -n - | \\
            samtools fastq -@2 -0 R.fastq - &
        perl -ple 'if ((\$. % 4) == 1) { s/\$/ 1:N:0:NNNNNNNN/; }' R.fastq | gzip -1 > ${prefix}_R1.fastq.gz

        wait
        touch NO_FILE

        rm -f unmapped_bam mhc_mapped_bam R.fastq
        """
    else
        """
        rm -f unmapped_bam mhc_mapped_bam R1.fastq R2.fastq
        mkfifo unmapped_bam
        mkfifo mhc_mapped_bam
        mkfifo R1.fastq
        mkfifo R2.fastq

        samtools  view -@4 -h -b -u -f 4 ${tumor_BAM_aligned_sort_mkdp[0]} > unmapped_bam &
        samtools  view -@4 -h -b -u ${tumor_BAM_aligned_sort_mkdp[0]} ${mhc_region} > mhc_mapped_bam &

        samtools merge -@4 -u - mhc_mapped_bam unmapped_bam | \\
            samtools sort -@4 -n - | \\
            samtools fastq -@2 -1 R1.fastq -2 R2.fastq -s /dev/null -0 /dev/null - &
        perl -ple 'if ((\$. % 4) == 1) { s/\$/ 1:N:0:NNNNNNNN/; }' R1.fastq | gzip -1 > ${prefix}_R1.fastq.gz &
        perl -ple 'if ((\$. % 4) == 1) { s/\$/ 2:N:0:NNNNNNNN/; }' R2.fastq | gzip -1 > ${prefix}_R2.fastq.gz &

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

if (run_OptiType) {

    process 'pre_map_hla' {

        label 'nextNEOpiENV'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/Optitype/processing/",
            mode: publishDirMode,
            enabled: params.fullOutput

        input:
        tuple(
            val(meta),
            path(reads)
        ) from reads_tumor_hla_ch

        path yaraIdx_files from Channel.value(reference.YaraIndexDNA)
        val yaraIdx from Channel.value(reference.YaraIndexDNA[0].simpleName)

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path("dna_mapped_{1,2}.bam")
        ) into fished_reads

        script:
        def yara_cpus = task.cpus
        def samtools_cpus = 1
        if(meta.libType == "SE") {
            yara_cpus = ((task.cpus - 2).compareTo(2) == -1) ? 2 : (task.cpus - 2)
            samtools_cpus = ((task.cpus - yara_cpus).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus)
        } else {
            yara_cpus = ((task.cpus - 6).compareTo(2) == -1) ? 2 : (task.cpus - 6)
            samtools_cpus =  (((task.cpus - yara_cpus).div(3)).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus).div(3)
        }

        if (meta.libType == "SE")
            """
            yara_mapper -e 3 -t $yara_cpus -f bam ${yaraIdx} ${reads} | \\
                samtools view -@ $samtools_cpus -h -F 4 -b1 -o dna_mapped_1.bam
            """
        else
            """
            rm -f R1 R2
            mkfifo R1 R2
            yara_mapper -e 3 -t $yara_cpus -f bam ${yaraIdx} ${reads} | \\
                samtools view -@ $samtools_cpus -h -F 4 -b1 | \\
                tee R1 R2 > /dev/null &
                samtools view -@ $samtools_cpus -h -f 0x40 -b1 R1 > dna_mapped_1.bam &
                samtools view -@ $samtools_cpus -h -f 0x80 -b1 R2 > dna_mapped_2.bam &
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

        label 'nextNEOpiENV'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/Optitype/",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(reads)
        ) from fished_reads

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}_optitype_result.tsv")
        ) into optitype_output
        path("${meta.sampleName}_optitype_coverage_plot.pdf")

        script:
        """
        OPTITYPE="\$(readlink -f \$(which OptiTypePipeline.py))"
        \$OPTITYPE -i ${reads} -e 1 -b 0.009 --dna -o ./tmp && \\
        mv ./tmp/*/*_result.tsv ./${meta.sampleName}_optitype_result.tsv && \\
        mv ./tmp/*/*_coverage_plot.pdf ./${meta.sampleName}_optitype_coverage_plot.pdf && \\
        rm -rf ./tmp/
        """
    }

    if (! have_RNA_tag_seq) {
        process 'pre_map_hla_RNA' {

            label 'nextNEOpiENV'

            tag "${meta.sampleName}"

            publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/Optitype/processing/",
                mode: publishDirMode,
                enabled: params.fullOutput

            input:
            tuple(
                val(meta),
                path(reads_RNA)
            ) from reads_tumor_optitype_ch

            path yaraIdx_files from Channel.value(reference.YaraIndexRNA)
            val yaraIdx from Channel.value(reference.YaraIndexRNA[0].simpleName)

            val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

            output:
            tuple(
                val(meta),
                path("rna_mapped_{1,2}.bam")
            ) into fished_reads_RNA

            script:
            // check if single end
            def yara_cpus = task.cpus
            def samtools_cpus = 1
            if(meta.libType == "SE") {
                yara_cpus = ((task.cpus - 2).compareTo(2) == -1) ? 2 : (task.cpus - 2)
                samtools_cpus = ((task.cpus - yara_cpus).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus)
            } else {
                yara_cpus = ((task.cpus - 6).compareTo(2) == -1) ? 2 : (task.cpus - 6)
                samtools_cpus =  (((task.cpus - yara_cpus).div(3)).compareTo(1) == -1) ? 1 : (task.cpus - yara_cpus).div(3)
            }

            // check if single end
            if (meta.libType == "SE")
                """
                yara_mapper -e 3 -t $yara_cpus -f bam ${yaraIdx} ${reads_RNA} | \\
                    samtools view -@ $samtools_cpus -h -F 4 -b1 -o rna_mapped_1.bam
                """
            else
                """
                rm -f R1 R2
                mkfifo R1 R2
                yara_mapper -e 3 -t $yara_cpus -f bam ${yaraIdx} ${reads_RNA} | \\
                    samtools view -@ $samtools_cpus -h -F 4 -b1 | \\
                    tee R1 R2 > /dev/null &
                    samtools view -@ $samtools_cpus -h -f 0x40 -b1 R1 > rna_mapped_1.bam &
                    samtools view -@ $samtools_cpus -h -f 0x80 -b1 R2 > rna_mapped_2.bam &
                wait
                rm -f R1 R2
                """
        }

        process 'OptiType_RNA' {

            label 'nextNEOpiENV'

            tag "${meta.sampleName}"

            publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/Optitype/",
                mode: publishDirMode

            input:
            tuple(
                val(meta),
                path(reads)
            ) from fished_reads_RNA

            val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

            output:
            tuple(
                val(meta),
                path("${meta.sampleName}_optitype_RNA_result.tsv")
            ) into optitype_RNA_output
            path("${meta.sampleName}_optitype_RNA_coverage_plot.pdf")

            script:
            def read_count_R1 = "samtools view -c ${reads[0]}"
            def read_count_R2 = (meta.libType == "PE") ? "samtools view -c ${reads[1]}" : ""
            """
            OPTITYPE="\$(readlink -f \$(which OptiTypePipeline.py))"
            MHC_MAPPED_R1=\$($read_count_R1)
            MHC_MAPPED_R2=\$($read_count_R2)
            MHC_MAPPED=\$((MHC_MAPPED_R1+MHC_MAPPED_R2))
            if [ "\$MHC_MAPPED" != "0" ]; then
                \$OPTITYPE -i ${reads} -e 1 -b 0.009 --rna -o ./tmp && \\
                mv ./tmp/*/*_result.tsv ./${meta.sampleName}_optitype_RNA_result.tsv && \\
                mv ./tmp/*/*_coverage_plot.pdf ./${meta.sampleName}_optitype_RNA_coverage_plot.pdf && \\
                rm -rf ./tmp/
            else
                touch ${meta.sampleName}_optitype_RNA_result.tsv
                echo "No result" >  ${meta.sampleName}_optitype_RNA_coverage_plot.pdf
            fi
            """

        }

    } else if (have_RNA_tag_seq) {

        log.info "INFO: will not run HLA typing on RNAseq from tag libraries"

        optitype_RNA_output = reads_tumor_optitype_ch
                                .map{ it -> tuple(it[0], [])}

    }
}  else { // End if run_OptiType
    log.info "INFO: will not run HLA typing with OptiType"

    optitype_output = reads_tumor_hla_ch
                            .map{ it -> tuple(it[0], [])}

    optitype_RNA_output = reads_tumor_optitype_ch
                            .map{ it -> tuple(it[0], [])}
}

/*
*********************************************
**             H L A - H D                 **
*********************************************
*/

if (have_HLAHD) {
    process 'run_hla_hd' {

        label 'HLAHD'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/HLA_HD/",
            saveAs: { fileName -> fileName.endsWith("_final.result.txt") ? file(fileName).getName() : null },
            mode: publishDirMode

        if(params.HLAHD_module != "") {
            module = params.HLAHD_module
        }

        input:
        tuple(
            val(meta),
            path(reads)
        ) from reads_tumor_hlaHD_ch

        path frData from Channel.value(reference.HLAHDFreqData)
        path gSplit from Channel.value(reference.HLAHDGeneSplit)
        path dict from Channel.value(reference.HLAHDDict)

        output:
        tuple(
            val(meta),
            path("**/*_final.result.txt")
        ) into (
            hlahd_output,
            hlahd_mixMHC2_pred_ch0
        )

        script:
        def hlahd_p = Channel.value(HLAHD_PATH).getVal()
        def in_reads = (meta.libType == "PE") ? reads : reads[0] + " " + reads[0]

        """
        export PATH=\$PATH:$hlahd_p
        $HLAHD -t ${task.cpus} \\
            -m 50 \\
            -f ${frData} ${in_reads} \\
            ${gSplit} ${dict} ${meta.sampleName}_${meta.sampleType} .
        """
    }
} else {
    // fill channels
    hlahd_output = reads_tumor_hlaHD_ch
                        .map{ it -> tuple(it[0], [])}

    (hlahd_mixMHC2_pred_ch0, hlahd_output) = hlahd_output.into(2)

}

if (have_HLAHD && ! have_RNA_tag_seq && params.run_HLAHD_RNA) {
    process 'run_hla_hd_RNA' {

        label 'HLAHD_RNA'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/10_HLA_typing/HLA_HD/",
            saveAs: { fileName -> fileName.endsWith("_final.result.txt") ? file(fileName).getName().replace(".txt", ".RNA.txt") : null },
            mode: publishDirMode

        if(params.HLAHD_module != "") {
            module = params.HLAHD_module
        }

        input:
        tuple(
            val(meta),
            path(reads)
        ) from reads_tumor_hlahd_RNA_ch

        path frData from Channel.value(reference.HLAHDFreqData)
        path gSplit from Channel.value(reference.HLAHDGeneSplit)
        path dict from Channel.value(reference.HLAHDDict)

        output:
        tuple(
            val(meta),
            path("**/*_final.result.txt")
        ) into (
            hlahd_output_RNA
        )

        script:
        def hlahd_p = Channel.value(HLAHD_PATH).getVal()
        def in_reads = (meta.libType == "PE") ? reads : reads[0] + " " + reads[0]

        """
        export PATH=\$PATH:$hlahd_p
        $HLAHD -t ${task.cpus} \\
            -m 50 \\
            -f ${frData} ${in_reads} \\
            ${gSplit} ${dict} ${meta.sampleName}_${meta.sampleType} .
        """
    }
} else if ((! have_HLAHD) || (have_RNA_tag_seq) || (! params.run_HLAHD_RNA)) {

    if(have_RNA_tag_seq) {
        log.info "INFO: will not run HLA typing on RNAseq from tag libraries"
    }

    hlahd_output_RNA = reads_tumor_hlahd_RNA_ch
                            .map{ it -> tuple(it[0], [])}

}

/*
Get the HLA types from OptiType and HLA-HD ouput as a "\n" seperated list.
To be used as input for pVACseq

From slack discussion on 20200425 (FF, GF, DR):

    * We run Optitype RNA and WES in the main pipeline (only on tumor)
    * We consider the WES class I HLA
    * IFF the WES-based HLA are homo and are a subset of RNA-based HLA, then we consider also the second HLA alleles predicted from RNA
    * These class I HLA are used for the somatic pipeline and also for the embedded NeoFuse
*/

optitype_output = optitype_output.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

optitype_RNA_output = optitype_RNA_output.mix(no_RNA_0).map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

hlahd_output = hlahd_output.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

hlahd_output_RNA = hlahd_output_RNA.mix(no_RNA_1).map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

batch_custom_HLA_data_ch = batch_custom_HLA_data_ch.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

process get_vhla {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/neoantigens/${meta.sampleName}/Final_HLAcalls/",
    mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(opti_out),
        path(opti_out_rna),
        path(hlahd_out),
        path(hlahd_out_rna),
        path(custom_hlas)
    ) from optitype_output
        .join(optitype_RNA_output, by: 0)
        .join(hlahd_output, by: 0)
        .join(hlahd_output_RNA, by: 0)
        .join(batch_custom_HLA_data_ch, by: 0)

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_hlas.txt")
    ) into (hlas, hlas_neoFuse)

    script:
    def optitype_hlas = (opti_out.size() > 0) ? "--opti_out $opti_out" : ''
    def user_hlas = (custom_hlas.size() > 0) ? "--custom $custom_hlas" : ''
    def rna_hlas = (meta.have_RNA && ! have_RNA_tag_seq && (opti_out_rna.size() > 0)) ? "--opti_out_RNA $opti_out_rna" : ''
    rna_hlas = (meta.have_RNA && have_HLAHD && ! have_RNA_tag_seq && params.run_HLAHD_RNA) ? rna_hlas + " --hlahd_out_RNA $hlahd_out_rna" : rna_hlas
    def force_seq_type = ""

    def force_RNA = (params.HLA_force_DNA || have_RNA_tag_seq) ? false : params.HLA_force_RNA

    if(force_RNA && ! meta.have_RNA) {
        log.warn "WARNING: Can not force RNA data for HLA typing: no RNAseq data provided!"
    } else if (force_RNA && meta.have_RNA) {
        force_seq_type = "--force_RNA"
    } else if (params.HLA_force_DNA) {
        force_seq_type = "--force_DNA"
    }

    hlahd_opt = (have_HLAHD && (hlahd_out.size() > 0)) ? "--hlahd_out ${hlahd_out}" : ""

    def pVACseqAlleles = baseDir.toRealPath()  + "/assets/pVACseqAlleles.txt"
    """
    # merging script
    HLA_parser.py \\
        ${optitype_hlas} \\
        ${hlahd_opt} \\
        ${rna_hlas} \\
        ${user_hlas} \\
        ${force_seq_type} \\
        --ref_hlas ${pVACseqAlleles} \\
        > ./${meta.sampleName}_hlas.txt
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

reads_tumor_neofuse_ch = reads_tumor_neofuse_ch.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, meta_ori, out_file]
}
MantaSomaticIndels_out_NeoFuse_in_ch0 = MantaSomaticIndels_out_NeoFuse_in_ch0.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, meta_ori, out_file]
}

process Neofuse {

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/",
        saveAs: {
            fileName ->
                if(fileName.indexOf("Arriba") >= 0) {
                    targetFile = "11_Fusions/Arriba/" + file(fileName).getName()
                } else if(fileName.indexOf("Custom_HLAs") >= 0) {
                    targetFile = params.fullOutput ? "11_Fusions/Custom_HLAs/" + file(fileName).getName() : ""
                } else if(fileName.indexOf("LOGS") >= 0) {
                    targetFile = params.fullOutput ? "11_Fusions/LOGS/" + file(fileName).getName() : ""
                } else if(fileName.indexOf("_NeoFuse_MHC_Class_I_") >= 0) {
                    targetFile = "11_Fusions/NeoFuse/Class_I/" + file(fileName).getName()
                } else if(fileName.indexOf("_NeoFuse_MHC_Class_II_") >= 0) {
                    if(fileName.indexOf("_mixMHC2pred_conf.txt") >= 0) {
                        targetFile = params.fullOutput ? "11_Fusions/NeoFuse/Class_II/" + file(fileName).getName() : ""
                    } else {
                        targetFile = "11_Fusions/NeoFuse/Class_II/" + file(fileName).getName()
                    }
                } else if(fileName.indexOf("STAR") >= 0) {
                    if(fileName.indexOf("Aligned.sortedByCoord.out.bam") >= 0) {
                        targetFile = "02_alignments/" + file(fileName).getName().replace(".Aligned.sortedByCoord.out", "_RNA.Aligned.sortedByCoord.out")
                    }
                } else if(fileName.indexOf("TPM") >= 0) {
                    targetFile = "04_expression/" + file(fileName).getName()
                } else {
                    targetFile = "11_Fusions/" + fileName
                }
                return "$targetFile"
        },
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        val(meta_RNA_ori),
        path(reads_RNA),
        path(hla_types),
        val(meta_DNA_ori),
        path(SVvcf)
    ) from reads_tumor_neofuse_ch
        .join(hlas_neoFuse, by: 0)
        .join(MantaSomaticIndels_out_NeoFuse_in_ch0, by: 0)

    path STARidx from file(reference.STARidx)
    path RefFasta from file(reference.RefFasta)
    path AnnoFile from file(reference.AnnoFile)

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_NeoFuse_MHC_Class_I_filtered.tsv"),
        path("${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_NeoFuse_MHC_Class_I_unfiltered.tsv"),
        path("${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_NeoFuse_MHC_Class_II_filtered.tsv"),
        path("${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_NeoFuse_MHC_Class_II_unfiltered.tsv")
    ) into (
        Neofuse_results_ch0,
        Neofuse_results_ch1
    )
    tuple(
        val(meta),
        path("${meta.sampleName}/TPM/${meta.sampleName}.tpm.txt")
    ) into tpm_file
    tuple(
        val(meta),
        path("${meta.sampleName}/STAR/${meta.sampleName}.Aligned.sortedByCoord.out.{bam,bam.bai}")
    ) into star_bam_file
    path("${meta.sampleName}/Arriba/*")
    path("${meta.sampleName}/Custom_HLAs/*"), optional: true
    path("${meta.sampleName}/LOGS/*")


    script:
    def reads = (meta_RNA_ori.libType == "SE") ? "-1 " + reads_RNA[0] : "-1 " + reads_RNA[0] + " -2 " + reads_RNA[1]
    def sv_options = (meta_DNA_ori.libType == "SE") ? "" : "-v ${SVvcf[0]}"
    """
    NeoFuse_single ${reads} \\
        -d ${meta.sampleName} \\
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
        ${sv_options} \\
        -k true

    mv ${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_MHCI_filtered.tsv ${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_NeoFuse_MHC_Class_I_filtered.tsv
    mv ${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_MHCI_unfiltered.tsv ${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_NeoFuse_MHC_Class_I_unfiltered.tsv
    mv ${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_unsupported.txt ${meta.sampleName}/NeoFuse/MHC_I/${meta.sampleName}_NeoFuse_MHC_Class_I_unsupported.txt
    mv ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_MHCII_filtered.tsv ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_NeoFuse_MHC_Class_II_filtered.tsv
    mv ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_MHCII_unfiltered.tsv ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_NeoFuse_MHC_Class_II_unfiltered.tsv
    mv ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_unsupported.txt ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_NeoFuse_MHC_Class_II_unsupported.txt
    mv ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_mixMHC2pred_conf.txt ${meta.sampleName}/NeoFuse/MHC_II/${meta.sampleName}_NeoFuse_MHC_Class_II_mixMHC2pred_conf.txt

    """
}

process publish_NeoFuse {
    tag "${meta.sampleName}"

    publishDir "$params.outputDir/neoantigens/${meta.sampleName}/",
    saveAs: {
        fileName ->
            if(fileName.indexOf("_MHCI_") >= 0) {
                targetFile = "Class_I/Fusions/" + file(fileName).getName()
            } else if(fileName.indexOf("_MHCII_") >= 0) {
                targetFile = "Class_II/Fusions/" + file(fileName).getName()
            } else {
                targetFile = fileName
            }
            return targetFile
    },
    mode: "copy",
    enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(MHC_I_filtered),
        path(MHC_I_unfiltered),
        path(MHC_II_filtered),
        path(MHC_II_unfiltered)
    ) from Neofuse_results_ch0

    output:
    path(MHC_I_filtered)
    path(MHC_I_unfiltered)
    path(MHC_II_filtered)
    path(MHC_II_unfiltered)

    script:
    """
    echo "Done"
    """
}

/*
Add the gene ID (required by vcf-expression-annotator) to the TPM file
*/
process add_geneID {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    input:
    tuple(
        val(meta),
        path(tpm)
    ) from tpm_file
    path AnnoFile from file(reference.AnnoFile)

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("*.tpm_final.txt")
    ) into final_gene_expression_file

    script:
    """
    NameToID.py -i ${tpm} -a ${AnnoFile} -o .
    """
}

/*
Add gene expression info to the VEP annotated, phased VCF file
*/

// branch into samples with and without RNA
VEPvcf_out_ch2 = VEPvcf_out_ch2.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}.branch {
    RNA   : it[0].have_RNA == true
    no_RNA: it[0].have_RNA == false
}

process gene_annotator {

    tag "${meta.sampleName}"

    label 'pVACtools'

    input:
    tuple(
        val(meta),
        path(vep_somatic_vcf_gz),
        path(final_tpm),
        path(RNA_bam),
    ) from VEPvcf_out_ch2.RNA
        .join(final_gene_expression_file, by: 0)
        .join(star_bam_file, by: 0)

    tuple(
        path(RefFasta),
        path(RefIdx),
    ) from Channel.value(
        [ reference.RefFasta,
        reference.RefIdx ]
    )

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_vep_somatic_gx.{vcf.gz,vcf.gz.tbi}")
    ) into (
        vcf_vep_ex_gz,
        gene_annotator_out_mixMHC2pred_ch0,
        generate_protein_fasta_tumor_vcf_ch0
    )

    script:
    """
    vcf-expression-annotator \\
        -i GeneID \\
        -e TPM \\
        -s ${meta.sampleName}_tumor \\
        ${vep_somatic_vcf_gz[0]} ${final_tpm} \\
        custom gene \\
        -o ./${meta.sampleName}_vep_somatic_gx_tmp.vcf
    bgzip -f ${meta.sampleName}_vep_somatic_gx_tmp.vcf
    tabix -p vcf ${meta.sampleName}_vep_somatic_gx_tmp.vcf.gz

    vt decompose \\
        -s ${meta.sampleName}_vep_somatic_gx_tmp.vcf.gz \\
        -o ${meta.sampleName}_vep_somatic_gx_dec_tmp.vcf.gz

    bam_readcount_helper.py \\
        ${meta.sampleName}_vep_somatic_gx_dec_tmp.vcf.gz \\
        ${meta.sampleName}_tumor \\
        ${RefFasta} \\
        ${RNA_bam[0]} \\
        ./

    vcf-readcount-annotator \\
        -s ${meta.sampleName}_tumor \\
        -t snv \\
        -o ${meta.sampleName}_vep_somatic_gx_dec_snv_rc_tmp.vcf \\
        ${meta.sampleName}_vep_somatic_gx_dec_tmp.vcf.gz \\
        ${meta.sampleName}_tumor_bam_readcount_snv.tsv \\
        RNA

    vcf-readcount-annotator \\
        -s ${meta.sampleName}_tumor \\
        -t indel \\
        -o ${meta.sampleName}_vep_somatic_gx.vcf \\
        ${meta.sampleName}_vep_somatic_gx_dec_snv_rc_tmp.vcf \\
        ${meta.sampleName}_tumor_bam_readcount_indel.tsv \\
        RNA

    bgzip -f ${meta.sampleName}_vep_somatic_gx.vcf
    tabix -p vcf ${meta.sampleName}_vep_somatic_gx.vcf.gz
    """
}

// re-add un-annotated vcfs for samples without RNAseq data
(vcf_vep_ex_gz, generate_protein_fasta_tumor_vcf_ch0, gene_annotator_out_mixMHC2pred_ch0) = vcf_vep_ex_gz.mix(VEPvcf_out_ch2.no_RNA).into(3)


// IEDB installation
iedb_chck_file_name = ".iedb_install_ok.chck"
iedb_chck_file = file(params.databases.IEDB_dir + "/" + iedb_chck_file_name)

// TODO: check if unittests for IEDB_MHC_II work in future version, and uncomment call of configure.py
if(!iedb_chck_file.exists() || iedb_chck_file.isEmpty()) {

    log.warn "WARNING: IEDB yet not installed, starting installation. This may take a while..."

    process install_IEDB {

        tag "Install IEDB"

        publishDir "${params.databases.IEDB_dir}",
            mode: "copy"

        label 'pVACtools'

        input:
        val(iedb_MHCI_url) from Channel.value(params.IEDB_MHCI_url)
        val(iedb_MHCII_url) from Channel.value(params.IEDB_MHCII_url)

        output:
        path("${iedb_chck_file_name}") into iedb_install_out_ch

        script:
        def mhci_file = iedb_MHCI_url.split("/")[-1]
        def mhcii_file = iedb_MHCII_url.split("/")[-1]
        """
        export TMPDIR=${params.tmpDir}

        CWD=`pwd`
        cd /opt/iedb/
        rm -f $mhci_file
        wget $iedb_MHCI_url
        tar -xzvf $mhci_file
        cd mhc_i
        bash -c "./configure"
        cd /opt/iedb/
        rm -f $mhci_file

        rm -f $mhcii_file
        wget $iedb_MHCII_url
        tar -xzvf $mhcii_file
        #### ATTENTION: IEDB_MHC_II-3.1.8.tar.gz "python configure.py"
        ####            returns an assertion error in the unittest needs
        ####            to be fixed, skip unittests for now
        # cd mhc_ii
        # bash -c "python ./configure.py"
        cd /opt/iedb/
        rm $mhcii_file

        export MHCFLURRY_DATA_DIR=/opt/mhcflurry_data
        mhcflurry-downloads fetch

        cd \$CWD
        echo "OK" > ${iedb_chck_file_name}
        """
    }
} else {

    iedb_install_out_ch = Channel.fromPath(iedb_chck_file)

}

/*
Run pVACseq
*/
mkPhasedVCF_out_pVACseq_ch0 = mkPhasedVCF_out_pVACseq_ch0.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

process 'pVACseq' {

    tag "${meta.sampleName}"

    label 'pVACtools'

    input:
    tuple(
        val(meta),
        path(vep_phased_vcf_gz),
        path(anno_vcf),
        val(hla_types),
        val(tumor_purity),
        path(iedb_install_ok)
    ) from mkPhasedVCF_out_pVACseq_ch0
        .join(vcf_vep_ex_gz, by: [0])
        .combine(hlas.splitText(), by: 0)
        .combine(purity_estimate_ch1, by: 0)
        .combine(iedb_install_out_ch)

    output:
    tuple(
        val(meta),
        val("Class_I"),
        path(pvacseq_class_I_out)
    ) optional true into mhcI_out_f

    tuple(
        val(meta),
        val("Class_II"),
        path(pvacseq_class_II_out)
    ) optional true into mhcII_out_f


    script:
    def hla_type = (hla_types - ~/\n/)
    def NetChop = params.use_NetChop ? "--net-chop-method cterm" : ""
    def NetMHCstab = params.use_NetMHCstab ? "--netmhc-stab" : ""
    def phased_vcf_opt = (have_GATK3) ? "-p " + vep_phased_vcf_gz[0] : ""

    def filter_set = params.pVACseq_filter_sets[ "standard" ]

    if (params.pVACseq_filter_sets[ params.pVACseq_filter_set ] != null) {
        filter_set = params.pVACseq_filter_sets[ params.pVACseq_filter_set ]
    } else {
        log.warn "WARNING: pVACseq_filter_set must be one of: standard, relaxed, custom\n" +
            "using standard"
        filter_set = params.pVACseq_filter_sets[ "standard" ]
    }


    if(!have_GATK3) {

        log.warn "WARNING: GATK3 not installed! Have no readbacked phased VCF:\n" +
            "You should manually review the sequence data for all candidates (e.g. in IGV) for proximal variants and\n" +
            " either account for these manually, or eliminate these candidates. Failure to do so may lead to inclusion\n" +
            " of incorrect peptide sequences."

    }

    if (have_RNA_tag_seq) {
        filter_set = filter_set.replaceAll(/--trna-vaf\s+\d+\.{0,1}\d*/, "--trna-vaf 0.0")
        filter_set = filter_set.replaceAll(/--trna-cov\s+\d+/, "--trna-cov 0")
    }

    pvacseq_class_I_out = [ "MHC_Class_I/" + meta.sampleName + "_tumor_" + hla_type + ".filtered.tsv",
                            "MHC_Class_I/" + meta.sampleName + "_tumor_" + hla_type + ".all_epitopes.tsv"]
    pvacseq_class_II_out = [ "MHC_Class_II/" + meta.sampleName + "_tumor_" + hla_type + ".filtered.tsv",
                             "MHC_Class_II/" + meta.sampleName + "_tumor_" + hla_type + ".all_epitopes.tsv"]

    """
    pvacseq run \\
        --iedb-install-directory /opt/iedb \\
        -t ${task.cpus} \\
        ${phased_vcf_opt} \\
        -e1 ${params.mhci_epitope_len} \\
        -e2 ${params.mhcii_epitope_len} \\
        --normal-sample-name ${meta.sampleName}_normal \\
        --tumor-purity ${tumor_purity} \\
        ${NetChop} \\
        ${NetMHCstab} \\
        ${filter_set} \\
        ${anno_vcf[0]} ${meta.sampleName}_tumor ${hla_type} ${params.epitope_prediction_tools} ./

    if [ -e ./MHC_Class_I/${meta.sampleName}_tumor.filtered.tsv ]; then
        mv ./MHC_Class_I/${meta.sampleName}_tumor.filtered.tsv ./MHC_Class_I/${meta.sampleName}_tumor_${hla_type}.filtered.tsv
    fi
    if [ -e ./MHC_Class_I/${meta.sampleName}_tumor.all_epitopes.tsv ]; then
        mv ./MHC_Class_I/${meta.sampleName}_tumor.all_epitopes.tsv ./MHC_Class_I/${meta.sampleName}_tumor_${hla_type}.all_epitopes.tsv
    fi
    if [ -e ./MHC_Class_II/${meta.sampleName}_tumor.filtered.tsv ]; then
        mv ./MHC_Class_II/${meta.sampleName}_tumor.filtered.tsv ./MHC_Class_II/${meta.sampleName}_tumor_${hla_type}.filtered.tsv
    fi
    if [ -e ./MHC_Class_II/${meta.sampleName}_tumor.all_epitopes.tsv ]; then
        mv ./MHC_Class_II/${meta.sampleName}_tumor.all_epitopes.tsv ./MHC_Class_II/${meta.sampleName}_tumor_${hla_type}.all_epitopes.tsv
    fi
    """
}

// combine class_I and class_II, and make list with filtered and list with all
pVACseq_out = mhcI_out_f.mix(mhcII_out_f).groupTuple(by:[0,1]).map{
    meta, mhc_class, files -> {
        def flat_files = files.flatten()
        def filtered = []
        def all = []
        flat_files.eachWithIndex{ v, ix -> ( ix& 1 ? all : filtered ) << v}
        return [meta, mhc_class, filtered, all]
    }
}

process concat_pVACseq_files {

    tag "${meta.sampleName}"

    label 'pVACtools'

    publishDir "$params.outputDir/analyses/${meta.sampleName}/12_pVACseq/MHC_${mhc_class}/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        val(mhc_class),
        path("*pVACseq_filtered_files"),
        path("*pVACseq_all_files")
    ) from pVACseq_out

    output:
    tuple(
        val(meta),
        val(mhc_class),
        path(out_files)
    ) into pVACseq_out_concat

    script:
    out_files = [ meta.sampleName + "_MHC_" + mhc_class + "_filtered.tsv",
                  meta.sampleName + "_MHC_" + mhc_class + "_all_epitopes.tsv" ]
    """
    concat_pvacseq.py --pattern "pVACseq_filtered_files" --output ${meta.sampleName}_MHC_${mhc_class}_filtered.tsv
    concat_pvacseq.py --pattern "pVACseq_all_files" --output ${meta.sampleName}_MHC_${mhc_class}_all_epitopes.tsv
    """
}
// pVACseq_out_concat.view()
(pVACseq_out_concat_ch0, pVACseq_out_concat_ch1, pVACseq_out_concat_ch2) = pVACseq_out_concat.into(3)

add_CCF_ch = pVACseq_out_concat_ch2.map{
    meta, mhc_class, files -> {
        return [meta, files ]
    }
}.transpose()


igs_ch = pVACseq_out_concat_ch0.map{
    meta, mhc_class, files -> {
        if (mhc_class == "Class_I") {
            return [meta, files[0]]
        }
    }
}

(aggregated_reports_ch0, csin_ch0) = pVACseq_out_concat_ch1.map{
    meta, mhc_class, files -> {
        return [meta, mhc_class, files[1]]
    }
}.into(2)


process aggregated_reports {

    tag "${meta.sampleName}"

    label 'pVACtools'

    publishDir "$params.outputDir/analyses/${meta.sampleName}/12_pVACseq/MHC_${mhc_class}/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        val(mhc_class),
        path(epitope_file),
        val(tumor_purity)
    ) from aggregated_reports_ch0
        .combine(purity_estimate_ch2, by: 0)

    output:
    path("${meta.sampleName}_MHC_${mhc_class}_all_aggregated.tsv")

    script:
    """
    pvacseq generate_aggregated_report \\
        $epitope_file \\
        ${meta.sampleName}_MHC_${mhc_class}_all_aggregated.tsv
    """
}

generate_protein_fasta_phased_vcf_ch0 = generate_protein_fasta_phased_vcf_ch0.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}


process 'pVACtools_generate_protein_seq' {

    tag "${meta.sampleName}"

    label 'pVACtools'

    publishDir "$params.outputDir/analyses/${meta.sampleName}/06_proteinseq/",
    mode: publishDirMode,
    enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(vep_phased_vcf_gz),
        path(vep_tumor_vcf_gz)
    ) from generate_protein_fasta_phased_vcf_ch0
        .combine(generate_protein_fasta_tumor_vcf_ch0, by:[0])

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_long_peptideSeq.fasta")
    ) optional true into pVACtools_generate_protein_seq

    script:
    def phased_vcf_opt = (have_GATK3) ? "-p " + vep_phased_vcf_gz[0] : ""

    if(!have_GATK3) {

        log.warn "WARNING: GATK3 not installed! Have no readbacked phased VCF:\n" +
        "You should manually review the sequence data for all candidates (e.g. in IGV) for proximal variants and\n" +
        " either account for these manually, or eliminate these candidates. Failure to do so may lead to inclusion\n" +
        " of incorrect peptide sequences."
    }

    """
    pvacseq generate_protein_fasta \\
        ${phased_vcf_opt} \\
        -s ${meta.sampleName}_tumor \\
        -d full \\
        ${vep_tumor_vcf_gz[0]} \\
        31 \\
        ${meta.sampleName}_long_peptideSeq.fasta
    """
}

hlahd_mixMHC2_pred_ch0 = hlahd_mixMHC2_pred_ch0.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}


if(have_HLAHD) {
    process 'pepare_mixMHC2_seq' {

        label 'nextNEOpiENV'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/13_mixMHC2pred/processing/",
            mode: publishDirMode,
            enabled: params.fullOutput

        input:
        tuple(
            val(meta),
            path(long_peptideSeq_fasta),
            path(hlahd_allel_file)
        ) from pVACtools_generate_protein_seq
            .join(hlahd_mixMHC2_pred_ch0, by:0)

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        tuple(
            val(meta),
            path("${meta.sampleName}_peptides.fasta")
        ) optional true into pepare_mixMHC2_seq_out_ch0
        path("${meta.sampleName}_mixMHC2pred.txt") optional true into pepare_mixMHC2_seq_out_ch1
        path("${meta.sampleName}_unsupported.txt") optional true
        path("${meta.sampleName}_mixMHC2pred_conf.txt") optional true

        script:
        def supported_list = baseDir.toRealPath() + "/assets/hlaii_supported.txt"
        def model_list     = baseDir.toRealPath() + "/assets/hlaii_models.txt"
        """
        pepChopper.py \\
            --pep_len ${params.mhcii_epitope_len.split(",").join(" ")} \\
            --fasta_in ${long_peptideSeq_fasta} \\
            --fasta_out ${meta.sampleName}_peptides.fasta
        HLAHD2mixMHC2pred.py \\
            --hlahd_list ${hlahd_allel_file} \\
            --supported_list ${supported_list} \\
            --model_list ${model_list} \\
            --output_dir ./ \\
            --sample_name ${meta.sampleName}
        """
    }

    mixmhc2pred_chck_file = file(workflow.workDir + "/.mixmhc2pred_install_ok.chck")
    mixmhc2pred_target = workflow.workDir + "/MixMHC2pred"
    if(( ! mixmhc2pred_chck_file.exists() || mixmhc2pred_chck_file.isEmpty()) && params.MiXMHC2PRED == "") {
        process install_mixMHC2pred {

            tag 'install mixMHC2pred'

            // do not cache
            cache false

            output:
            path(".mixmhc2pred_install_ok.chck") into mixmhc2pred_chck_ch

            script:
            """
            curl -sLk ${params.MiXMHC2PRED_url} -o mixmhc2pred.zip && \\
            unzip -o mixmhc2pred.zip -d ${mixmhc2pred_target} && \\
            echo "OK" > .mixmhc2pred_install_ok.chck && \\
            cp -f .mixmhc2pred_install_ok.chck ${mixmhc2pred_chck_file}
            """
        }
    } else if (( ! mixmhc2pred_chck_file.exists() || mixmhc2pred_chck_file.isEmpty()) && params.MiXMHC2PRED != "") {
        process link_mixMHC2pred {

            tag 'link mixMHC2pred'

            // do not cache
            cache false

            output:
            path(".mixmhc2pred_install_ok.chck") into mixmhc2pred_chck_ch

            script:
            """
            ln -s ${params.MiXMHC2PRED} ${mixmhc2pred_target} && \\
            echo "OK" > .mixmhc2pred_install_ok.chck && \\
            cp -f .mixmhc2pred_install_ok.chck ${mixmhc2pred_chck_file}
            """
        }
    } else {
        mixmhc2pred_chck_ch = Channel.fromPath(mixmhc2pred_chck_file)
    }

    process mixMHC2pred {

        label 'nextNEOpiENV'

        tag "${meta.sampleName}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/13_mixMHC2pred",
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(mut_peps),
            path(vep_somatic_gx_vcf_gz),
            path(mixmhc2pred_chck_file)
        ) from pepare_mixMHC2_seq_out_ch0
            .join(gene_annotator_out_mixMHC2pred_ch0, by: [0])
            .combine(mixmhc2pred_chck_ch)
        val allelesFile from pepare_mixMHC2_seq_out_ch1

        val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

        output:
        path("${meta.sampleName}_mixMHC2pred_all.tsv") optional true
        path("${meta.sampleName}_mixMHC2pred_filtered.tsv") optional true

        script:
        def alleles = file(allelesFile).readLines().join(" ")

        if(alleles.length() > 0)
            """
            ${mixmhc2pred_target}/MixMHC2pred_unix \\
                -i ${mut_peps} \\
                -o ${meta.sampleName}_mixMHC2pred.tsv \\
                -a ${alleles}
            parse_mixMHC2pred.py \\
                --vep_vcf ${vep_somatic_gx_vcf_gz[0]} \\
                --pep_fasta ${mut_peps} \\
                --mixMHC2pred_result ${meta.sampleName}_mixMHC2pred.tsv \\
                --out ${meta.sampleName}_mixMHC2pred_all.tsv \\
                --sample_name ${meta.sampleName}_tumor \\
                --normal_name ${meta.sampleName}_normal
            awk \\
                '{
                    if (\$0 ~ /\\#/) { print }
                    else { if (\$18 <= 2) { print } }
                }' ${meta.sampleName}_mixMHC2pred_all.tsv > ${meta.sampleName}_mixMHC2pred_filtered.tsv
            """
        else
            """
            true
            """
    }
}

// add CCF clonality to neoepitopes result files

// adjust channel
Clonality_out_ch0 = Clonality_out_ch0.map{
    meta_ori, out_file, ascatOK, sequenzaOK ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file, ascatOK, sequenzaOK]
}

process addCCF {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/neoantigens/${meta.sampleName}/",
        saveAs: {
            fileName ->
                targetFile = fileName
                if(fileName.indexOf("_MHC_Class_I_") >= 0) {
                    targetFile = "Class_I/" + file(fileName).getName()
                } else if(fileName.indexOf("_MHC_Class_II_") >= 0) {
                    targetFile = "Class_II/" + file(fileName).getName()
                }
                return "$targetFile"
        },
        mode: publishDirMode,
        enabled: params.fullOutput

    input:
    tuple(
        val(meta),
        path(epitopes),
        path(CCF),
        val(ascatOK),
        val(sequenzaOK)
    ) from add_CCF_ch
        .combine(Clonality_out_ch0, by: 0)

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        val("pVACseq"),
        val(f_type),
        val(mhc_class),
        path(outfile)
    ) into addCCF_out_ch
    path("INFO.txt") optional true

    script:
    outfile = (ascatOK || sequenzaOK) ? epitopes.baseName + "_ccf.tsv" : epitopes
    f_type = (epitopes.baseName.indexOf("filtered") >= 0) ? "filtered" : "unfiltered"
    mhc_class = (epitopes.baseName.indexOf("MHCII") >= 0) ? "II" : "I"
    if (ascatOK || sequenzaOK)
        """
        add_CCF.py \\
            --neoepitopes ${epitopes} \\
            --ccf ${CCF} \\
            --outfile ${outfile}
        """
    else
        """
        echo "WARNING: neither ASCAT nor Sequenza produced results: clonality information missing" > INFO.txt
        """
}

// convert channel in to (id, file) tuples
Neofuse_results_ch1 = Neofuse_results_ch1.map{ it ->
                                                meta = it[0]
                                                l = []
                                                for ( f in it[1..-1] ) {
                                                    f_type = (f.getName().indexOf("unfiltered") >= 0) ? "unfiltered" : "filtered"
                                                    mhc_class = (f.getName().indexOf("MHC_Class_II") >= 0) ? "II" : "I"
                                                    l.add(tuple(meta, "NeoFuse", f_type, mhc_class, f))
                                                }
                                                return l
                                            }
                                            .flatten()
                                            .collate(5)

epitopes_fasta_in_ch = addCCF_out_ch.mix(Neofuse_results_ch1)
(epitopes_fasta_in_ch, epitopes_protein_match_in_ch) = epitopes_fasta_in_ch.into(2)

process make_epitopes_fasta {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    cache 'lenient'

    input:
    tuple(
        val(meta),
        val(caller),
        val(f_type),
        val(mhc_class),
        path(epitopes)
    ) from epitopes_fasta_in_ch

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        val(caller),
        val(f_type),
        val(mhc_class),
        path(outfile)
    ) into make_epitopes_fasta_out_ch

    script:
    outfile = epitopes.baseName + "_epitopes.fasta"
    """
    make_peptide_fasta.py \\
        --epitope_caller ${caller} \\
        --fasta ${outfile} \\
        --epitope_file ${epitopes}
    """
}

process blast_epitopes {

    label 'Blast'

    tag "${meta.sampleName}"

    cache 'lenient'

    input:
    tuple(
        val(meta),
        val(caller),
        val(f_type),
        val(mhc_class),
        path(epitopes_fasta)
    ) from make_epitopes_fasta_out_ch

    path(blastdb) from Channel.value(reference.ProteinBlastDBdir)

    output:
    tuple(
        val(meta),
        val(caller),
        val(f_type),
        val(mhc_class),
        path(outfile)
    ) into blast_epitopes_out_ch

    script:
    outfile = epitopes_fasta.baseName + "_blast.tsv"
    """
    blastp -task blastp-short \\
        -db ${blastdb}/${params.references.ProteinBlastDBname} \\
        -query ${epitopes_fasta} \\
        -out ${outfile} \\
        -outfmt "6 qseqid sseqid qlen length nident qseq" \\
        -num_alignments 5 \\
        -num_threads ${task.cpus} \\
        -comp_based_stats 0 \\
        -ungapped \\
        -seg no
    """
}

process add_blast_hits {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/neoantigens/${meta.sampleName}/",
        saveAs: {
            fileName ->
                targetFile = fileName
                if(fileName.indexOf("NeoFuse_MHC_Class_I_") >= 0) {
                    targetFile = "Class_I/Fusions/" + file(fileName).getName()
                } else if(fileName.indexOf("NeoFuse_MHC_Class_II_") >= 0) {
                    targetFile = "Class_II/Fusions/" + file(fileName).getName()
                } else if(fileName.indexOf("${meta.sampleName}_MHC_Class_I_") >= 0) {
                    targetFile = "Class_I/" + file(fileName).getName()
                } else if(fileName.indexOf("${meta.sampleName}_MHC_Class_II_") >= 0) {
                    targetFile = "Class_II/" + file(fileName).getName()
                }

                return "$targetFile"
        },
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        val(caller),
        val(f_type),
        val(mhc_class),
        path(blast_result),
        path(epitopes)
    ) from blast_epitopes_out_ch
        .join(epitopes_protein_match_in_ch, by: [0,1,2,3])

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    path("*_ref_match.tsv")

    script:
    """
    parse_blast_result.py \\
        --blast_result ${blast_result} \\
        --epitope_file ${epitopes} \\
        --epitope_caller ${caller}
    """
}


/*
  Immunogenicity scoring
*/

process csin {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "$params.outputDir/analyses/${meta.sampleName}/14_CSiN/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        val(mhc_class),
        path(all_epitopes),
    ) from csin_ch0
        .groupTuple()

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    tuple(
        val(meta),
        path("${meta.sampleName}_CSiN.tsv")
    ) into sample_info_csin

    script:
    def mhc_i_idx = mhc_class.indexOf("Class_I")
    def mhc_ii_idx = mhc_class.indexOf("Class_II")
    """
    CSiN.py --MHCI_tsv ${all_epitopes[mhc_i_idx]} \\
        --MHCII_tsv ${all_epitopes[mhc_ii_idx]} \\
        --rank $params.csin_rank \\
        --ic50 $params.csin_ic50 \\
        --gene_exp $params.csin_gene_exp \\
        --output ./${meta.sampleName}_CSiN.tsv
    """
}

igs_chck_file = file(workflow.workDir + "/.igs_install_ok.chck")
igs_target = workflow.workDir + "/IGS"
if(( ! igs_chck_file.exists() || igs_chck_file.isEmpty()) && params.IGS == "") {
    process install_IGS {

        tag 'install IGS'

        // do not cache
        cache false

        output:
        val("OK") into igs_chck_ch
        path(".igs_install_ok.chck")

        script:
        """
        mkdir -p ${igs_target} && \\
        curl -sLk ${params.IGS_script_url} -o ${igs_target}/NeoAg_immunogenicity_predicition_GBM.R && \\
        curl -sLk ${params.IGS_model_url} -o ${igs_target}/Final_gbm_model.rds && \\
        patch -p0 ${igs_target}/NeoAg_immunogenicity_predicition_GBM.R ${baseDir}/assets/NeoAg_immunogenicity_predicition_GBM.patch && \\
        chmod +x ${igs_target}/NeoAg_immunogenicity_predicition_GBM.R  && \\
        echo "OK" > .igs_install_ok.chck && \\
        cp -f .igs_install_ok.chck ${igs_chck_file}
        """
    }
} else if (( ! igs_chck_file.exists() || igs_chck_file.isEmpty()) && params.IGS != "") {
    process link_IGS {

        tag 'link IGS'

        // do not cache
        cache false

        output:
        val("OK") into igs_chck_ch
        path(".igs_install_ok.chck")

        script:
        """
        ln -s ${params.IGS} ${igs_target} && \\
        echo "OK" > .igs_install_ok.chck && \\
        cp -f .igs_install_ok.chck ${igs_chck_file}
        """
    }
} else {
    igs_chck_ch = Channel.value("OK")
}


process immunogenicity_scoring {

    label 'IGS'

    tag "${meta.sampleName}"

    // TODO: check why sometimes this fails: workaround ignore errors
    errorStrategy 'ignore'

    publishDir "$params.outputDir/analyses/${meta.sampleName}/14_IGS/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(pvacseq_file)
    ) from igs_ch
    val (igs_install_chck) from igs_chck_ch

    output:
    path("${meta.sampleName}_Class_I_immunogenicity.tsv")

    script:
    """
    get_epitopes.py \\
        --pvacseq_out $pvacseq_file \\
        --sample_id $meta.sampleName \\
        --output ./${meta.sampleName}_epitopes.tsv
    NR_EPI=`wc -l ./${meta.sampleName}_epitopes.tsv | cut -d" " -f 1`
    if [ \$NR_EPI -gt 1 ]; then
        ${igs_target}/NeoAg_immunogenicity_predicition_GBM.R \\
            ./${meta.sampleName}_epitopes.tsv ./${meta.sampleName}_temp_immunogenicity.tsv \\
            ${igs_target}/Final_gbm_model.rds
        immuno_score.py \\
            --pvacseq_tsv $pvacseq_file \\
            --score_tsv ${meta.sampleName}_temp_immunogenicity.tsv \\
            --output ${meta.sampleName}_Class_I_immunogenicity.tsv
    fi
    """
}

if(params.TCR) {

    mixcr_chck_file = file(baseDir + "/bin/.mixcr_install_ok.chck")
    mixcr_target = baseDir + "/bin/"
    if(!mixcr_chck_file.exists() && params.MIXCR == "") {
        process install_mixcr {

            tag 'install mixcr'

            // do not cache
            cache false

            output:
            path(".mixcr_install_ok.chck") into mixcr_chck_ch

            script:
            """
            curl -sLk ${params.MIXCR_url} -o mixcr.zip && \\
            unzip -o mixcr.zip && \\
            chmod +x mixcr && \\
            cp -f mixcr ${mixcr_target} && \\
            cp -f mixcr.jar ${mixcr_target} && \\
            cp -f ${params.MIXCR_lic} ${mixcr_target}/mi.license && \\
            touch .mixcr_install_ok.chck && \\
            cp -f .mixcr_install_ok.chck ${mixcr_chck_file}
            """
        }
    } else if (!mixcr_chck_file.exists() && params.MIXCR != "") {
        process link_mixcr {

            tag 'link mixcr'

            // do not cache
            cache false

            output:
            path(".mixcr_install_ok.chck") into mixcr_chck_ch

            script:
            """
            ln -s ${params.MIXCR}/mixcr ${mixcr_target} && \\
            ln -s ${params.MIXCR}/mixcr.jar ${mixcr_target} && \\
            ln -s ${params.MIXCR_lic} ${mixcr_target}/mi.license && \\
            touch .mixcr_install_ok.chck && \\
            cp -f .mixcr_install_ok.chck ${mixcr_chck_file}
            """
        }
    } else {
        mixcr_chck_ch = Channel.fromPath(mixcr_chck_file)
    }

    process mixcr {

        tag "${meta.sampleName} : ${meta.sampleType}"

        publishDir "$params.outputDir/analyses/${meta.sampleName}/15_BCR_TCR",
            saveAs: {
                filename ->
                    if (filename.indexOf(".tsv") == -1 && params.fullOutput) {
                        return "extra_output/$filename"
                    } else if (filename.indexOf(".tsv") == -1 && ! params.fullOutput) {
                        return null
                    } else {
                        return "$filename"
                    }
            },
            mode: publishDirMode

        input:
        tuple(
            val(meta),
            path(reads),
            path(mixcr_chck_file)
        ) from reads_mixcr_ch
            .combine(mixcr_chck_ch)

        output:
        tuple(
            val(meta),
            path("${procSampleName}*.tsv"),
        )

        script:
        def starting_material = (meta.sampleType == "tumor_RNA") ? "rna" : "dna"
        def libtype = (meta.sampleType == "tumor_RNA") ? "rna-seq" : "exome-seq"
        procSampleName = meta.sampleName + "_" + meta.sampleType + "_mixcr"
        """
        ${baseDir}/bin/mixcr analyze ${libtype} \\
            --threads ${task.cpus} \\
            --species hsa \\
            --${starting_material} \\
            $reads \\
            ${procSampleName}
        """
    }
}


// adjust channel
sample_info_tmb = sample_info_tmb.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}
sample_info_tmb_coding = sample_info_tmb_coding.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

process collectSampleInfo {

    label 'nextNEOpiENV'

    tag "${meta.sampleName}"

    publishDir "${params.outputDir}/neoantigens/${meta.sampleName}/",
        mode: publishDirMode

    input:
    tuple(
        val(meta),
        path(csin),
        path(tmb),
        path(tmb_coding)
    ) from sample_info_csin
        .join(sample_info_tmb, by: 0)
        .join(sample_info_tmb_coding, by: 0)

    val (nextNEOpiENV_setup) from nextNEOpiENV_setup_ch0

    output:
    path("${meta.sampleName}_sample_info.tsv")

    script:
    """
    mkSampleInfo.py \\
        --sample_name ${meta.sampleName} \\
        --csin ${csin} \\
        --tmb ${tmb} \\
        --tmb_coding ${tmb_coding} \\
        --out ${meta.sampleName}_sample_info.tsv
    """
}

/*
***********************************
*  Generate final multiQC output  *
***********************************
*/

ch_fastqc = ch_fastqc.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

ch_fastp = ch_fastp.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

ch_fastqc_trimmed = ch_fastqc_trimmed.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}

alignmentMetrics_ch = alignmentMetrics_ch.map{
    meta_ori, out_file ->
        def meta = meta_ori.clone()
        meta.keySet().removeAll(['sampleType', 'libType'])
        return [meta, out_file]
}


process multiQC {

    label 'multiqc'

    tag "${meta.sampleName}"

    publishDir "${params.outputDir}/analyses/${meta.sampleName}/QC",
        mode: publishDirMode

   input:
    tuple(
        val(meta),
        path(qcfiles)
    ) from ch_fastqc
        .mix(ch_fastp)
        .mix(ch_fastqc_trimmed)
        .mix(alignmentMetrics_ch)
        .transpose()
        .groupTuple()

    output:
    path("multiqc_data/*")
    path("multiqc_report.html")

    script:
    def set_locale = ""
    if(! params.enable_conda && workflow.containerEngine == 'singularity' ) {
        set_locale = "export LC_ALL=C.UTF-8; export LC_ALL=C.UTF-8"
    }
    """
    ${set_locale}
    multiqc .
    """

}


/*
________________________________________________________________________________

                            F U N C T I O N S
________________________________________________________________________________

*/

def checkDir(d, description) {
    myDir = file(d)
    result = myDir.mkdirs()
    if (result) {
        println(description + ": " + myDir.toRealPath())
    } else {
        exit 1, "Cannot create directory: " + myDir
    }
    return myDir.toRealPath()
}

def mkTmpDir(d) {
   return checkDir(d, "tmpDir")
}

def check_iedb_dir(d) {
    checkDir(d, "IEDB_dir")
}

def check_mhcflurry_dir(d) {
    checkDir(d, "MHCFLURRY_dir")
}

def setExomeCaptureKit(captureKit) {
    resources = ['BaitsBed', 'RegionsBed']
    for (r in resources) {
        if (params.exomeCaptureKits[ captureKit ][r] == null) {
            exit 1, "ERROR: \"" + r + "\" file not set for: " + captureKit + "\nPlease check the \"exomeCaptureKits\" resource file settings in conf/resources.config"
        }
    }

    params.references.BaitsBed = params.exomeCaptureKits[ captureKit ].BaitsBed
    params.references.RegionsBed = params.exomeCaptureKits[ captureKit ].RegionsBed
}

def check_resource(resource, resource_type) {

    resource_file = params[resource_type][resource]

    if (resource_file == null) {
        exit 1, "ERROR: Resource file not set for: " + resource + "\nPlease check the \"" + resource_type + "\" resource file settings in conf/resources.config"
    }

    rp = file(resource_file)

    if (rp instanceof java.util.LinkedList) {
        rf = rp
    } else {
        rf = [rp]
    }

    err = 0
    if (rf[0] != null) {
        for (f in rf) {
            err += file(f).exists() ? 0 : 1
        }
    }

    if (err == 0) {
        return(rp)
    } else {
        exit 1, "ERROR: Resource file does not exist: " + resource_file + "\nPlease check the " + resource_type + " resource file settings in conf/resources.config"
    }
}

def defineResources(resource_type, wes, hlahd) {

    // vep check file
    vep_cache_chck_file_name = "." + params.vep_species + "_" + params.vep_assembly + "_" + params.vep_cache_version + "_cache_ok.chck"
    vep_cache_chck_file = file(params.databases.vep_cache + "/" + vep_cache_chck_file_name)

    // define references
    references = ['RefFasta', 'RefIdx', 'RefDict', 'RefChrLen', 'RefChrDir', 'BwaRef',
                  'YaraIndexDNA', 'YaraIndexRNA', 'STARidx', 'AnnoFile', 'ExonsBED', 'acLoci',
                  'acLociGC', 'SequenzaGC', 'ProteinBlastDBdir']
    if (wes) {
        references.addAll(['BaitsBed', 'RegionsBed'])
    }
    if (hlahd != "") {
        references.addAll(['HLAHDFreqData', 'HLAHDGeneSplit', 'HLAHDDict'])
    }
    if(vep_cache_chck_file.exists() && !vep_cache_chck_file.isEmpty()) {
        references.addAll(['VepFasta'])
    }

    // define databases
    databases = ['MillsGold', 'MillsGoldIdx', 'hcSNPS1000G', 'hcSNPS1000GIdx', 'HapMap',
                 'HapMapIdx', 'DBSNP', 'DBSNPIdx', 'GnomAD', 'GnomADIdx', 'GnomADfull', 'GnomADfullIdx',
                 'KnownIndels', 'KnownIndelsIdx', 'vep_cache', 'IEDB_dir', 'MHCFLURRY_dir']

    resources = ['references' : references, 'databases' : databases ]
    resources_files = [:]

    for (r in resources[resource_type]) {
        resources_files[r] = check_resource(r, resource_type)
    }

    return(resources_files)
}

def checkToolAvailable(tool, check, errMode, module=false) {
    def checkResult = false
    def res = ""
    def chckCmd = ""
    def envlist = [];

    if (check == "inPath") {
        if(module) {
            chckCmd = 'module load ' + module + ' && which ' + tool + ' && module unload ' + module
        } else {
            chckCmd = "which " + tool
        }
        def processBuilder = new ProcessBuilder(['/bin/bash', '-c', chckCmd])
        processBuilder.environment().putAll(System.getenv())

        def process = processBuilder.start()
        def out = new StringBuilder()

        process.inputStream.eachLine { line ->
            out.append(line).append('\n')
        }

        process.waitFor()
        res = out.toString().trim()
    }

    if (check == "exists") {
        if (file(tool).exists()) {
            res = tool
        }
    }

    if (res == "") {
        def msg = tool + " not found, please make sure " + tool + " is installed"

        if(errMode == "err") {
            msg = "ERROR: " + msg
            msg = (check == "inPath") ? msg + " and in your \$PATH" : msg
            exit(1, msg)
        } else {
            msg = "Warning: " + msg
            msg = (check == "inPath") ? msg + " and in your \$PATH" : msg
            println("Warning: " + msg)
        }
    } else {
        println("Found " + tool + " at: " + res)
        checkResult = true
    }

    return checkResult
}

def checkCondaChannels() {
    Yaml parser = new Yaml()
    def channels = []
    try {
        def config = parser.load("conda config --show channels".execute().text)
        channels = config.channels
    } catch(NullPointerException | IOException e) {
        log.warn "Could not verify conda channel configuration."
        return
    }

    // Check that all channels are present
    def required_channels = ['conda-forge', 'bioconda', 'defaults']
    def conda_check_failed = !required_channels.every { ch -> ch in channels }

    // Check that they are in the right order
    conda_check_failed |= !(channels.indexOf('conda-forge') < channels.indexOf('bioconda'))
    conda_check_failed |= !(channels.indexOf('bioconda') < channels.indexOf('defaults'))

    if (conda_check_failed) {
        log.warn "=============================================================================\n" +
            "  There is a problem with your Conda configuration!\n\n" +
            "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
            "  Please refer to https://bioconda.github.io/#usage\n" +
            "  NB: The order of the channels matters!\n" +
            "==================================================================================="

        exit 1
    }
}

def showLicense() {

    licenseFile = file(baseDir + "/LICENSE")
    log.info licenseFile.text

    log.info ""
    log.warn "To accept the licence terms, please rerun with '--accept_license'"
    log.info ""

    exit 1
}

def acceptLicense() {
    log.info ""
    log.warn "I have read and accept the licence terms"
    log.info ""

    licenseChckFile = file(baseDir + "/.license_accepted.chck")
    licenseChckFile.text = "License accepted by " + workflow.userName + " on "  + workflow.start

    return true
}

def checkLicense() {
    licenseChckFile = file(baseDir + "/.license_accepted.chck")

    if(!licenseChckFile.exists()) {
        showLicense()
    } else {
        return true
    }
}


def check_seqLibTypes_ok(seqLib_ch) {
    def seqLibs = seqLib_ch.toList().get()
    def pe_count = 0
    def se_count = 0

    def lt_map = [:]

    for (seqLib in seqLibs) {
        if (seqLib[0].sampleType != "tumor_RNA") {
            if (! lt_map.containsKey(seqLib[0].sampleType)) {
                lt_map[seqLib[0].sampleName] = seqLib[2]
            } else {
                if (lt_map[seqLib[0].sampleName] != seqLib[2]) {
                    exit 1, "Please do not mix pe and se for tumor/normal pairs: " + seqLib[0].sampleName + " - Not supported"
                }
            }
        }
    }
    return "OK"
}

// This function removes all keys in key list from meta object at idx 0.
// If keep is true then the original meta object is kept at idx 1
// of the channel values. Note with keep true all values [1..-1] will be
// right shifted by 1
/* NOT working or now: need to check
def remove_from_meta(ch, keys=[], keep=false) {
    ch = ch.map {
            meta, f ->
            def meta_new = meta.clone()
            keys.each{ k ->
                meta_new.remove(k)
            }
            if(keep) {
                return [meta_new, meta.clone(), f]
            } else {
                return [meta_new, f]
            }
        }
    return ch
}
*/

def get_publishMode(d, mode) {
    def req_mode = mode

    if (req_mode != "auto" && req_mode != "link") {
        return mode
    }

    // default to copy
    mode = "copy"

    file(d).mkdirs()

    testFile = file(workflow.workDir + "/.test")
    testFile.write("test")

    testLink = file(d + "/.test")

    // let's see if we can create hard links
    try {
        Files.createLink(testLink, testFile);
        mode = "link"
    } catch (IOException e) {
        if (req_mode == "link") {
            System.err.println("WARNING: using copy as publish mode, reason: " + e)
        }
    }

    testLink.delete()
    testFile.delete()

    return mode
 }

def helpMessage() {
    log.info ""
    log.info "----------------------------"
    log.info "--        U S A G E       "
    log.info "----------------------------"
    log.info ""
    log.info ' nextflow run nextNEOpi.nf -config conf/params.config --batchFile <batchfile.csv> -profile [conda|singularity],[cluster] [-resume]'
    log.info ""
    log.info "-----------------------------------------------------------------------------------------------------------------------------------------"
    log.info ""
    log.info ""
    log.info " Mandatory arguments:"
    log.info " --------------------"
    log.info "--batchFile"
    log.info ""
    log.info "CSV-file, T/N reads, and optionally RNAseq reads:"

    log.info "sampleName,reads1,reads2,sampleType,HLAfile,sex"
    log.info "sample1,reads_s1_t_1.fastq.gz,reads_s1_t_2.fastq.gz,tumor_DNA,,female"
    log.info "sample1,reads_s1_n_1.fastq.gz,reads_s1_n_2.fastq.gz,normal_DNA,,female"
    log.info "sample1,reads_s1_r_1.fastq.gz,reads_s1_r_2.fastq.gz,tumor_RNA,,female"
    log.info "sample2,reads_s2_t_1.fastq.gz,reads_s2_t_2.fastq.gz,tumor_DNA,/data/sample2_hla.txt,male"
    log.info "sample2,reads_s2_n_1.fastq.gz,reads_s2_n_2.fastq.gz,normal_DNA,,male"
    log.info "sample2,reads_s2_r_1.fastq.gz,,tumor_RNA,,male"

    log.info "Note: You can not use samples that have mixed single-end and paired-end DNA reads in tumor and normal."
    log.info "Both, tumor and normal DNA library types need to be either SE or PE for a given sample."

    log.info "Note: in the HLAfile coulumn a user suppiled HLA types file may be specified for a given sample"

    log.info "Note: sex can be XX, female or Female, XY, male or Male. If not specified or \"NA\" the gender is inferred from the data"
    log.info ""
    log.info ""
    log.info ""
    log.info ""
    log.info " All references, databases, software should be edited in the resources.config file"
    log.info ""
    log.info " For further options e.g. BAM input see the README.md, the params.config and process.config files"
    log.info "-----------------------------------------------------------------------------------------------------------------------------------------"
}

// workflow complete
workflow.onComplete {
    // Set up the e-mail variables
    def subject = "[icbi/nextNEOpi] Successful: $workflow.runName"
    if(!workflow.success){
        subject = "[icbi/nextNEOpi] FAILED: $workflow.runName"
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
            log.info "[icbi/nextNEOpi] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << email_txt
            log.info "[icbi/nextNEOpi] Sent summary e-mail to $params.email (mail)"
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

    log.info "[icbi/nextNEOpi] Pipeline Complete! You can find your results in ${params.outputDir}"
}
