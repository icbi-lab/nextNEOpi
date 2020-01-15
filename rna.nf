#!/usr/bin/env nextflow

log.info ""
log.info " NEXTFLOW ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
log.info "-------------------------------------------------------------------------"
log.info "         NeoFuse + pVACseq   P I P E L I N E             "
log.info "-------------------------------------------------------------------------"
log.info ""
log.info " Predict Neo-antigens in tumor only or tumor + matched control samples"
log.info " using two different tools: \n \t * NeoFuse \n \t * pVACseq \n"
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
if  (params.reads1 != "NO_FILE") log.info " Reads 1: \t\t " + params.reads1
if  (params.reads2 != "NO_FILE") log.info "Reads 2: \t\t " + params.reads2
log.info ""
log.info "Please check --help for further instruction"
log.info "-------------------------------------------------------------------------"

/*
________________________________________________________________________________

                           C O N F I G U R A T I O N
________________________________________________________________________________
*/
if (params.help) exit 0, helpMessage()

// RefFASTA = params.RefFASTA

// if (params.neofuse) {
// 	if (params.reads1 == "NO_FILE") && (params.tsv_file == "NO_FILE") {
// 		exit 1, "No input files defined"
// 	}
// }

reference = defineReference()

if (params.out_ID == "") {
	out_ID = ""
} else {
	out_ID = "-d " + params.out_ID
}

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

if (params.tsv_file != "NO_FILE") {
	Channel
    .fromPath(params.tsv_file)
    .splitCsv(header:true)
    .map { row -> tuple(row.ID, 
    	file(row.Read1),
    	file(row.Read2)) }
    .set { samples_ch }


	process Neofuse_multi {
		publishDir "${params.outputDir}/NeoFuse/", mode: "link"
		
		input:
		tuple (val(ID), file(Read1), file(Read2)) from samples_ch
		file STARidx from file(params.references.STARidx)
		file RefFASTA from file(params.references.RefFASTA)
		file AnnoFile from file(params.references.AnnoFile)

		output:
		file("*${ID}/TPM/*.tpm.txt") into tpm_file
		file("**")

		script:
		"""
		NeoFuse -1 ${Read1} -2 ${Read2} -d ${ID} -o . -m ${params.pepMin_length} -M ${params.pepMax_length} \
		-n 10 -t ${params.IC50_Threshold} -T ${params.rank} -c ${params.conf_lvl} -s ${STARidx} -g ${RefFASTA} -a ${AnnoFile} \
		-N ${params.netMHCpan} -r .
		"""
	}

} else {
	process Neofuse_single {
		publishDir "${params.outputDir}/NeoFuse/", mode: "link"
		
		input:
		file read1 from Channel.fromPath(params.reads1) 
		file read2 from Channel.fromPath(params.reads2)
		file STARidx from file(params.references.STARidx)
		file RefFASTA from file(params.references.RefFASTA)
		file AnnoFile from file(params.references.AnnoFile)

		output:
		file("**/*.tpm.txt") into tpm_file
		file("**")

		script:
		if( params.reads2 != 'NO_FILE' ) 
			"""
			NeoFuse -1 ${read1} -2 ${read2} -o . -m ${params.pepMin_length} -M ${params.pepMax_length} \
			-n 10 -t ${params.IC50_Threshold} -T ${params.rank} -c ${params.conf_lvl} -s ${STARidx} -g ${RefFASTA} -a ${AnnoFile} \
			-N ${params.netMHCpan} -r .
			"""
		else
			script:
			"""
			NeoFuse -1 ${read1} -o . -m ${params.pepMin_length} -M ${params.pepMax_length} \
			-n 10 -t ${params.IC50_Threshold} -T ${params.rank} -c ${params.conf_lvl} -s ${STARidx} -g ${RefFASTA} -a ${AnnoFile} \
			-N ${params.netMHCpan} -r .
			"""
	}
}

process add_geneID {
	
	input:
	file tpm from tpm_file
	file AnnoFile from file(params.references.AnnoFile)

	output:
	file("*.tpm_final.txt") into final_file

	script:
	"""
	python /home/fotakis/myScratch/neoAG_pipeline/RNA_nf/bin/NameToID.py -i ${tpm} -a ${AnnoFile}
	"""
}

process gene_annotator {

	input:
	file vcf from file(params.vcf_file)
	file final_tpm from final_file

	output:
	file("*gx.vcf") into (vcf_vep_gx_ch1, vcf_vep_gx_ch2)

	script:
	"""
	vcf-expression-annotator -i GeneID -e TPM -s TUMOR \
	${vcf} ${final_tpm} custom gene
	"""
}

process bgzip {

	input:
	file(anno_vcf) from vcf_vep_gx_ch1

	output:
	file("*gx.vcf.gz") into vcf_vep_ex_gz
	file("*gx.vcf.gz.tbi") into vcf_vep_ex_gz_tbi

	script:
	"""
	bgzip -f ${anno_vcf}
	tabix -p vcf "${anno_vcf}.gz"
	"""

}

process pVACseq {
	publishDir"${params.outputDir}/pVACseq/${params.sample_name}", mode: "link"

	input:
	file(anno_vcf) from vcf_vep_ex_gz
	file(anno_vcf_tbi) from vcf_vep_ex_gz_tbi
	file(phased_vcf) from file(params.phased_vcf)
	file(phased_vcf_tbi) from file(params.phased_vcf_tbi)

	output:
	file("**")

	script:
	"""
	pvacseq run --normal-sample-name NORMAL \
	${anno_vcf} ${params.sample_name} ${params.HLA_types} \
	${params.baff_tools} \
	. -e ${params.epitope_len} --iedb-install-directory /opt/iedb -t 10 \
	-p ${phased_vcf}
	"""
}

/*
________________________________________________________________________________

                            F U N C T I O N S
________________________________________________________________________________

*/

def checkParamReturnFileReferences(item) {
  params."${item}" = params.references."${item}"
  return Channel.fromPath(params.references."${item}")
}

def defineReference() {
  if (params.references.size() != 3) exit 1, """
  ERROR: Not all References needed found in configuration
  Please check if genome file, genome index file, genome dict file, bwa reference files, vep reference file and interval file is given.
  """
  return [
    'STARidx'   : checkParamReturnFileReferences("STARidx"),
    'AnnoFile'       : checkParamReturnFileReferences("AnnoFile"),
    'RefFASTA'    : checkParamReturnFileReferences("RefFASTA")
  ]
}

def helpMessage() {
  log.info ""
  log.info "----------------------------"
  log.info "--        U S A G E       --"
  log.info "----------------------------"
  log.info ""
  log.info ' nextflow run rna.nf "--reads1" "[--reads2]" "[--tsv_file]" "--outputDir"'
  log.info ""
  log.info "-------------------------------------------------------------------------"
  log.info ""
  log.info "Single-end reads:"
  log.info "---------------------------"
  log.info "--single_end \t\t sets parameter to TRUE (default false)"
  log.info ""
  log.info " Mandatory arguments:"
  log.info " --------------------"
  log.info "--reads1 \t\t reads_{1,2}.fastq \t\t paired-end reads; FASTQ file (can be zipped)"
  log.info ""
  log.info " Optional argument:"
  log.info " ------------------"
  log.info ""
  log.info " Mandatory references:"
  log.info " ---------------------"
  log.info " STARidx \t\t reference.fa \t\t\t STAR index directory"
  log.info " AnnoFile \t\t\t reference.gtf \t\t\t Gene Annotation, GTF file"
  log.info " RefFASTA \t\t reference.fa \t\t Reference Genome; FASTA file"
  log.info ""
  log.info " Required software:"
  log.info " ------------------"
  log.info " pVACtools \t\t\t Version 1.5.4"
  log.info " NeoFuse \t\t\t Version 1.2"
  log.info "-------------------------------------------------------------------------"
}