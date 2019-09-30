nextflow.preview.dsl = 2

log.info """\
W E S - N F   P I P E L I N E
=============================
genome: ${params.RefFasta}
reads : ${params.reads}
outdir: ${params.outdir}
"""

// import modules
include '' params(params)

genome = file(params.transcritpome)
index


workflow {
  01_BwaSort(reads_ch)



}
