#!/usr/bin/env Rscript

# Load required library
library(ASCAT)
library(optparse)

# Define command-line options
option_list <- list(
    make_option(c("-t", "--tumourseqfile"),
        type = "character", default = "Tumour.bam",
        help = "Tumour BAM file [default=%default]"
    ),
    make_option(c("-n", "--normalseqfile"),
        type = "character", default = "Normal.bam",
        help = "Normal BAM file [default=%default]"
    ),
    make_option(c("--tumourname"),
        type = "character", default = "Tumour_name",
        help = "Tumour sample name [default=%default]"
    ),
    make_option(c("--normalname"),
        type = "character", default = "Normal_name",
        help = "Normal sample name [default=%default]"
    ),
    make_option(c("--allelecounter_exe"),
        type = "character", default = "/opt/conda/bin/alleleCounter",
        help = "Path to allelecounter executable [default=%default]"
    ),
    make_option(c("--alleles_prefix"),
        type = "character", default = "G1000_alleles_hg19_chr",
        help = "Alleles prefix [default=%default]"
    ),
    make_option(c("--loci_prefix"),
        type = "character", default = "G1000_loci_hg19_chr",
        help = "Loci prefix [default=%default]"
    ),
    make_option(c("--gender"),
        type = "character", default = "XX",
        help = "Gender (XX or XY) [default=%default]"
    ),
    make_option(c("--genomeVersion"),
        type = "character", default = "hg38",
        help = "Genome version (hg19 or hg38) [default=%default]"
    ),
    make_option(c("--nthreads"),
        type = "integer", default = 1,
        help = "Number of threads [default=%default]"
    ),
    make_option(c("--tumourLogR_file"),
        type = "character", default = "Tumor_LogR.txt",
        help = "Tumour LogR output file [default=%default]"
    ),
    make_option(c("--tumourBAF_file"),
        type = "character", default = "Tumor_BAF.txt",
        help = "Tumour BAF output file [default=%default]"
    ),
    make_option(c("--normalLogR_file"),
        type = "character", default = "Germline_LogR.txt",
        help = "Normal LogR output file [default=%default]"
    ),
    make_option(c("--normalBAF_file"),
        type = "character", default = "Germline_BAF.txt",
        help = "Normal BAF output file [default=%default]"
    ),
    make_option(c("--GCcontentfile"),
        type = "character", default = "GC_file.txt",
        help = "GC content file [default=%default]"
    ),
    make_option(c("--replictimingfile"),
        type = "character", default = "RT_file.txt",
        help = "Replication timing file [default=%default]"
    ),
    make_option(c("--gamma"),
        type = "numeric", default = 1,
        help = "Gamma parameter for ascat.runAscat [default=%default]"
    ),
    make_option(c("--output_prefix"),
        type = "character", default = "ASCAT_output",
        help = "Prefix for output files [default=%default]"
    )
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Extract options
tumourseqfile <- opt$tumourseqfile
normalseqfile <- opt$normalseqfile
tumourname <- opt$tumourname
normalname <- opt$normalname
allelecounter_exe <- opt$allelecounter_exe
alleles.prefix <- opt$alleles_prefix
loci.prefix <- opt$loci_prefix
gender <- opt$gender
genomeVersion <- opt$genomeVersion
nthreads <- opt$nthreads
tumourLogR_file <- opt$tumourLogR_file
tumourBAF_file <- opt$tumourBAF_file
normalLogR_file <- opt$normalLogR_file
normalBAF_file <- opt$normalBAF_file
GCcontentfile <- opt$GCcontentfile
replictimingfile <- opt$replictimingfile
gamma <- opt$gamma
output_prefix <- opt$output_prefix

# Run ASCAT pipeline
ascat.prepareHTS(
    tumourseqfile = tumourseqfile,
    normalseqfile = normalseqfile,
    tumourname = tumourname,
    normalname = normalname,
    allelecounter_exe = allelecounter_exe,
    alleles.prefix = alleles.prefix,
    loci.prefix = loci.prefix,
    gender = gender,
    genomeVersion = genomeVersion,
    nthreads = nthreads,
    tumourLogR_file = tumourLogR_file,
    tumourBAF_file = tumourBAF_file,
    normalLogR_file = normalLogR_file,
    normalBAF_file = normalBAF_file
)

ascat.bc <- ascat.loadData(
    Tumor_LogR_file = tumourLogR_file,
    Tumor_BAF_file = tumourBAF_file,
    Germline_LogR_file = normalLogR_file,
    Germline_BAF_file = normalBAF_file,
    gender = gender,
    genomeVersion = genomeVersion
)

ascat.plotRawData(ascat.bc, img.prefix = paste0(output_prefix, "_Before_correction_"))

ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = GCcontentfile, replictimingfile = replictimingfile)

ascat.plotRawData(ascat.bc, img.prefix = paste0(output_prefix, "_After_correction_"))

ascat.bc <- ascat.aspcf(ascat.bc)

ascat.plotSegmentedData(ascat.bc, img.prefix = paste0(output_prefix, "_Segmented_"))

ascat.output <- ascat.runAscat(ascat.bc, gamma = gamma, write_segments = TRUE)

QC <- ascat.metrics(ascat.bc, ascat.output)

# save(ascat.bc, ascat.output, QC, file = paste0(output_prefix, "_objects.Rdata"))

# Write out CNVs in bed format
cnvs <- ascat.output$segments[2:6]
write.table(cnvs, file = paste(tumourname, ".cnvs.txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

# Write out purity and ploidy info
summary <- tryCatch(
    {
        matrix(c(ascat.output$aberrantcellfraction, ascat.output$ploidy), ncol = 2, byrow = TRUE)
    },
    error = function(err) {
        # error handler picks up where error was generated
        print(paste("Could not find optimal solution:  ", err))
        return(matrix(c(0, 0), nrow = 1, ncol = 2, byrow = TRUE))
    }
)
colnames(summary) <- c("AberrantCellFraction", "Ploidy")
write.table(summary, file = paste(tumourname, ".purityploidy.txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)


# print QC metrics to standard out.
# print(QC)

# Write QC metrics to a file
write.table(t(as.matrix(QC)), file = paste0(output_prefix, "_QC_metrics.txt"), sep = "\t", quote = F, row.names = TRUE, col.names = FALSE)
