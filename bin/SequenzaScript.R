library("sequenza")


args = commandArgs(trailingOnly=TRUE)

seqzFile <- args[1]
subjectID <- args[2]
gender <- args[3]


if (gender == "XY") {
    chromosomes = paste0("chr", c(1:22,"X","Y"));
}
if (gender == "XX") {
    chromosomes = paste0("chr", c(1:22,"X"));
}

seqzData <- sequenza.extract(seqzFile,
                             chromosome.list = chromosomes);

for(i in 1:length(seqzData$chromosomes)) {
    pngOut <- paste(subjectID, "_MAF_BAF_DPR_", seqzData$chromosomes[i], ".png", sep = "")
    png(pngOut,width = 1800, height = 600, units = "px", pointsize = 12);
    chromosome.view(mut.tab = seqzData$mutations[[i]],
                    baf.windows = seqzData$BAF[[i]],
                    ratio.windows = seqzData$ratio[[i]],
                    min.N.ratio = 1,
                    segments = seqzData$segments[[i]],
                    main = seqzData$chromosomes[i])
    dev.off()
}

CP.subjectID <- sequenza.fit(seqzData)

sequenza.results(sequenza.extract = seqzData,
                 cp.table = CP.subjectID,
                 sample.id = subjectID,
                 out.dir = "./");
