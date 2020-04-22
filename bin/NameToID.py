#!/usr/bin/env python

import sys, argparse
import ntpath

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        usage="NameToID.py [-h] -i {/path/to/counts/input_file/} -a {/path/to/counts/annotation.gtf/}"
    )  # -o {/path/to/outDir}')
    parser.add_argument("-i", "--InFile", help="input file", required=True)
    parser.add_argument("-a", "--AnnotationFile", help="annotation GTF file", required=True)
    parser.add_argument("-o", "--outDir", help="Output dir", required=True)
    args = parser.parse_args()

    annoFile = args.AnnotationFile
    inFile = args.InFile
    outFile = ntpath.basename(inFile)
    outFile = args.outDir + "/" + outFile.replace("tpm.txt", "tpm_final.txt")

    gene_name = []
    gene_id = []
    attr = {}
    with open(annoFile) as in_file:
        for line in in_file:
            if line.startswith("#"):
                pass
            else:
                gene_id.append(line.split("\t")[8].strip("\n").split(" ")[1].strip('"').strip('";'))
                gene_name.append(line.split("\t")[8].strip("\n").split(" ")[7].strip('"').strip('";'))

    for i in range(0, len(gene_name)):
        attr[gene_name[i]] = gene_id[i]

    with open(inFile) as in_file2:
        header = in_file2.readline()
        with open(outFile, "w") as out_file:
            out_file.write("%s\t%s\t%s\t%s\n" % ("Gene_name", "GeneID", "Counts", "TPM"))
            for line in in_file2:
                out_file.write(
                    "%s\t%s\t%s\t%s"
                    % (
                        line.split("\t")[0],
                        attr["%s" % line.split("\t")[0]].split(".")[0],
                        line.split("\t")[1],
                        line.split("\t")[2],
                    )
                )
