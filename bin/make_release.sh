#!/bin/sh


VERSION=$1
ICBI=$2

if [ "$ICBI" != "" ]; then
    BASE_DIR=nextNEOpi_v${VERSION}_icbi
    BUILD_DIR=make_nextNEOpi_release/${BASE_DIR}/
else
    BASE_DIR=nextNEOpi_v${VERSION}
    BUILD_DIR=make_nextNEOpi_release/${BASE_DIR}/
fi
mkdir -p ${BUILD_DIR}/{conf,bin,assets,resources}

cp nextNEOpi.nf example_batchFile_FASTQ.csv example_batchFile_BAM.csv README.md README.html LICENSE ${BUILD_DIR}
cp assets/Alleles_list.txt ${BUILD_DIR}/assets
cp assets/Final_gbm_model.rds ${BUILD_DIR}/assets
cp assets/NeoAg_immunogenicity_predicition_GBM.patch ${BUILD_DIR}/assets
cp assets/Frameshift.pm ${BUILD_DIR}/assets
cp assets/Wildtype.pm ${BUILD_DIR}/assets
cp assets/email_template.html ${BUILD_DIR}/assets
cp assets/email_template.txt ${BUILD_DIR}/assets
cp assets/gatkPythonPackageArchive.zip ${BUILD_DIR}/assets
cp assets/hlahdenv.yml ${BUILD_DIR}/assets
cp assets/hlaii_models.txt ${BUILD_DIR}/assets
cp assets/hlaii_secondary.txt ${BUILD_DIR}/assets
cp assets/hlaii_supported.txt ${BUILD_DIR}/assets
cp assets/nextNEOpi.yml ${BUILD_DIR}/assets
cp assets/pVACseqAlleles.txt ${BUILD_DIR}/assets
cp assets/sendmail_template.txt ${BUILD_DIR}/assets

cp bin/CCF.R ${BUILD_DIR}/bin
cp bin/CSiN.py ${BUILD_DIR}/bin
cp bin/HLAHD2mixMHC2pred.py ${BUILD_DIR}/bin
cp bin/HLA_parser.py ${BUILD_DIR}/bin
cp bin/NameToID.py ${BUILD_DIR}/bin
cp bin/NeoAg_immunogenicity_predicition_GBM.R ${BUILD_DIR}/bin
cp bin/SequenzaScript.R ${BUILD_DIR}/bin
cp bin/add_CCF.py ${BUILD_DIR}/bin
cp bin/ascat.R ${BUILD_DIR}/bin
cp bin/assess_significance.R ${BUILD_DIR}/bin
cp bin/check_pe.py ${BUILD_DIR}/bin
cp bin/convertAlleleCounts.r ${BUILD_DIR}/bin
cp bin/freec2bed.pl ${BUILD_DIR}/bin
cp bin/freec2circos.pl ${BUILD_DIR}/bin
cp bin/get_epitopes.py ${BUILD_DIR}/bin
cp bin/immuno_score.py ${BUILD_DIR}/bin
cp bin/makeGraph.R ${BUILD_DIR}/bin
cp bin/make_hc_vcf.py ${BUILD_DIR}/bin
cp bin/mkCCF_input.py ${BUILD_DIR}/bin
cp bin/mutationalLoad.py ${BUILD_DIR}/bin
cp bin/parse_mixMHC2pred.py ${BUILD_DIR}/bin
cp bin/pepChopper.py ${BUILD_DIR}/bin
cp bin/run_ascat.r ${BUILD_DIR}/bin
cp bin/mkSampleInfo.py ${BUILD_DIR}/bin
cp bin/parse_blast_result.py ${BUILD_DIR}/bin
cp bin/make_peptide_fasta.py ${BUILD_DIR}/bin

if [ "$ICBI" != "" ]; then
    cp conf/params_icbi.config  ${BUILD_DIR}/conf/params.config
    # cp conf/resources_icbi.config  ${BUILD_DIR}/conf/resources.config
    cp conf/resources.config  ${BUILD_DIR}/conf/resources.config
    cp conf/process_icbi.config  ${BUILD_DIR}/conf/process.config
else
    cp conf/params.config  ${BUILD_DIR}/conf
    cp conf/resources.config  ${BUILD_DIR}/conf
    cp conf/process.config  ${BUILD_DIR}/conf
fi

cp conf/profiles.config  ${BUILD_DIR}/conf

cd $BUILD_DIR
cd ..
tar -cvzf ${BASE_DIR}.tar.gz $BASE_DIR
