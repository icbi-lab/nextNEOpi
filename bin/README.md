# Gradient boosting machine model for neoantigen immunogenicity prediction

## Introduction
This gradient boosting machine (GBM) model was described by Smith et al. CIR (2019) in the manuscript entitle "Machine-learning prediction of tumor antigen immunogenicity informs selection of therapeutic epitopes".  This model can be used to predict immunogenicity scores for MHC class I single nucleotide variant (SNV) neoantigens 8-11 amino acid residues in length.  To run the model, you will require the amino acid sequence of the variant and reference peptides, as well as the site of variation.

Current tumor antigen calling algorithms primarily rely on epitope/MHC binding affinity predictions to rank and select for potential epitope targets.  These algorithms do not predict for epitope immunogenicity using approaches modeled from tumor-specific antigen data.  In the above study, we describe peptide-intrinsic biochemical features associated with neoantigen and minor histocompatibility mismatch antigen (mHA) immunogenicity and present a machine-learning gradient boosting algorithm for predicting tumor antigen immunogenicity.  This algorithm is validated in two murine tumor models, demonstrating the capacity to inform selection of therapeutically active antigens.

# Prerequisites
You will require a working installation of R.  The original analysis was performed using R version 3.5.2.  In addition, the following packages are required for running the R code:

    1. caret #Original analysis run in v6.0-84
    2. Peptides #Original analysis run in v2.4
    3. data.table #Original analysis run in v1.12.0
    4. doParallel #Original analysis run in v1.0.14

# Included files
**1. NeoAg_immunogenicity_prediction_GBM.R:** The R script for running the GBM model.

**2. Final_gbm_model.rds:** The R data file of the GBM model.

**3. TCGA_neoAg_example.txt:** Example input file containing TCGA-derived neoantigens for immunogenicity prediction.
    
# Running the model
After ensuring above R and packages are installed, the R code entitle "NeoAg_immunogenicity_prediction_GBM.R" can be run directly.  Prior to running the R code, set the necessary input path variables within the file, including:

  **neo_tab_path**: Path to the input data containing 4 columns:
  
  	1) Sample_ID
	2) mut_peptide (neoantigen peptide sequence)
	3) Reference (respective reference epitope amino acid sequence)
	4) peptide_variant_position (numeric location of variant amino acid)
An example of this format can be found at the file entitle "TCGA_neoAg_example.txt": 

  **GBM_model_path**: Path to "Final_gbm_model.rds".
  
The output variable "TCGA_predict" will contain numerical values corresponding to the predicted immunogenicity score of each respective input peptide.
 
# Documentation
Documentation for the most recent version is available on the project website. A copy of this document is included in each version as Markdown-formatted text.
  
# Lastest version
The latest release can be found at https://github.com/vincentlaboratories/neoag

# License
The model and associated scripts are licensed for non-commercial research purposes only - see LICENSE.txt for full license

# Contact
GitHub: https://github.com/vincentlaboratories/neoag


