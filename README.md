
# Code repository for "ShaPRS: Leveraging shared genetic effects across traits or ancestries improves accuracy of polygenic scores"

DOI: XXX

This respository represents the last snapshot of the bash and R scripts used to generate our results and is provided as-is. As the analysis involved a lot of input/output operations on very large files, these were generated asynchronously on a cluster. Thus these scripts are meant to be executed on the command line manually, block-by-block, waiting for the remote jobs to finish and verifying the integrity of the resulting files at each step. To avoid code duplication, functions are defined only once across all files, however, they may be called from different scripts.

Code for most of the analyses in the paper can be found under /main/: 

**The simulation analyses used:**
1. shaPRS_ukbb_extra.sh: code relevant for the simulation analyses

**The inflammatory bowel disease analyses used:**
1. shaPRS_ibd_assoc.sh: commands for obtaining the initial GWAS results
3. shaPRS_ibd_final.sh: the final IBD analyses

**The cross-ancestry analyses used:**
1. shaPRS_crossAncestry_v4.sh: initial cross-ancestry analyses
2. shaPRS_crossAncestry_prscsx.sh: the addition of PRS-CSx as a comparison

Auxilliary analyses for Fig 4 are found in figure4.R and the code used to evaluate the PRS in the UKBB for CAD, T2D and BRCA can be found under /ukbb_validation/.
