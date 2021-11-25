
<!-- README.md is generated from README.Rmd. Please edit that file -->

## apply pipeline for shaPRS paper

This pipeline is designed to applied pre-built PGS to selected samples
of UKBB.

## Requirements

R versions \>= 3.4.0.

### To install the R libraries required for the pipeline do:

``` r
## Within R:
install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.4.tgz", repos = NULL, type = "source" )
```

# Workflow

The
[Snakefile](https://https://github.com/mkelcb/shaprs-paper/-/blob/master/UKBB_validation/Snakefile)
details the steps performed for analysis. The config.yaml file is not
inlcuded in the repo because contains paths to data directories. For
each PGS score an initial QC is performed checking if the SNPs in the
PGS are present in UKBB. Multiallelic SNPs or SNPs below a user defined
info score are discarded (rule score\_qc).

For the SNPs that passed the QC the PGS is applied to UKBB data (rules
run\_subscores and sum\_subscores). Relevant individuals to apply the
score are selected by rule subset\_score. As PRS\_CSx requires genotype
data for the target population to estimates parameters, we split the
relevant individuals to apply the score 50:50. We train the PRS-CSx in
rule prscsx\_score and we apply all scores to the ‘test’ data in rule
subset\_test. We finally assess the performance of all scores in rule
pgs\_auc\_r2\_test.
