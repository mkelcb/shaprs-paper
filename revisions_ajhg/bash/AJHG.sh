####################################################
#
# shaPRS revisions for AJHG

# IBD_reNewWave_AJHG.sh
# screen -r -D 23680.pts-566.farm5-head2

# EUR-EAS
# screen -r -D 17501.pts-627.farm5-head2

# shaPRS_revisionAFR_AJHG.sh
# screen -r -D 36827.pts-409.farm5-head2

# sims
# screen -r -D 42511.pts-58.farm5-head2


#########################
# static vars





#########################
# 3. Formal evaluation of shaPRS' performance over other methods
######


# IBD results: IBD_reNewWave_AJHG.sh ( use PRS-CS instead of LDpred2, IE NOT PRS-CSx!!)
# shaPRS
# correlation_sq 0.107 (sd: 0.000399 ) /  AUC 0.706 CI low: 0.6886  / high:  0.7239
# correlation_sq 0.0646 (sd: 0.000308 ) /   AUC 0.649 CI low: 0.6332  / high:  0.6651

# MTAG
# correlation_sq  0.0963 (sd: 0.000333 ) /  AUC  AUC 0.694 CI low: 0.6757  / high:  0.7116
# correlation_sq 0.0436 (sd: 0.000291 ) /  AUC 0.622 CI low: 0.6053  / high:  0.6379

# Meta
# correlation_sq 0.095 (sd: 0.000351 ) /   AUC 0.693 CI low: 0.6756  / high:  0.7112
# correlation_sq 0.0612 (sd: 0.000325 ) /    AUC 0.646 CI low: 0.63  / high:  0.662

# Primary only
# correlation_sq 0.103 (sd: 0.000344 ) /    AUC 0.701 CI low: 0.6834  / high:  0.7187
# correlation_sq 0.052 (sd: 0.00028 ) /     AUC 0.632 CI low: 0.6155  / high:  0.6478

#SMTPred
# correlation_sq 0.1 (sd: 0.000318 ) / AUC 0.702 CI low: 0.6844  / high:  0.7196 
# correlation_sq 0.0586 (sd: 0.000323 ) / AUC 0.642 CI low: 0.626  / high:  0.6581

# PRS1= shaPRS... so LRT p (PRS1) nested means baseline is shaPRS, and if we add another method (PRS2), does other method add to shaPRS. Whereas (LRT p (PRS2), means the baseline is the other method, and we add shaPRS to that, and tells us if shaPRS adds to that?

# shaPRS vs Meta
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_shaprs_het_PRSCS/gwas3_cd_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_cd_meta_het_PRSCS/gwas3_cd_meta_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_meta_het_vsShaPRS' '0' '1' #  r2redux diff p:  0.00047 ,  delong p-val 0.000699   / LRT p (PRS1):  0.49  / LRT p (PRS2):  5.04e-14
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_shaprs_het_PRSCS/gwas3_uc_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_uc_meta_het_PRSCS/gwas3_uc_meta_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_meta_het_vsShaPRS' '0' '1' #  r2redux diff p:  0.222,   delong p-val 0.368   / LRT p (PRS1):  0.019  / LRT p (PRS2):  1.67e-06 

# shaPRS vs Primary only
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_shaprs_het_PRSCS/gwas3_cd_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_cd_primary_het_PRSCS/gwas3_cd_primary_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_primary_het_vsShaPRS' '0' '1' # r2redux diff p:  0.341,   p-val 0.301  / LRT p (PRS1):  1.63e-06  / LRT p (PRS2):  2.72e-11
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_shaprs_het_PRSCS/gwas3_uc_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_uc_primary_het_PRSCS/gwas3_uc_primary_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_primary_het_vsShaPRS' '0' '1' # r2redux diff p:  0.000804 , p-val 0.00043 / LRT p (PRS1):  0.0919  / LRT p (PRS2):  7.17e-16 

# shaPRS vs SMTPred
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_shaprs_het_PRSCS/gwas3_cd_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_cd_het_SMTPRED/_SMTPRED.profile' $newResultsLoc$'gwas3_smtpred_cd_het_vsShaPRS' '0' '1' # r2redux diff p:  0.0409,   delong p-val 0.23  / LRT p (PRS1):  0.055  / LRT p (PRS2):  2.29e-09
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_shaprs_het_PRSCS/gwas3_uc_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_uc_het_SMTPRED/_SMTPRED.profile' $newResultsLoc$'gwas3_smtpred_uc_het_vsShaPRS' '0' '1' #  r2redux diff p:  0.0293,   delong p-val 0.0439 / LRT p (PRS1):  0.171  / LRT p (PRS2):  1.53e-08 

# shaPRS vs MTAG
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_shaprs_het_PRSCS/gwas3_cd_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_cd_mtag_het_PRSCS/gwas3_cd_mtag_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_mtag_het_vsShaPRS' '0' '1' #  r2redux diff p:  0.053,  delong p-val 0.0412   / LRT p (PRS1):  4.55e-07  / LRT p (PRS2):  3.69e-18
Compare_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_shaprs_het_PRSCS/gwas3_uc_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_uc_mtag_het_PRSCS/gwas3_uc_mtag_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_mtag_het_vsShaPRS' '0' '1' #  r2redux diff p:  2.44e-06,  delong p-val 5.31e-06  / LRT p (PRS1):  0.138  / LRT p (PRS2):  1.2e-24



# interesting, adding shaPRS improves the baseline more frequently than the other way around, suggesting that it improves things more


function Compare_PRS { 
allScoreLoc=$1
allScoreLoc2=$2
outResultLoc=$3
doAppend=$4
doAUC=$5

# we are working with PLINK2 profile scores, which are 6 columns, but the R scripts expects a 2 col one now
awk '{print $2"\t"$3"\t"$6 }' $allScoreLoc > $allScoreLoc$'_3col'
awk '{print $2"\t"$3"\t"$6 }' $allScoreLoc2 > $allScoreLoc2$'_3col'

arguments='/nfs/users/nfs_m/mk23/scripts/PRSComparator.R '$allScoreLoc$'_3col '$allScoreLoc2$'_3col '$outResultLoc$'_r2 '$doAUC$' 1 '$doAppend$' 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}
export -f Compare_PRS # this makes local functions executable when bsubbed



#########################
#########################
# EAS-EUR results: shaPRS_crossAncestry_v4.sh  shaPRS_revisionAFR_AJHG.sh and  shaPRS_crossAncestry_prscsx_AJHG.sh
# this had to change, when I cleanued up my Sanger storage

crossAncestryResults='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/results/'
crossAncestrySumStats='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/sumstats/'

PRS_array=( 'EUR_JAP_asthma' 'EUR_JAP_height' 'EUR_JAP_BRCA' 'EUR_JAP_CAD' 'EUR_JAP_T2D' )
EUR_array=( 'EUR_asthma_hm3' 'EUR_height_hm3' 'EUR_BRCA_hm3' 'EUR_CAD_hm3' 'EUR_T2D_hm3' )
JPT_array=( 'JP_asthma_hm3' 'JP_height_hm3' 'JP_BRCA_hm3' 'JP_CAD_hm3' 'JP_T2D_hm3' )
EUR_CUSTOMLDREFARRAY=( $asthmaEURRef $heightEURRef $brcaEURRef $cadEURRef $t2dEURRef )
JAP_CUSTOMLDREFARRAY=( $asthmaEASRef $heightEASRef $brcaEASRef $cadEASRef $t2dEASRef )
PRS_binary=( '1' '0' '1' '1' '1' )
famFiles=( $crossAncestrySumStats$'all_hm3.fam' $crossAncestrySumStats$'all_height_hm3.fam' )



PRSName=$pheno$'_EUR_PRSCSx.PRS'
calcAUC3=$binPheno
subsetFile=$crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'


# evaluates ready made profile scores
function evaluatePRS_r2redux {
famFile=$1
PRSName=$2
calcAUC3=$3
subsetFile=$4

# find out accuracy 
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $crossAncestryResults$PRSName$'.profile' FS=" " $famFile  > $outste$crossAncestryResults$PRSName$'.csv'

# need to subset the final PRS file to be the same as the PRS-CSx ones ( must be in the same order as the subset file, so that they are aligned
awk 'FNR == NR { file1[ $2 ] = $0; next; } FNR <= NR { { if ( $2 in file1 || FNR == 1 ) { print file1[ $2 ]  } } }' OFS='\t' $crossAncestryResults$PRSName$'.profile' $subsetFile   > $crossAncestryResults$PRSName$'_PRSCSxsubset.sscore'

#head $crossAncestryResults$PRSName$'_PRSCSxsubset.sscore'
#wc -l $crossAncestryResults$PRSName$'_PRSCSxsubset.sscore'

#head $subsetFile
#wc -l $subsetFile
arguments='/nfs/users/nfs_m/mk23/scripts/correlator_r2redux.R '$crossAncestryResults$PRSName$'.csv '$crossAncestryResults$PRSName$'_res_r2redux '$calcAUC3
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}



j=2
# PRS-CSx
arraylength=${#PRS_array[@]}
for (( j=1; j<${arraylength}+1; j++ )); do

#### Only for the first 2 traits that I have at Sanger: (the others are sent to Elena)

if [ "$j" -lt 3 ]; then
famFile=${famFiles[$j-1]} 
pheno=${PRS_array[$j-1]} 
binPheno=${PRS_binary[$j-1]} 

# PRS-CSx
# 3) Perform Stage 2: find the best linear combination of the EUR/EAS + evaluate
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb_r2redux.R '$crossAncestryResults$pheno$'_EUR_PRSCSx.PRS.profile '$crossAncestryResults$pheno$'_EAS_PRSCSx.PRS.profile '$famFile$' '$crossAncestrySumStats$'all_TEST '$binPheno$' '$crossAncestryResults$pheno$'_PRSCSx_res_r2redux'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
# this outputs $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'

# asthma: correlation_sq 0.0134 (CI_LOW: 0.0121 CI_HIGH: 0.0147 ) AUC 0.603 CI low: 0.5981  / high:  0.6079
# height:  correlation_sq 0.123 (CI_LOW: 0.12 CI_HIGH: 0.127 )

# PRS-CSx-unweighted    #correlation_sq 0.121 (sd: 5.48e-06 )
evaluatePRS_r2redux $famFile $pheno$'_EUR_PRSCSx.PRS' $binPheno $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'
# asthma:  correlation_sq 0.0113 (CI_LOW: 0.0105 CI_HIGH: 0.0122 )  AUC 0.595 CI low: 0.5916  / high:  0.5985
# height:  correlation_sq 0.121 (CI_LOW: 0.118 CI_HIGH: 0.123 )


# a) PRS-CS
evaluatePRS_r2redux $famFile $pheno$'_EUR_EUR_PRSCS.PRS' $binPheno $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'
# asthma: correlation_sq 0.00986 (CI_LOW: 0.0091 CI_HIGH: 0.0106 )   AUC 0.588 CI low: 0.5849  / high:  0.5918
# height: correlation_sq 0.116 (CI_LOW: 0.114 CI_HIGH: 0.118 )


# ShaPRS + PRS-CS
evaluatePRS_r2redux $famFile $pheno$shaPRS$'_EUR_PRSCS.PRS' $binPheno $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'
# asthma:  correlation_sq 0.0123 (CI_LOW: 0.0114 CI_HIGH: 0.0132 )   AUC 0.599 CI low: 0.5951  / high:  0.602 
# height:  correlation_sq 0.122 (CI_LOW: 0.119 CI_HIGH: 0.124 ) 


# b) LDpred2
# evaluate the LDpred2 PRS too
evaluatePRS_r2redux $famFile $pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS' $binPheno $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'
# asthma:  correlation_sq 0.0115 (CI_LOW: 0.0107 CI_HIGH: 0.0124 )    AUC 0.595 CI low: 0.5919  / high:  0.5989 
# height: correlation_sq 0.0976 (CI_LOW: 0.0954 CI_HIGH: 0.0998 )


# shaPRS+LDpred2
evaluatePRS_r2redux $famFile $pheno$'_shaPRS_LDpred2_CSxsubset.PRS' $binPheno $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore'
# asthma:  correlation_sq 0.0136 (CI_LOW: 0.0128 CI_HIGH: 0.0146 )    AUC 0.604 CI low: 0.6004  / high:  0.6073
# height:  correlation_sq 0.122 (CI_LOW: 0.119 CI_HIGH: 0.124 )

# Investigate why Height had 'same' results for both shaPRS LDpred2/PRS-CS?:
# head /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/results/EUR_JAP_height_shaPRS_EUR_PRSCS.PRS_res_r2redux
#correlation     correlation^2   corr_p  r2redux_p       CI_LOW  CI_HIGH
#0.348682628012399       0.121579575077633       0       0       0.119190069286422       0.123986577115622

# head /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/results/EUR_JAP_height_shaPRS_LDpred2_CSxsubset.PRS_res_r2redux
# correlation     correlation^2   corr_p  r2redux_p       CI_LOW  CI_HIGH
#0.348986381403989       0.12179149440545        0       0       0.119400474248912       0.12420000236731
# they are NOT exactly the same, just very similar!

###################################
# Now compare PRS methods against shaPRS:

#$crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore' # PRS-CSx
#$crossAncestryResults$pheno$'_EUR_PRSCSx.PRS_PRSCSxsubset.sscore' # PRS-CSx-stage1
#$crossAncestryResults$pheno$'_EUR_EUR_PRSCS.PRS_PRSCSxsubset.sscore' # PRS-CS
#$crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' # ldpred2 baseline#

#$crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS_PRSCSxsubset.sscore' # shaPRS + PRS-CS
#$crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' # shaPRS + LDpred2



# PRS-CSx 
#vs shaPRS+PRS-CS
Compare_PRS $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore' $crossAncestryResults$pheno$'PRS_CSx_vsShaPRS_PRSCS' '0' $binPheno 
# asthma:  r2redux diff p:  9.42e-07, delong p-val 3.29e-06  / LRT p (PRS1):  2.11e-37  / LRT p (PRS2):  0.00296
# height: r2redux diff p:  0.0011  / LRT p (PRS1):  2.66e-77  / LRT p (PRS2):  7.05e-34
# vs shaPRS+LDpred2
Compare_PRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_PRSCSx_res_r2redux.sscore' $crossAncestryResults$pheno$'PRS_CSx_vsShaPRS_ldpred2' '0' $binPheno 
#  asthma: r2redux diff p:  0.047, delong p-val 0.0947  / LRT p (PRS1):  2.02e-16  / LRT p (PRS2):  9.32e-34 
# height:  r2redux diff p:  0.821  / LRT p (PRS1):  0  / LRT p (PRS2):  0

# PRS-CSx-stage1
#vs shaPRS+PRS-CS
Compare_PRS $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_EUR_PRSCSx.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'PRSCSx_stage1_vsShaPRS_PRSCS' '0' $binPheno 
# asthma:  r2redux diff p:  0.0026, delong p-val 0.00257  / LRT p (PRS1):  7.42e-09  / LRT p (PRS2):  7.12e-32
# height: r2redux diff p:  0.387  / LRT p (PRS1):  1.87e-56  / LRT p (PRS2):  6.5e-69
# vs shaPRS+LDpred2
Compare_PRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_EUR_PRSCSx.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'PRSCSx_stage1_vsShaPRS_ldpred2' '0' $binPheno 
#  asthma: r2redux diff p:  4.83e-13, delong p-val 5.64e-12  / LRT p (PRS1):  3.26e-07  / LRT p (PRS2):  3.28e-81
# height:  r2redux diff p:  0.127  / LRT p (PRS1):  2.6e-297  / LRT p (PRS2):  0

# PRS-CS
#vs shaPRS+PRS-CS
Compare_PRS $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_EUR_EUR_PRSCS.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'PRSCS_vsShaPRS_PRSCS' '0' $binPheno 
# asthma:  r2redux diff p:  8.78e-16, delong p-val 1.58e-15  / LRT p (PRS1):  0.0357  / LRT p (PRS2):  5.72e-70
# height: r2redux diff p:  1.13e-20  / LRT p (PRS1):  1.34e-28  / LRT p (PRS2):  5.38e-191
# vs shaPRS+LDpred2
Compare_PRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_EUR_EUR_PRSCS.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'PRSCS_vsShaPRS_ldpred2' '0' $binPheno 
#  asthma:  r2redux diff p:  3.84e-27,  delong p-val 1.91e-25  / LRT p (PRS1):  0.0027  / LRT p (PRS2):  8.29e-122
# height:  r2redux diff p:  1.33e-09  / LRT p (PRS1):  4.26e-251  / LRT p (PRS2):  0

# LDpred2
#vs shaPRS+PRS-CS
Compare_PRS $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'ldpred2_vsShaPRS_PRSCS' '0' $binPheno 
# asthma:  r2redux diff p:  0.428, delong p-val 0.384  / LRT p (PRS1):  5.98e-39  / LRT p (PRS2):  9.28e-47
# height: r2redux diff p:  4.35e-81  / LRT p (PRS1):  1.46e-196  / LRT p (PRS2):  0
# vs shaPRS+LDpred2
Compare_PRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $crossAncestryResults$pheno$'ldpred2_vsShaPRS_ldpred2' '0' $binPheno 
#  asthma:  r2redux diff p:  3.96e-11,  delong p-val 6.07e-10  / LRT p (PRS1):  2.22e-05  / LRT p (PRS2):  1.89e-64 
# height:   r2redux diff p:  1.21e-244  / LRT p (PRS1):  1.37e-08  / LRT p (PRS2):  0

fi



$crossAncestryResults$PRSName$'.profile'
$crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS.profile'


#############
# EUR-AFR results : shaPRS_revisionAFR_AJHG.sh

# these paths changed due to having had to move my stuff at Sanger
asthmaCrossAncestry='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/crossAncestry/'
crossAncestrySumStats=$asthmaCrossAncestry$'sumstats/'
crossAncestryResults=$asthmaCrossAncestry$'results/'
crossAncestryRaw=$asthmaCrossAncestry$'raw/'
prscssorig=$crossAncestryRaw$'prscss/'

shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"
prscsrefs=$crossAncestryRaw$'prscsrefs/'
prscsscript=$crossAncestryRaw$'prscsscript/'
crossAncestryScratch=$asthmaCrossAncestry$'scratch/'



# evaluates ready made profile scores
function evaluatePRS_AFR_r2redux {
phenoFile=$1
PRSName=$2
calcAUC3=$3
subsetFile=$4

# find out accuracy 
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $revisionsResults$PRSName$'.profile' FS="\t" $phenoFile  > $revisionsResults$PRSName$'.csv'

# need to subset the final PRS file to be the same as the PRS-CSx ones ( must be in the same order as the subset file, so that they are aligned
awk 'FNR == NR { file1[ $2 ] = $0; next; } FNR <= NR { { if ( $2 in file1 || FNR == 1 ) { print file1[ $2 ]  } } }' OFS='\t' $revisionsResults$PRSName$'.profile' $subsetFile   > $revisionsResults$PRSName$'_PRSCSxsubset.sscore'

#head $crossAncestryResults$PRSName$'_PRSCSxsubset.sscore'
#wc -l $crossAncestryResults$PRSName$'_PRSCSxsubset.sscore'

#head $subsetFile
#wc -l $subsetFile
arguments='/nfs/users/nfs_m/mk23/scripts/correlator_r2redux.R '$revisionsResults$PRSName$'.csv '$revisionsResults$PRSName$'_res_r2redux '$calcAUC3
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}








PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' 'AFR_EUR_LDL' ) 
EUR_array=( 'bmi_EUR_HM3' 'height_EUR_HM3' 'LDL_EUR_HM3' )
AFR_array=( 'bmi_AFR_HM3' 'height_AFR_HM3' 'LDL_AFR_HM3' )

phenoFiles=( $revisionsRaw$'UKBB_BMI_AFR.phe' $revisionsRaw$'UKBB_height_AFR.phe' $revisionsRaw$'UKBB_LDL_AFR.phe' )



arraylength=${#PRS_array[@]}
# PRS-CS: baseline PRS-CS, with the built in AFR LDpanel
# ShaPRS-PRS-CS: PRS-CS  (IE not CSx) applied onto the post shaPRS data

j=2
for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno='0'
AFRPheno=${AFR_array[$j-1]} 
phenoFile=${phenoFiles[$j-1]} 


# PRS-CSx
# 3) Perform Stage 2: find the best linear combination of the AFR/EUR + evaluate
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb_r2redux.R '$revisionsResults$pheno$'_AFR_PRSCSx.PRS.profile '$revisionsResults$pheno$'_EUR_PRSCSx.PRS.profile '$phenoFile$' '$revisionsRaw$'all_AFR_hm3_TEST '$binPheno$' '$revisionsResults$pheno$'_PRSCSx_res_r2redux'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
# this outputs $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore'
# BMI:  correlation_sq 0.021 (CI_LOW: 0.0119 CI_HIGH: 0.0317 )
# Height:  correlation_sq 0.0282 (CI_LOW: 0.0176 CI_HIGH: 0.0404 )
# LDL:  correlation_sq 0.0563 (CI_LOW: 0.0411 CI_HIGH: 0.0732 )


# PRS-CSx-unweighted
evaluatePRS_AFR_r2redux $phenoFile $pheno$'_AFR_PRSCSx.PRS' $binPheno $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore'
# AFR_EUR_BMI correlation_sq 0.0053 (CI_LOW: 0.00217 CI_HIGH: 0.00931 )
# AFR_EUR_height  correlation_sq 0.00911 (CI_LOW: 0.00489 CI_HIGH: 0.0142 )
# AFR_EUR_LDL correlation_sq 0.0296 (CI_LOW: 0.0216 CI_HIGH: 0.0385 )


#PRS-CS
evaluatePRS_AFR_r2redux $phenoFile $pheno$'_AFR_AFR_PRSCS.PRS' $binPheno $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore'
# AFR_EUR_BMI correlation_sq 0.0037 (CI_LOW: 0.00116 CI_HIGH: 0.00714 )
# AFR_EUR_height correlation_sq 0.00387 (CI_LOW: 0.00126 CI_HIGH: 0.00737 )
# AFR_EUR_LDL  correlation_sq 0.00816 (CI_LOW: 0.00409 CI_HIGH: 0.0132 ) 

#SHAPRS PRS-CS
evaluatePRS_AFR_r2redux $phenoFile $pheno$shaPRS$'_AFR_PRSCS.PRS' $binPheno $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore'
# AFR_EUR_BMI  correlation_sq 0.0249 (CI_LOW: 0.0178 CI_HIGH: 0.033 )
# AFR_EUR_height correlation_sq 0.019 (CI_LOW: 0.0127 CI_HIGH: 0.026 )
# AFR_EUR_LDL correlation_sq 0.0187 (CI_LOW: 0.0124 CI_HIGH: 0.0259 )

# LDpred2: baseline primary via LDpred2
evaluatePRS_AFR_r2redux $phenoFile $pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS' $binPheno $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore'
# AFR_EUR_BMI  correlation_sq 0.00254 (CI_LOW: 0.000506 CI_HIGH: 0.00547 ) 
# AFR_EUR_height correlation_sq 0.00349 (CI_LOW: 0.00103 CI_HIGH: 0.00684 )
# AFR_EUR_LDL  correlation_sq 0.0559 (CI_LOW: 0.045 CI_HIGH: 0.0677 )

# ShaPRS+ LDpred2: shaPRS using the custom LD panel via LDpred2
evaluatePRS_AFR_r2redux $phenoFile $pheno$'_shaPRS_LDpred2_CSxsubset.PRS' $binPheno $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore'
# AFR_EUR_BMI  correlation_sq 0.0196 (CI_LOW: 0.0132 CI_HIGH: 0.0268 )
# AFR_EUR_height  correlation_sq 0.0126 (CI_LOW: 0.00757 CI_HIGH: 0.0185 )
# AFR_EUR_LDL  correlation_sq 0.0208 (CI_LOW: 0.0141 CI_HIGH: 0.0284 )



###################################
# Now compare PRS methods against shaPRS:

# $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore' # PRS-CSx
# $revisionsResults$pheno$'_AFR_PRSCSx.PRS_PRSCSxsubset.sscore' # PRS-CSx-stage1
# $revisionsResults$pheno$'_AFR_AFR_PRSCS.PRS_PRSCSxsubset.sscore' # PRS-CS
# $revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' # ldpred2 baseline#
# $revisionsResults$pheno$shaPRS$'_AFR_PRSCS.PRS_PRSCSxsubset.sscore' # shaPRS + PRS-CS
# $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' # shaPRS + LDpred2

# PRS-CSx 
#vs shaPRS+PRS-CS
Compare_PRS $revisionsResults$pheno$shaPRS$'_AFR_PRSCS.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore' $revisionsResults$pheno$'PRS_CSx_vsShaPRS_PRSCS' '0' $binPheno 
# BMI:    r2redux diff p:  0.159  / LRT p (PRS1):  0.00711  / LRT p (PRS2):  1.49e-07
# height: r2redux diff p:  0.00981  / LRT p (PRS1):  4.16e-09  / LRT p (PRS2):  0.269
# LDL:    r2redux diff p:  1.13e-11  / LRT p (PRS1):  1.83e-30  / LRT p (PRS2):  0.794 

# vs shaPRS+LDpred2
Compare_PRS $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_PRSCSx_res_r2redux.sscore' $revisionsResults$pheno$'PRS_CSx_vsShaPRS_ldpred2' '0' $binPheno 
# BMI:    r2redux diff p:  0.915  / LRT p (PRS1):  2.25e-06  / LRT p (PRS2):  8.98e-07
# height:   r2redux diff p:  0.000104  / LRT p (PRS1):  1.47e-12  / LRT p (PRS2):  0.996 
# LDL:  r2redux diff p:  9.81e-09  / LRT p (PRS1):  1.01e-27  / LRT p (PRS2):  0.297 


# PRS-CSx-stage1
#vs shaPRS+PRS-CS
Compare_PRS $revisionsResults$pheno$shaPRS$'_AFR_PRSCS.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_AFR_PRSCSx.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'PRSCSx_stage1_vsShaPRS_PRSCS' '0' $binPheno 
# BMI:    r2redux diff p:  0.00164  / LRT p (PRS1):  0.00129  / LRT p (PRS2):  2.62e-16 
# height: r2redux diff p:  0.25  / LRT p (PRS1):  2.32e-06  / LRT p (PRS2):  8.41e-11 
# LDL: r2redux diff p:  0.00494  / LRT p (PRS1):  4.04e-15  / LRT p (PRS2):  0.00117

# vs shaPRS+LDpred2
Compare_PRS $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_AFR_PRSCSx.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'PRSCSx_stage1_vsShaPRS_ldpred2' '0' $binPheno 
# BMI:    r2redux diff p:  0.0207  / LRT p (PRS1):  0.000881  / LRT p (PRS2):  2.33e-12
# height:  r2redux diff p:  0.847  / LRT p (PRS1):  3.81e-07  / LRT p (PRS2):  7.73e-08
# LDL: r2redux diff p:  0.0475  / LRT p (PRS1):  1.62e-13  / LRT p (PRS2):  4.11e-05

# PRS-CS
#vs shaPRS+PRS-CS
Compare_PRS $revisionsResults$pheno$shaPRS$'_AFR_PRSCS.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_AFR_AFR_PRSCS.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'PRSCS_vsShaPRS_PRSCS' '0' $binPheno 
# BMI:    r2redux diff p:  0.00164  / LRT p (PRS1):  0.000209  / LRT p (PRS2):  1.6e-17
# height:  r2redux diff p:  0.0233  / LRT p (PRS1):  0.000376  / LRT p (PRS2):  3.22e-12
# LDL:   r2redux diff p:  0.289  / LRT p (PRS1):  0.000185  / LRT p (PRS2):  6.69e-08
# vs shaPRS+LDpred2
Compare_PRS $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_AFR_AFR_PRSCS.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'PRSCS_vsShaPRS_ldpred2' '0' $binPheno 
# BMI:     r2redux diff p:  0.0173  / LRT p (PRS1):  0.000349  / LRT p (PRS2):  3.29e-13
# height:   r2redux diff p:  0.177  / LRT p (PRS1):  0.00025  / LRT p (PRS2):  1.16e-08
# LDL:  r2redux diff p:  0.0589  / LRT p (PRS1):  0.00096  / LRT p (PRS2):  3.06e-10 


# LDpred2
#vs shaPRS+PRS-CS
Compare_PRS $revisionsResults$pheno$shaPRS$'_AFR_PRSCS.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'ldpred2_vsShaPRS_PRSCS' '0' $binPheno 
# BMI:    r2redux diff p:  0.000783  / LRT p (PRS1):  0.000502  / LRT p (PRS2):  8.24e-18
# height:  r2redux diff p:  0.00638  / LRT p (PRS1):  0.00251  / LRT p (PRS2):  1.26e-12
# LDL: r2redux diff p:  6.68e-11  / LRT p (PRS1):  6.42e-43  / LRT p (PRS2):  0.000152

# vs shaPRS+LDpred2
Compare_PRS $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS_PRSCSxsubset.sscore' $revisionsResults$pheno$'ldpred2_vsShaPRS_ldpred2' '0' $binPheno 
# BMI:     r2redux diff p:  0.00879  / LRT p (PRS1):  0.0012  / LRT p (PRS2):  2.36e-13
# height: r2redux diff p:  0.0676  / LRT p (PRS1):  0.00217  / LRT p (PRS2):  5.8e-09
# LDL:  r2redux diff p:  1.08e-08  / LRT p (PRS1):  3.89e-42  / LRT p (PRS2):  8.2e-07 


done



#########################
# 2. Double check the AFR LDL results by comparing its SD against the other traits:
# investigate why LDL for LDpred2 is so good, by checking the CIs 
# double check the results for the AFR LDL LDpred2 results, particularly the CIs
######

# shaPRS_revisionAFR_AJHG.sh
 $revisionsScratch$'LDL_AFR_hm3' # 1362368
wc -l $revisionsRaw$'LDL_EUR_HM3' # 849698

j=3
pheno=${PRS_array[$j-1]} 
binPheno='0'
AFRPheno=${AFR_array[$j-1]} 

#evaluatePRS_AFR_r2redux $phenoFile $pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS' $binPheno

# LDPred2
# AFR_EUR_BMI
# --score: 638552 valid predictors loaded  correlation_sq 0.00254 (sd: 0.000216 )
# AFR_EUR_height
# --score: 680312 valid predictors loaded  correlation_sq 0.00349 (sd: 0.000213 )
# AFR_EUR_LDL
# --score: 660233 valid predictors loaded. correlation_sq 0.0559 (sd: 0.000222 )

# so the LDL r2 SD calculated by Elena's bootstrapping method, so it is comparable to the rest
# there doesn't seem to be anything out of the ordinary
# My initial feeling was that this is an LD blending failure, where the AFR-EUR LD is blended in a disadvantegous way, and as LDpred2 doesn't do this, hence it is unaffected




#######################################
# Simulation results:  shaPRS_revisionSim2_AJHG.sh

baseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/'
shaprsLoc=$baseLoc$'shaprs_ukbb/'
shaprsScratch=$shaprsLoc$'scratch/'
shaprsRaw=$shaprsLoc$'raw/'
shaprsResults=$shaprsLoc$'results/'
#shaprsSimPhe=$shaprsLoc$'simPhe/'
shaprsSimPhe='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/shaprs_ukbb/simPhe/'

pheno_splitA=( 50 40 20 ) # 1= 50/50, 2=20/80
pheno_splitB=( 50 60 80 ) # 1= 50/50, 2=20/80
h2s=( 0.25 0.5 0.75 )
numCausals=( 1000 3000 5000 ) 


#currentPhe=1
#outputDir=$outste$currentPhe
# Reusable GWAS function
function perform_lrt_evaluation { 
outste=$1
PRS_FileStem=$2
pheA_testSet=$3
shaprsSimPhe=$4
pheste=$5
resultsDir=$6

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
plinkMem=2000  # how much RAM PLINK would use
# go through each pheno replicate


rm -rf $resultsDir$'/lrt_combined'
rm -rf $resultsDir$'/lrt_primary'
rm -rf $resultsDir$'/lrt_smtpred_r2'
rm -rf $resultsDir$'/lrt_mtag_r2'

for ((currentPhe=1; currentPhe<=20; currentPhe++)); do #replicates
outputDir=$outste$currentPhe


# convert these into dummy 6 col format to preserve compatibility with 'Compare_PRS'
awk '{print "0\t"$1"\t"$2"\t0\t0\t"$3 }' FS=',' $outste$currentPhe$'_subpheno.csv' > $outste$currentPhe$'_subpheno_6col'
awk '{print "0\t"$1"\t"$2"\t0\t0\t"$3 }' FS=',' $outste$currentPhe$'_combined.csv' > $outste$currentPhe$'_combined_6col'
awk '{print "0\t"$1"\t"$2"\t0\t0\t"$3 }' FS=',' $outste$currentPhe$'_mtag.csv' > $outste$currentPhe$'_mtag_6col'
awk '{print "0\t"$1"\t"$2"\t0\t0\t"$3 }' FS=','  $outste$currentPhe$'_shaPRS_meta.csv' > $outste$currentPhe$'_shaPRS_meta_6col'
awk '{print "0\t"$1"\t"$2"\t0\t0\t"$3 }' FS=','  $outste$currentPhe$'_SMTPred.csv' > $outste$currentPhe$'_SMTPred.csv_6col'


# Combined
Compare_PRS $outste$currentPhe$'_shaPRS_meta_6col' $outste$currentPhe$'_combined_6col' $resultsDir$'/lrt_combined' '1' '0'


# Primary
Compare_PRS $outste$currentPhe$'_shaPRS_meta_6col' $outste$currentPhe$'_subpheno_6col' $resultsDir$'/lrt_primary' '1' '0'


# SMTPRED
Compare_PRS $outste$currentPhe$'_shaPRS_meta_6col' $outste$currentPhe$'_SMTPred.csv_6col' $resultsDir$'/lrt_smtpred' '1' '0'

#MTAG
Compare_PRS $outste$currentPhe$'_shaPRS_meta_6col' $outste$currentPhe$'_mtag_6col' $resultsDir$'/lrt_mtag' '1' '0'

done # end of loop num replicates

# now aggregate the results and get their median -log10p 
arguments='/nfs/users/nfs_m/mk23/scripts/medianator_lrt.R '$resultsDir$' '$resultsDir$'/lrt_combined_r2 '$resultsDir$'/lrt_primary_r2 '$resultsDir$'/lrt_smtpred_r2 '$resultsDir$'/lrt_mtag_r2'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

}

export -f perform_lrt_evaluation # this makes local functions executable when bsubbed


##############

# go through the representative
j=1
r=2
e=1
p=2 # sample size 'full' ie ''
sample_size=${samplesizes[$p-1]}
sample_size_text="_full"
causals_all=${numCausals[$j-1]}
h2=${h2s[$r-1]}
h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point
current_het=${heterogeneity[$e-1]}
echo -e "combined,subpheno,SMTPred,MTAG" > $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_r2diff_eval'
echo -e "combined,subpheno,SMTPred,MTAG" > $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_lrt_eval'
echo -e "#causals,rG,p_current,shared_corr,heterogeneity,splitA:splitB,sample_size" > $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors_eval'
# go through each phenoA/B split scenario
arraylength_split=${#pheno_splitA[@]}
for (( s=1; s<${arraylength_split}+1; s++ )); do
splitA=${pheno_splitA[$s-1]}
splitB=${pheno_splitB[$s-1]}


# go through each genetic correlation
arraylength_rG=${#rGs[@]}
for (( k=1; k<${arraylength_rG}+1; k++ )); do
rG=${rGs[$k-1]}
rGText=$(awk -v rG=$rG 'BEGIN{ print rG*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point

# based on RG, select the heterogeneity architecture
case $rG in
0.1)
  Ps=( "${rG01_ps[@]}" )
  shared_corrs=( "${rG01_shared_corrs[@]}" )
  echo "rG01"
  ;;
0.25)
  Ps=( "${rG025_ps[@]}" )
  shared_corrs=( "${rG025_shared_corrs[@]}" )
  echo "rG025"
  ;;
*)
  Ps=( "${rG05_ps[@]}" )
  shared_corrs=( "${rG05_shared_corrs[@]}" )
  echo "rG05"
  ;;
esac

# go through each heterogeneity architecture
arraylength_Ps=${#Ps[@]}
for (( b=1; b<${arraylength_Ps}+1; b++ )); do
p_current=${Ps[$b-1]}
shared_corr=${shared_corrs[$b-1]}
echo "rG: "$rG$" | p: "$p_current$" / shared_corr: "$shared_corr



# go through the 2 extra/regular
arraylength_heterogeneity=${#heterogeneity[@]}
for (( e=1; e<${arraylength_heterogeneity}+1; e++ )); do
current_het=${heterogeneity[$e-1]}


resultsDir=$shaprsResults$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
outste=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'
pheste=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text

#bsub -q normal -m "modern_hardware" -G team152 -n1 -J TEST${sample_size}_${splitA}${splitB}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $outste$'_TEST.out' -e $outste$'_TEST.err'  "perform_lrt_evaluation $outste $PRS_FileStem $pheA_testSet $shaprsSimPhe $pheste $resultsDir"
perform_lrt_evaluation $outste $PRS_FileStem $pheA_testSet $shaprsSimPhe $pheste $resultsDir


# append results to ongoing files to
cat $resultsDir$'/r2diffs' >> $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_r2diff_eval'
cat $resultsDir$'/lrts' >> $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_lrt_eval'

# signature of predictors: #causals, rg, p_current, shared_corr, current_het, splitA:splitB, sample_size
echo -e "$causals_all,$rG,$p_current,$shared_corr,$current_het,$splitA:$splitB,$sample_size_text" > $resultsDir$'/predictors_r2diff' 
cat $resultsDir$'/predictors_r2diff'  >> $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors_eval'


done # end of extra/regular

done # end of heterogeneity architecture loop

done # end of rG loop

done # end of split loop




###############
# create heatmaps of the above
j=1
r=2
causals_all=${numCausals[$j-1]}
h2=${h2s[$r-1]}

arguments='/nfs/users/nfs_m/mk23/scripts/sims_evaluator_ajhg.R '$shaprsResults$causals_all$'_h2_'$h2$'sims_eval_ukbb_ajhg '$shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_r2diff_eval '$shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors_eval 0'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


arguments='/nfs/users/nfs_m/mk23/scripts/sims_evaluator_ajhg.R '$shaprsResults$causals_all$'_h2_'$h2$'sims_eval_ukbb_ajhg_lrt '$shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_lrt_eval '$shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors_eval 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments



